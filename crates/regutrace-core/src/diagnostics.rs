//! Structured input review and machine-readable diagnostics.
//!
//! Every check emits a stable `code` plus a severity, a human-readable message,
//! optional file/line context, and an optional suggestion. The report is
//! serialized to JSON so a frontend (or a RAG assistant) can display or explain
//! the result without parsing CLI text.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

use crate::fasta::Fasta;
use crate::gff::{GeneModel, GffStats};
use crate::io::{EdgeRecord, TfMotif};
use crate::promoter::Promoter;
use crate::pwm::Pwm;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Severity {
    Error,
    Warning,
    Info,
}

impl Severity {
    pub fn label(self) -> &'static str {
        match self {
            Severity::Error => "ERROR",
            Severity::Warning => "WARN",
            Severity::Info => "INFO",
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Diagnostic {
    pub severity: Severity,
    pub code: String,
    pub message: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub file: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub line: Option<usize>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub context: BTreeMap<String, String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub suggestion: Option<String>,
}

impl Diagnostic {
    pub fn new(
        severity: Severity,
        code: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            severity,
            code: code.into(),
            message: message.into(),
            file: None,
            line: None,
            context: BTreeMap::new(),
            suggestion: None,
        }
    }

    pub fn error(code: impl Into<String>, message: impl Into<String>) -> Self {
        Self::new(Severity::Error, code, message)
    }

    pub fn warning(code: impl Into<String>, message: impl Into<String>) -> Self {
        Self::new(Severity::Warning, code, message)
    }

    pub fn info(code: impl Into<String>, message: impl Into<String>) -> Self {
        Self::new(Severity::Info, code, message)
    }

    pub fn with_file(mut self, path: &Path) -> Self {
        self.file = Some(path.display().to_string());
        self
    }

    pub fn with_line(mut self, line: usize) -> Self {
        self.line = Some(line);
        self
    }

    pub fn with_context(mut self, key: &str, value: impl ToString) -> Self {
        self.context.insert(key.to_string(), value.to_string());
        self
    }

    pub fn with_suggestion(mut self, suggestion: impl Into<String>) -> Self {
        self.suggestion = Some(suggestion.into());
        self
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Summary {
    pub errors: usize,
    pub warnings: usize,
    pub infos: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct Report {
    pub tool: String,
    pub version: String,
    pub command: String,
    pub inputs: BTreeMap<String, String>,
    pub summary: Summary,
    pub diagnostics: Vec<Diagnostic>,
}

impl Report {
    pub fn new(command: &str) -> Self {
        Self {
            tool: "regutrace".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            command: command.to_string(),
            inputs: BTreeMap::new(),
            summary: Summary::default(),
            diagnostics: Vec::new(),
        }
    }

    pub fn add_input(&mut self, name: &str, path: &Path) {
        self.inputs
            .insert(name.to_string(), path.display().to_string());
    }

    pub fn push(&mut self, diagnostic: Diagnostic) {
        match diagnostic.severity {
            Severity::Error => self.summary.errors += 1,
            Severity::Warning => self.summary.warnings += 1,
            Severity::Info => self.summary.infos += 1,
        }
        self.diagnostics.push(diagnostic);
    }

    pub fn error(&mut self, code: impl Into<String>, message: impl Into<String>) {
        self.push(Diagnostic::error(code, message));
    }

    pub fn warning(&mut self, code: impl Into<String>, message: impl Into<String>) {
        self.push(Diagnostic::warning(code, message));
    }

    pub fn info(&mut self, code: impl Into<String>, message: impl Into<String>) {
        self.push(Diagnostic::info(code, message));
    }

    pub fn has_errors(&self) -> bool {
        self.summary.errors > 0
    }

    pub fn write_json(&self, path: &Path) -> Result<()> {
        let text = serde_json::to_string_pretty(self).context("serialize diagnostics report")?;
        fs::write(path, format!("{text}\n"))
            .with_context(|| format!("write diagnostics report {}", path.display()))
    }

    pub fn human_lines(&self) -> Vec<String> {
        let mut lines = vec![format!(
            "{} {}: {} error(s), {} warning(s), {} info",
            self.tool, self.version, self.summary.errors, self.summary.warnings, self.summary.infos
        )];
        for d in &self.diagnostics {
            let mut line = format!("[{}] {}: {}", d.severity.label(), d.code, d.message);
            let mut loc = d.file.clone().unwrap_or_default();
            if let Some(line_no) = d.line {
                if !loc.is_empty() {
                    loc.push(':');
                }
                loc.push_str(&line_no.to_string());
            }
            if !loc.is_empty() {
                line.push_str(&format!(" ({loc})"));
            }
            if let Some(suggestion) = &d.suggestion {
                line.push_str(&format!(" — {suggestion}"));
            }
            lines.push(line);
        }
        lines
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct FastaStats {
    pub records: usize,
    pub total_bases: u64,
    pub n_bases: u64,
    pub non_acgtn_bases: u64,
    pub min_length: u64,
    pub max_length: u64,
    pub chromosomes: Vec<String>,
}

pub fn review_fasta(fa: &Fasta, path: &Path, report: &mut Report) -> FastaStats {
    let mut total = 0u64;
    let mut n_bases = 0u64;
    let mut other = 0u64;
    let mut min_len = u64::MAX;
    let mut max_len = 0u64;
    for name in &fa.order {
        if let Some(seq) = fa.seqs.get(name) {
            let len = seq.len() as u64;
            total += len;
            min_len = min_len.min(len);
            max_len = max_len.max(len);
            for &b in seq {
                match b {
                    b'N' => n_bases += 1,
                    b'A' | b'C' | b'G' | b'T' => {}
                    _ => other += 1,
                }
            }
        }
    }
    if min_len == u64::MAX {
        min_len = 0;
    }
    let stats = FastaStats {
        records: fa.order.len(),
        total_bases: total,
        n_bases,
        non_acgtn_bases: other,
        min_length: min_len,
        max_length: max_len,
        chromosomes: fa.order.clone(),
    };
    report.push(
        Diagnostic::info(
            "FASTA_SUMMARY",
            format!(
                "{} records, {} bases; lengths {}..{}",
                stats.records, stats.total_bases, stats.min_length, stats.max_length
            ),
        )
        .with_file(path)
        .with_context("records", stats.records)
        .with_context("total_bases", stats.total_bases),
    );
    if total > 0 && n_bases * 20 > total {
        let pct = n_bases as f64 * 100.0 / total as f64;
        report.push(
            Diagnostic::warning(
                "FASTA_HIGH_N",
                format!("{pct:.1}% of bases are N; assembly gaps can hide binding sites"),
            )
            .with_file(path)
            .with_context("n_bases", n_bases)
            .with_context("total_bases", total)
            .with_suggestion("check assembly gaps; N-containing windows are skipped by the scanner"),
        );
    }
    if other > 0 {
        report.push(
            Diagnostic::warning(
                "FASTA_AMBIGUOUS_BASES",
                format!("{other} bases are not A/C/G/T/N (ambiguity codes)"),
            )
            .with_file(path)
            .with_suggestion("only A/C/G/T are scored; ambiguity codes are treated as non-matching"),
        );
    }
    report.push(
        Diagnostic::info(
            "FASTA_CHROMOSOMES",
            format!("{} sequence(s)", stats.chromosomes.len()),
        )
        .with_file(path)
        .with_context("names", stats.chromosomes.join(",")),
    );
    stats
}

pub fn review_gff(genes: &[GeneModel], stats: &GffStats, path: &Path, report: &mut Report) {
    if stats.skipped_records > 0 {
        report.push(
            Diagnostic::warning(
                "GFF_SKIPPED_RECORDS",
                format!(
                    "{} feature record(s) lacked a required ID/Parent and were skipped",
                    stats.skipped_records
                ),
            )
            .with_file(path)
            .with_context("count", stats.skipped_records)
            .with_suggestion("check GFF ID/Parent attributes; skipped records do not contribute models"),
        );
    }
    if genes.is_empty() {
        report.push(
            Diagnostic::error("GFF_NO_GENES", "no gene/transcript models were found")
                .with_file(path)
                .with_suggestion("check the GFF feature types and ID/Parent attributes"),
        );
        return;
    }
    let features = stats
        .feature_counts
        .iter()
        .map(|(k, v)| format!("{k}={v}"))
        .collect::<Vec<_>>()
        .join(", ");
    report.push(
        Diagnostic::info(
            "GFF_SUMMARY",
            format!(
                "{} genes, {} transcripts; features: {features}",
                stats.genes, stats.transcripts
            ),
        )
        .with_file(path)
        .with_context("genes", stats.genes)
        .with_context("transcripts", stats.transcripts),
    );
    if stats.inferred_genes > 0 || stats.inferred_transcripts > 0 {
        report.push(
            Diagnostic::info(
                "GFF_INFERRED_MODELS",
                format!(
                    "inferred {} gene(s) and {} transcript(s) from CDS/exon records",
                    stats.inferred_genes, stats.inferred_transcripts
                ),
            )
            .with_file(path)
            .with_suggestion(
                "a GFF with explicit gene/mRNA records is preferred; the inferred models use CDS/exon spans",
            ),
        );
    }
    if stats.genes_without_transcript > 0 {
        report.push(
            Diagnostic::warning(
                "GFF_GENES_WITHOUT_TRANSCRIPT",
                format!(
                    "{} gene(s) have no mRNA/transcript record; TSS uses gene bounds",
                    stats.genes_without_transcript
                ),
            )
            .with_file(path)
            .with_context("count", stats.genes_without_transcript)
            .with_suggestion("check whether the GFF uses CDS/exon features without mRNA records"),
        );
    }
    if stats.transcripts_without_features > 0 {
        report.push(
            Diagnostic::warning(
                "GFF_TRANSCRIPTS_WITHOUT_FEATURES",
                format!(
                    "{} transcript(s) have no CDS/exon/UTR record; representative length falls back to genomic span",
                    stats.transcripts_without_features
                ),
            )
            .with_file(path)
            .with_context("count", stats.transcripts_without_features),
        );
    }
    if stats.multi_isoform_genes > 0 {
        report.push(
            Diagnostic::info(
                "GFF_MULTI_ISOFORM",
                format!(
                    "{} gene(s) have multiple isoforms; one representative transcript is selected",
                    stats.multi_isoform_genes
                ),
            )
            .with_file(path)
            .with_context("count", stats.multi_isoform_genes)
            .with_suggestion("representative transcript: longest exon/CDS+UTR length, ties prefer the upstream TSS"),
        );
    }
    let unusual: Vec<String> = stats
        .strands
        .iter()
        .filter(|(strand, _)| strand.as_str() != "+" && strand.as_str() != "-")
        .map(|(strand, count)| format!("{strand}={count}"))
        .collect();
    if !unusual.is_empty() {
        report.push(
            Diagnostic::warning(
                "GFF_UNUSUAL_STRAND",
                format!("unusual strand values: {}", unusual.join(", ")),
            )
            .with_file(path)
            .with_suggestion("strand '.' is treated as '+' unless a transcript supplies a strand"),
        );
    }
    report.push(
        Diagnostic::info(
            "GFF_CHROMOSOMES",
            format!("{} chromosome name(s)", stats.chromosomes.len()),
        )
        .with_file(path)
        .with_context("names", stats.chromosomes.join(",")),
    );
}

#[derive(Clone, Debug, Serialize)]
pub struct ChromosomeStats {
    pub total: usize,
    pub exact: usize,
    pub alias: usize,
    pub missing: Vec<String>,
}

pub fn review_chromosomes(
    genes: &[GeneModel],
    fa: &Fasta,
    gff_path: &Path,
    fasta_path: &Path,
    report: &mut Report,
) -> ChromosomeStats {
    let names: BTreeSet<String> = genes.iter().map(|g| g.chrom.clone()).collect();
    let mut exact = 0usize;
    let mut alias = 0usize;
    let mut missing = Vec::new();
    let mut alias_examples = Vec::new();
    for name in &names {
        if fa.seqs.contains_key(name) {
            exact += 1;
        } else if fa.seq(name).is_some() {
            alias += 1;
            if alias_examples.len() < 10 {
                alias_examples.push(name.clone());
            }
        } else {
            missing.push(name.clone());
        }
    }
    let stats = ChromosomeStats {
        total: names.len(),
        exact,
        alias,
        missing: missing.clone(),
    };
    if !missing.is_empty() && missing.len() == names.len() {
        report.push(
            Diagnostic::error(
                "CHROMOSOME_NAMES_NO_MATCH",
                format!(
                    "none of the {} GFF chromosome names match the genome FASTA",
                    names.len()
                ),
            )
            .with_file(gff_path)
            .with_context("gff_chromosomes", names.iter().cloned().collect::<Vec<_>>().join(","))
            .with_context("fasta_chromosomes", fa.order.join(","))
            .with_suggestion("check chromosome naming (chr1 vs 1, scaffold/contig prefixes, zero padding)"),
        );
    } else if !missing.is_empty() {
        report.push(
            Diagnostic::warning(
                "CHROMOSOME_NAMES_PARTIAL",
                format!("{} chromosome name(s) are missing from the genome FASTA", missing.len()),
            )
            .with_file(fasta_path)
            .with_context("missing", missing.iter().take(20).cloned().collect::<Vec<_>>().join(","))
            .with_suggestion("genes on missing chromosomes will be skipped"),
        );
    }
    if alias > 0 {
        report.push(
            Diagnostic::info(
                "CHROMOSOME_ALIASES_RESOLVED",
                format!(
                    "{} chromosome name(s) matched through the naming fallback",
                    alias
                ),
            )
            .with_file(gff_path)
            .with_context("examples", alias_examples.join(",")),
        );
    }
    stats
}

#[derive(Clone, Debug, Serialize)]
pub struct MotifStats {
    pub motifs: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub duplicate_ids: Vec<String>,
    pub unnamed: usize,
}

pub fn review_motifs(motifs: &[Pwm], path: &Path, report: &mut Report) -> MotifStats {
    let mut min_len = usize::MAX;
    let mut max_len = 0usize;
    let mut seen: HashMap<&str, usize> = HashMap::new();
    let mut unnamed = 0usize;
    for m in motifs {
        min_len = min_len.min(m.matrix.len());
        max_len = max_len.max(m.matrix.len());
        *seen.entry(m.id.as_str()).or_insert(0) += 1;
        if m.name.trim().is_empty() {
            unnamed += 1;
        }
    }
    if min_len == usize::MAX {
        min_len = 0;
    }
    let mut duplicate_ids: Vec<String> = seen
        .iter()
        .filter(|(_, n)| **n > 1)
        .map(|(id, n)| format!("{id}({n})"))
        .collect();
    duplicate_ids.sort();
    let stats = MotifStats {
        motifs: motifs.len(),
        min_length: min_len,
        max_length: max_len,
        duplicate_ids: duplicate_ids.clone(),
        unnamed,
    };
    report.push(
        Diagnostic::info(
            "MOTIF_SUMMARY",
            format!(
                "{} motifs; lengths {}..{}",
                stats.motifs, stats.min_length, stats.max_length
            ),
        )
        .with_file(path)
        .with_context("motifs", stats.motifs),
    );
    if !duplicate_ids.is_empty() {
        report.push(
            Diagnostic::warning(
                "MOTIF_DUPLICATE_IDS",
                format!("duplicate motif IDs: {}", duplicate_ids.join(",")),
            )
            .with_file(path)
            .with_suggestion("only the first matrix for each ID is used by scan"),
        );
    }
    if unnamed > 0 {
        report.push(
            Diagnostic::info(
                "MOTIF_UNNAMED",
                format!("{unnamed} motif(s) have no header name; auto tf-map uses the motif ID"),
            )
            .with_file(path),
        );
    }
    if min_len < 4 {
        report.push(
            Diagnostic::warning(
                "MOTIF_SHORT",
                format!("shortest motif is {min_len} bp; specificity is low"),
            )
            .with_file(path)
            .with_suggestion("short motifs produce many false-positive hits"),
        );
    }
    if max_len > 30 {
        report.push(
            Diagnostic::info(
                "MOTIF_LONG",
                format!("longest motif is {max_len} bp"),
            )
            .with_file(path),
        );
    }
    stats
}

#[derive(Clone, Debug, Serialize)]
pub struct TfMapStats {
    pub rows: usize,
    pub missing_motifs: Vec<String>,
    pub duplicate_rows: usize,
    pub empty_symbols: usize,
}

pub fn review_tf_map(
    tf_map: &[TfMotif],
    motifs: &[Pwm],
    path: &Path,
    report: &mut Report,
) -> TfMapStats {
    let available: HashSet<&str> = motifs.iter().map(|m| m.id.as_str()).collect();
    let mut missing: BTreeSet<String> = BTreeSet::new();
    let mut seen: HashSet<(&str, &str)> = HashSet::new();
    let mut duplicate_rows = 0usize;
    let mut empty_symbols = 0usize;
    for row in tf_map {
        if !available.contains(row.motif_id.as_str()) {
            missing.insert(row.motif_id.clone());
        }
        if !seen.insert((row.tf_id.as_str(), row.motif_id.as_str())) {
            duplicate_rows += 1;
        }
        if row.tf_symbol.trim().is_empty() {
            empty_symbols += 1;
        }
    }
    let missing_motifs: Vec<String> = missing.into_iter().collect();
    let stats = TfMapStats {
        rows: tf_map.len(),
        missing_motifs: missing_motifs.clone(),
        duplicate_rows,
        empty_symbols,
    };
    report.push(
        Diagnostic::info(
            "TFMAP_SUMMARY",
            format!(
                "{} TF–motif row(s); {} duplicate row(s), {} empty symbol(s)",
                stats.rows, stats.duplicate_rows, stats.empty_symbols
            ),
        )
        .with_file(path)
        .with_context("rows", stats.rows),
    );
    if !missing_motifs.is_empty() {
        let all_missing = missing_motifs.len() == tf_map.iter().map(|r| r.motif_id.as_str()).collect::<HashSet<_>>().len();
        let message = format!(
            "{} motif ID(s) are not in the motif library: {}",
            missing_motifs.len(),
            missing_motifs.iter().take(20).cloned().collect::<Vec<_>>().join(",")
        );
        if all_missing {
            report.push(
                Diagnostic::error("TFMAP_NO_MATCHING_MOTIFS", message)
                    .with_file(path)
                    .with_suggestion("check that tf-map motif_id values match the motif headers"),
            );
        } else {
            report.push(
                Diagnostic::warning("TFMAP_MISSING_MOTIFS", message)
                    .with_file(path)
                    .with_suggestion("rows with missing motifs are skipped"),
            );
        }
    }
    if duplicate_rows > 0 {
        report.push(
            Diagnostic::warning(
                "TFMAP_DUPLICATE_ROWS",
                format!("{duplicate_rows} duplicate (tf_id, motif_id) row(s)"),
            )
            .with_file(path)
            .with_suggestion("duplicates produce duplicate evidence rows; the store deduplicates on import"),
        );
    }
    if empty_symbols > 0 {
        report.push(
            Diagnostic::info(
                "TFMAP_EMPTY_SYMBOLS",
                format!("{empty_symbols} row(s) have an empty tf_symbol; tf_id is used"),
            )
            .with_file(path),
        );
    }
    stats
}

#[derive(Clone, Debug, Serialize)]
pub struct AnnotationStats {
    pub entries: usize,
    pub genes_total: usize,
    pub genes_with_symbol: usize,
}

pub fn review_annotation(
    annotation: &HashMap<String, String>,
    genes: &[GeneModel],
    path: &Path,
    report: &mut Report,
) -> AnnotationStats {
    let genes_total = genes.len();
    let genes_with_symbol = genes
        .iter()
        .filter(|g| annotation.contains_key(&g.gene_id))
        .count();
    let stats = AnnotationStats {
        entries: annotation.len(),
        genes_total,
        genes_with_symbol,
    };
    let coverage = if genes_total > 0 {
        genes_with_symbol as f64 * 100.0 / genes_total as f64
    } else {
        0.0
    };
    report.push(
        Diagnostic::info(
            "ANNOTATION_COVERAGE",
            format!(
                "{} annotation entries; {genes_with_symbol}/{} genes have a symbol ({coverage:.1}%)",
                annotation.len(),
                genes_total
            ),
        )
        .with_file(path)
        .with_context("entries", annotation.len()),
    );
    stats
}

#[derive(Clone, Debug, Serialize)]
pub struct UserEvidenceStats {
    pub rows: usize,
    pub by_level: BTreeMap<String, usize>,
    pub by_type: BTreeMap<String, usize>,
}

pub fn review_user_list(rows: &[EdgeRecord], path: &Path, report: &mut Report) -> UserEvidenceStats {
    let mut by_level: BTreeMap<String, usize> = BTreeMap::new();
    let mut by_type: BTreeMap<String, usize> = BTreeMap::new();
    for row in rows {
        *by_level.entry(row.evidence_level.clone()).or_insert(0) += 1;
        *by_type.entry(row.evidence_type.clone()).or_insert(0) += 1;
    }
    let stats = UserEvidenceStats {
        rows: rows.len(),
        by_level: by_level.clone(),
        by_type: by_type.clone(),
    };
    report.push(
        Diagnostic::info(
            "USER_EVIDENCE_SUMMARY",
            format!("{} user evidence row(s)", rows.len()),
        )
        .with_file(path)
        .with_context("levels", format_levels(&by_level))
        .with_context("types", format_levels(&by_type)),
    );
    if rows.is_empty() {
        report.push(
            Diagnostic::warning("USER_EVIDENCE_EMPTY", "the user evidence list is empty")
                .with_file(path),
        );
    }
    stats
}

fn format_levels(map: &BTreeMap<String, usize>) -> String {
    map.iter()
        .map(|(k, v)| format!("{k}={v}"))
        .collect::<Vec<_>>()
        .join(",")
}

#[derive(Clone, Debug, Serialize)]
pub struct PromoterStats {
    pub promoters: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub truncated: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn report_serializes_and_counts() {
        let mut report = Report::new("validate");
        report.add_input("genome", std::path::Path::new("/tmp/g.fa"));
        report.push(
            Diagnostic::error("E1", "bad input")
                .with_file(std::path::Path::new("/tmp/g.fa"))
                .with_line(3)
                .with_context("count", 2)
                .with_suggestion("fix it"),
        );
        report.warning("W1", "check this");
        report.info("I1", "context");
        assert_eq!(report.summary.errors, 1);
        assert_eq!(report.summary.warnings, 1);
        assert_eq!(report.summary.infos, 1);
        assert!(report.has_errors());
        let path = std::env::temp_dir().join(format!(
            "regutrace_report_{}.json",
            std::process::id()
        ));
        report.write_json(&path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let value: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(value["summary"]["errors"], 1);
        assert_eq!(value["diagnostics"][0]["code"], "E1");
        assert_eq!(value["diagnostics"][0]["severity"], "error");
        assert_eq!(value["diagnostics"][0]["context"]["count"], "2");
        assert!(value["diagnostics"][0]["suggestion"]
            .as_str()
            .unwrap()
            .contains("fix"));
        assert!(report
            .human_lines()
            .iter()
            .any(|line| line.contains("[ERROR] E1")));
        let _ = std::fs::remove_file(path);
    }
}

pub fn review_promoters(
    promoters: &[Promoter],
    expected_len: Option<usize>,
    report: &mut Report,
) -> PromoterStats {
    let mut min_len = usize::MAX;
    let mut max_len = 0usize;
    let mut truncated = 0usize;
    for p in promoters {
        min_len = min_len.min(p.seq.len());
        max_len = max_len.max(p.seq.len());
        if let Some(expected) = expected_len {
            if p.seq.len() < expected {
                truncated += 1;
            }
        }
    }
    if min_len == usize::MAX {
        min_len = 0;
    }
    let stats = PromoterStats {
        promoters: promoters.len(),
        min_length: min_len,
        max_length: max_len,
        truncated,
    };
    report.push(
        Diagnostic::info(
            "PROMOTER_SUMMARY",
            format!(
                "{} promoter(s); lengths {}..{}",
                stats.promoters, stats.min_length, stats.max_length
            ),
        )
        .with_context("promoters", stats.promoters),
    );
    if promoters.is_empty() {
        report.push(
            Diagnostic::error("PROMOTER_EMPTY", "no promoters could be extracted")
                .with_suggestion("check chromosome names and TSS coordinates"),
        );
    } else if truncated > 0 {
        report.push(
            Diagnostic::warning(
                "PROMOTER_TRUNCATED",
                format!(
                    "{truncated} promoter(s) are shorter than the requested window",
                ),
            )
            .with_context("truncated", truncated)
            .with_suggestion("TSS near a chromosome boundary; reduce --promoter-len/--downstream or check the annotation"),
        );
    }
    stats
}
