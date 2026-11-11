//! Minimal GFF3 reader for gene/transcript models.
//!
//! Coordinates are kept 1-based inclusive as in GFF3. `tss` is the 1-based
//! coordinate of the transcription start site on the gene's strand.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::fasta::open_reader;

/// Summary of the features and inferred models in a GFF3 file. Used by the
/// input review and the JSON diagnostics report.
#[derive(Clone, Debug, Default)]
pub struct GffStats {
    pub feature_counts: BTreeMap<String, usize>,
    pub genes: usize,
    pub transcripts: usize,
    pub genes_without_transcript: usize,
    pub transcripts_without_features: usize,
    pub multi_isoform_genes: usize,
    pub inferred_genes: usize,
    pub inferred_transcripts: usize,
    pub skipped_records: usize,
    pub chromosomes: Vec<String>,
    pub strands: BTreeMap<String, usize>,
}

#[derive(Clone, Debug)]
pub struct GeneModel {
    pub gene_id: String,
    pub symbol: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub transcript_id: Option<String>,
    pub tss: u64,
}

#[derive(Clone, Debug)]
struct GeneInfo {
    symbol: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: char,
}

#[derive(Clone, Debug)]
struct TxInfo {
    transcript_id: String,
    start: u64,
    end: u64,
    chrom: String,
    strand: char,
}

fn attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let mut kv = part.splitn(2, '=');
        let k = kv.next()?.trim();
        if k != key {
            continue;
        }
        let v = kv.next()?.trim().trim_matches('"');
        if !v.is_empty() {
            return Some(v.to_string());
        }
    }
    None
}

fn clean_id(value: &str) -> String {
    let v = value.trim();
    for prefix in ["gene:", "transcript:", "mRNA:", "rna:"] {
        if let Some(rest) = v.strip_prefix(prefix) {
            return rest.to_string();
        }
    }
    v.to_string()
}

/// Representative-transcript length: exon sum when exons are present, else
/// CDS + UTR sum (spliced length), else the genomic span.
fn transcript_metric(
    tx: &TxInfo,
    exon_len: &HashMap<String, u64>,
    cds_len: &HashMap<String, u64>,
    utr_len: &HashMap<String, u64>,
) -> u64 {
    if let Some(&n) = exon_len.get(&tx.transcript_id) {
        if n > 0 {
            return n;
        }
    }
    let cds = cds_len.get(&tx.transcript_id).copied().unwrap_or(0);
    let utr = utr_len.get(&tx.transcript_id).copied().unwrap_or(0);
    if cds + utr > 0 {
        cds + utr
    } else {
        tx.end.saturating_sub(tx.start) + 1
    }
}

pub fn read_gff(path: &Path) -> Result<Vec<GeneModel>> {
    read_gff_with_stats(path).map(|(genes, _)| genes)
}

pub fn read_gff_with_stats(path: &Path) -> Result<(Vec<GeneModel>, GffStats)> {
    let reader = BufReader::new(open_reader(path)?);
    let mut genes: HashMap<String, GeneInfo> = HashMap::new();
    let mut tx: HashMap<String, Vec<TxInfo>> = HashMap::new();
    let mut exon_len: HashMap<String, u64> = HashMap::new();
    let mut cds_len: HashMap<String, u64> = HashMap::new();
    let mut utr_len: HashMap<String, u64> = HashMap::new();
    let mut feature_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut feature_span: HashMap<String, (String, char, u64, u64)> = HashMap::new();
    let mut feature_coding: HashSet<String> = HashSet::new();
    let mut inferred_genes = 0usize;
    let mut inferred_transcripts = 0usize;
    let mut skipped_records = 0usize;

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("read GFF line {}", lineno + 1))?;
        let line = line.trim_end_matches(|c| c == '\r' || c == '\n');
        if line.is_empty() || line.starts_with('#') {
            if line.starts_with("##FASTA") {
                break;
            }
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            bail!("GFF line {} has {} fields", lineno + 1, fields.len());
        }
        let chrom = fields[0].to_string();
        let kind = fields[2];
        let start: u64 = fields[3]
            .parse()
            .with_context(|| format!("GFF start at line {}", lineno + 1))?;
        let end: u64 = fields[4]
            .parse()
            .with_context(|| format!("GFF end at line {}", lineno + 1))?;
        if start == 0 || end < start {
            bail!(
                "GFF line {} has invalid coordinates {}-{}",
                lineno + 1,
                start,
                end
            );
        }
        let strand = fields[6].chars().next().unwrap_or('+');
        let attrs = fields[8];
        *feature_counts.entry(kind.to_string()).or_insert(0) += 1;

        if kind == "gene" {
            let Some(gene_id_raw) = attr(attrs, "ID").or_else(|| attr(attrs, "Name")) else {
                skipped_records += 1;
                continue;
            };
            let gene_id = clean_id(&gene_id_raw);
            let symbol = attr(attrs, "Name")
                .or_else(|| attr(attrs, "gene_name"))
                .unwrap_or_else(|| gene_id.clone());
            genes.insert(
                gene_id.clone(),
                GeneInfo {
                    symbol,
                    chrom,
                    start,
                    end,
                    strand,
                },
            );
        } else if kind == "mRNA" || kind == "transcript" {
            let Some(parent_raw) = attr(attrs, "Parent") else {
                skipped_records += 1;
                continue;
            };
            for parent_raw in parent_raw.split(',') {
                let parent = clean_id(parent_raw);
                if parent.is_empty() {
                    skipped_records += 1;
                    continue;
                }
                let transcript_id = attr(attrs, "ID")
                    .map(|v| clean_id(&v))
                    .unwrap_or_else(|| {
                        format!("{parent}.t{}", tx.get(&parent).map_or(1, |v| v.len() + 1))
                    });
                tx.entry(parent).or_default().push(TxInfo {
                    transcript_id,
                    start,
                    end,
                    chrom: chrom.clone(),
                    strand,
                });
            }
        } else if matches!(kind, "exon" | "CDS" | "five_prime_UTR" | "three_prime_UTR") {
            let Some(parent_raw) = attr(attrs, "Parent") else {
                skipped_records += 1;
                continue;
            };
            let len = end - start + 1;
            for parent_raw in parent_raw.split(',') {
                let transcript_id = clean_id(parent_raw);
                if transcript_id.is_empty() {
                    skipped_records += 1;
                    continue;
                }
                let entry = feature_span
                    .entry(transcript_id.clone())
                    .or_insert((chrom.clone(), strand, start, end));
                entry.2 = entry.2.min(start);
                entry.3 = entry.3.max(end);
                if kind == "exon" || kind == "CDS" {
                    feature_coding.insert(transcript_id.clone());
                }
                match kind {
                    "exon" => *exon_len.entry(transcript_id).or_insert(0) += len,
                    "CDS" => *cds_len.entry(transcript_id).or_insert(0) += len,
                    _ => *utr_len.entry(transcript_id).or_insert(0) += len,
                }
            }
        }
    }

    // Some GFFs have no gene/mRNA features and only CDS/exon records, or have
    // genes whose CDS records point directly at the gene ID. Infer one
    // transcript per coding parent so promoter extraction still works.
    // Reconcile every coding parent that is not already represented by an
    // explicit transcript. This covers CDS/exon-only GFFs, gene + CDS GFFs,
    // and partially annotated files where some genes have explicit models and
    // others only have CDS records.
    let represented: HashSet<String> = tx
        .values()
        .flatten()
        .map(|t| t.transcript_id.clone())
        .collect();
    for (parent, (chrom, strand, start, end)) in &feature_span {
        if !feature_coding.contains(parent) || represented.contains(parent) {
            continue;
        }
        let tx_info = TxInfo {
            transcript_id: parent.clone(),
            start: *start,
            end: *end,
            chrom: chrom.clone(),
            strand: *strand,
        };
        if genes.contains_key(parent) {
            tx.entry(parent.clone()).or_default().push(tx_info);
            inferred_transcripts += 1;
        } else {
            genes.insert(
                parent.clone(),
                GeneInfo {
                    symbol: parent.clone(),
                    chrom: chrom.clone(),
                    start: *start,
                    end: *end,
                    strand: *strand,
                },
            );
            tx.entry(parent.clone()).or_default().push(tx_info);
            inferred_genes += 1;
            inferred_transcripts += 1;
        }
    }

    // GFF3 without explicit gene features: use transcript parent as the gene.
    for (parent, txs) in &tx {
        if genes.contains_key(parent) {
            continue;
        }
        if let Some(first) = txs.first() {
            genes.insert(
                parent.clone(),
                GeneInfo {
                    symbol: parent.clone(),
                    chrom: first.chrom.clone(),
                    start: first.start,
                    end: first.end,
                    strand: first.strand,
                },
            );
            inferred_genes += 1;
        }
    }

    let mut chromosomes: BTreeSet<String> = BTreeSet::new();
    let mut strands: BTreeMap<String, usize> = BTreeMap::new();
    for g in genes.values() {
        chromosomes.insert(g.chrom.clone());
        *strands.entry(g.strand.to_string()).or_insert(0) += 1;
    }
    let transcript_count: usize = tx.values().map(|v| v.len()).sum();
    let stats = GffStats {
        feature_counts,
        genes: genes.len(),
        transcripts: transcript_count,
        genes_without_transcript: genes.keys().filter(|g| !tx.contains_key(*g)).count(),
        transcripts_without_features: tx
            .values()
            .flatten()
            .filter(|t| {
                !exon_len.contains_key(&t.transcript_id)
                    && !cds_len.contains_key(&t.transcript_id)
                    && !utr_len.contains_key(&t.transcript_id)
            })
            .count(),
        multi_isoform_genes: tx.values().filter(|v| v.len() > 1).count(),
        inferred_genes,
        inferred_transcripts,
        skipped_records,
        chromosomes: chromosomes.into_iter().collect(),
        strands,
    };

    let mut out = Vec::with_capacity(genes.len());
    for (gene_id, info) in genes {
        let strand = if info.strand == '+' || info.strand == '-' {
            info.strand
        } else {
            tx.get(&gene_id)
                .and_then(|v| v.first())
                .map(|t| t.strand)
                .unwrap_or('+')
        };
        let best = tx.get(&gene_id).and_then(|v| {
            v.iter().max_by_key(|t| {
                let metric = transcript_metric(t, &exon_len, &cds_len, &utr_len);
                // Ties prefer the most upstream TSS: smallest start on '+',
                // largest end on '-'.
                let tie = if strand == '-' {
                    t.end
                } else {
                    u64::MAX.saturating_sub(t.start)
                };
                (metric, tie)
            })
        });
        let (transcript_id, tss) = match best {
            Some(t) => (
                Some(t.transcript_id.clone()),
                if strand == '-' { t.end } else { t.start },
            ),
            None => (None, if strand == '-' { info.end } else { info.start }),
        };
        out.push(GeneModel {
            gene_id,
            symbol: info.symbol,
            chrom: info.chrom,
            start: info.start,
            end: info.end,
            strand,
            transcript_id,
            tss,
        });
    }
    out.sort_by(|a, b| a.gene_id.cmp(&b.gene_id));
    if out.is_empty() {
        bail!("no gene features found in {}", path.display());
    }
    Ok((out, stats))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_plus_and_minus_genes() {
        let path = std::env::temp_dir().join(format!("regutrace_gff_{}.gff3", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "##gff-version 3").unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene:AT1G00010;Name=PLUS").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=transcript:AT1G00010.1;Parent=gene:AT1G00010").unwrap();
            writeln!(f, "Chr2\ttest\tgene\t900\t1200\t.\t-\t.\tID=gene:AT2G00020;Name=MINUS").unwrap();
            writeln!(f, "Chr2\ttest\tmRNA\t900\t1100\t.\t-\t.\tID=transcript:AT2G00020.1;Parent=gene:AT2G00020").unwrap();
        }
        let genes = read_gff(&path).unwrap();
        let plus = genes.iter().find(|g| g.gene_id == "AT1G00010").unwrap();
        assert_eq!(plus.strand, '+');
        assert_eq!(plus.tss, 100);
        let minus = genes.iter().find(|g| g.gene_id == "AT2G00020").unwrap();
        assert_eq!(minus.strand, '-');
        assert_eq!(minus.tss, 1100);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn parses_transcripts_without_gene_features() {
        let path =
            std::env::temp_dir().join(format!("regutrace_gff_nogene_{}.gff3", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=geneA").unwrap();
            writeln!(f, "Chr2\ttest\ttranscript\t900\t1200\t.\t-\t.\tID=tx2;Parent=geneB")
                .unwrap();
        }
        let genes = read_gff(&path).unwrap();
        assert_eq!(genes.len(), 2);
        let a = genes.iter().find(|g| g.gene_id == "geneA").unwrap();
        assert_eq!(a.chrom, "Chr1");
        assert_eq!(a.strand, '+');
        assert_eq!(a.tss, 100);
        assert_eq!(a.transcript_id.as_deref(), Some("tx1"));
        let b = genes.iter().find(|g| g.gene_id == "geneB").unwrap();
        assert_eq!(b.chrom, "Chr2");
        assert_eq!(b.strand, '-');
        assert_eq!(b.tss, 1200);
        assert_eq!(b.transcript_id.as_deref(), Some("tx2"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn handles_multiple_parents() {
        let path =
            std::env::temp_dir().join(format!("regutrace_gff_multi_{}.gff3", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=geneA").unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=geneB").unwrap();
            writeln!(
                f,
                "Chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=geneA,geneB"
            )
            .unwrap();
        }
        let genes = read_gff(&path).unwrap();
        assert_eq!(genes.len(), 2);
        for id in ["geneA", "geneB"] {
            let g = genes.iter().find(|g| g.gene_id == id).unwrap();
            assert_eq!(g.transcript_id.as_deref(), Some("tx1"));
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn selects_longest_spliced_transcript() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_gff_isoform_{}.gff3",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=G1").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=G1").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t150\t300\t.\t+\t.\tID=tx2;Parent=G1").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=tx1").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t150\t300\t.\t+\t0\tParent=tx2").unwrap();
        }
        let genes = read_gff(&path).unwrap();
        let g = genes.iter().find(|g| g.gene_id == "G1").unwrap();
        assert_eq!(g.transcript_id.as_deref(), Some("tx2"));
        assert_eq!(g.tss, 150);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn selects_longest_exon_transcript_on_minus() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_gff_isoform_minus_{}.gff3",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t-\t.\tID=G2").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t500\t.\t-\t.\tID=tx1;Parent=G2").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t300\t.\t-\t.\tID=tx2;Parent=G2").unwrap();
            writeln!(f, "Chr1\ttest\texon\t100\t150\t.\t-\t.\tParent=tx1").unwrap();
            writeln!(f, "Chr1\ttest\texon\t400\t500\t.\t-\t.\tParent=tx1").unwrap();
            writeln!(f, "Chr1\ttest\texon\t100\t300\t.\t-\t.\tParent=tx2").unwrap();
        }
        let genes = read_gff(&path).unwrap();
        let g = genes.iter().find(|g| g.gene_id == "G2").unwrap();
        assert_eq!(g.transcript_id.as_deref(), Some("tx2"));
        assert_eq!(g.tss, 300);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn uses_cds_and_utr_when_exons_are_absent() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_gff_isoform_utr_{}.gff3",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=G3").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=G3").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t150\t300\t.\t+\t.\tID=tx2;Parent=G3").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=tx1").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t150\t250\t.\t+\t0\tParent=tx2").unwrap();
            writeln!(f, "Chr1\ttest\tfive_prime_UTR\t251\t300\t.\t+\t.\tParent=tx2").unwrap();
        }
        let genes = read_gff(&path).unwrap();
        let g = genes.iter().find(|g| g.gene_id == "G3").unwrap();
        assert_eq!(g.transcript_id.as_deref(), Some("tx2"));
        assert_eq!(g.tss, 150);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn infers_genes_from_cds_only_gff() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_gff_cdsonly_{}.gff3",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=geneA").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t300\t400\t.\t+\t0\tParent=geneA").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t500\t600\t.\t-\t0\tParent=geneB").unwrap();
        }
        let (genes, stats) = read_gff_with_stats(&path).unwrap();
        assert_eq!(stats.inferred_genes, 2);
        assert_eq!(stats.inferred_transcripts, 2);
        let a = genes.iter().find(|g| g.gene_id == "geneA").unwrap();
        assert_eq!(a.strand, '+');
        assert_eq!(a.tss, 100);
        let b = genes.iter().find(|g| g.gene_id == "geneB").unwrap();
        assert_eq!(b.strand, '-');
        assert_eq!(b.tss, 600);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn infers_pseudo_transcript_for_gene_with_cds_only() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_gff_gene_cds_{}.gff3",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=G4").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=G4").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t300\t400\t.\t+\t0\tParent=G4").unwrap();
        }
        let (genes, stats) = read_gff_with_stats(&path).unwrap();
        assert_eq!(stats.inferred_transcripts, 1);
        assert_eq!(stats.genes_without_transcript, 0);
        let g = genes.iter().find(|g| g.gene_id == "G4").unwrap();
        assert_eq!(g.transcript_id.as_deref(), Some("G4"));
        assert_eq!(g.tss, 100);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn reconciles_mixed_explicit_and_inferred_genes() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_gff_mixed_{}.gff3",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t100\t200\t.\t+\t0\tParent=t1").unwrap();
            writeln!(f, "Chr1\ttest\tCDS\t600\t700\t.\t+\t0\tParent=g2").unwrap();
        }
        let (genes, stats) = read_gff_with_stats(&path).unwrap();
        assert_eq!(genes.len(), 2);
        assert_eq!(stats.inferred_genes, 1);
        assert_eq!(stats.inferred_transcripts, 1);
        let g2 = genes.iter().find(|g| g.gene_id == "g2").unwrap();
        assert_eq!(g2.tss, 600);
        assert_eq!(g2.transcript_id.as_deref(), Some("g2"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn counts_transcript_only_inferred_gene() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_gff_txonly_{}.gff3",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1").unwrap();
        }
        let (genes, stats) = read_gff_with_stats(&path).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(stats.inferred_genes, 1);
        assert_eq!(stats.transcripts, 1);
        assert_eq!(genes[0].gene_id, "g1");
        assert_eq!(genes[0].tss, 100);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn counts_skipped_gene_records() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_gff_skipped_{}.gff3",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1").unwrap();
            writeln!(f, "Chr1\ttest\tgene\t600\t700\t.\t+\t.\tName=").unwrap();
        }
        let (genes, stats) = read_gff_with_stats(&path).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(stats.skipped_records, 1);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn transcript_metric_falls_back_to_span() {
        let tx = TxInfo {
            transcript_id: "t1".into(),
            start: 100,
            end: 300,
            chrom: "Chr1".into(),
            strand: '+',
        };
        assert_eq!(
            transcript_metric(&tx, &HashMap::new(), &HashMap::new(), &HashMap::new()),
            201
        );
    }

    #[test]
    fn rejects_invalid_coordinates() {
        let path =
            std::env::temp_dir().join(format!("regutrace_gff_bad_{}.gff3", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Chr1\ttest\tgene\t500\t100\t.\t+\t.\tID=geneA").unwrap();
        }
        let err = read_gff(&path).unwrap_err().to_string();
        assert!(err.contains("invalid coordinates"), "{err}");
        let _ = std::fs::remove_file(path);
    }
}
