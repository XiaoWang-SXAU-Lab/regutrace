use super::*;
use sha2::{Digest, Sha256};

// BLAST outfmt 6: qseqid sseqid bitscore pident length qlen slen evalue.
type BestHits = HashMap<String, (f64, BTreeSet<String>, bool)>;

/// Per-query alignment counts used by the Phase 0 diagnostics. Legacy mapping
/// consumes only `BestHits`; this audit distinguishes "no BLAST alignment" from
/// "alignments existed but none were family-compatible" without changing strict
/// mapping semantics.
#[derive(Clone, Debug, Default)]
struct HitAudit {
    all_alignments: usize,
    family_alignments: usize,
    best_eligible_score: Option<f64>,
}

#[derive(Clone, Debug)]
struct LocalDomainInfo {
    gene_id: String,
    family: String,
    models: Vec<String>,
    mode: &'static str,
    domain_count: usize,
}

fn choose_family_label<'a>(
    family_field: &'a str,
    model_by_family: &BTreeMap<String, DomainArchitecture>,
    evidence: &[EvidenceProtein],
) -> Option<&'a str> {
    let mut candidates: Vec<(&str, f64)> = family_field
        .split(';')
        .map(str::trim)
        .filter(|label| !label.is_empty())
        .filter_map(|label| {
            let architecture = model_by_family.get(label)?;
            let best = evidence
                .iter()
                .flat_map(|protein| protein.domains.iter())
                .filter(|domain| architecture.models().contains(&domain.model_name))
                .map(|domain| domain.domain_score)
                .fold(f64::NEG_INFINITY, f64::max);
            Some((label, best))
        })
        .collect();
    candidates.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(b.0))
    });
    candidates.into_iter().next().map(|(label, _)| label)
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct MappingRow {
    tf_id: String,
    tf_symbol: String,
    family: String,
    motif_id: String,
    reference_tf_id: String,
    mapping_method: String,
    domain_mode: String,
    extraction_family: String,
    extraction_models: String,
    domain_count: usize,
}

/// A catalogue row may carry several `;`-separated family labels. Match the
/// donor against any one of them; the old exact-string comparison silently
/// dropped multi-family rows such as `Hap4;bZIP`.
fn family_matches(catalogue: &str, donor: &str) -> bool {
    let donor = donor.trim();
    !donor.is_empty()
        && catalogue
            .split(';')
            .map(str::trim)
            .any(|family| !family.is_empty() && family == donor)
}

#[allow(dead_code)]
fn best_hits(path: &Path) -> Result<BestHits> {
    Ok(best_hits_with_audit(path)?.0)
}

fn best_hits_with_audit(path: &Path) -> Result<(BestHits, HashMap<String, HitAudit>)> {
    best_hits_by_audit(path, None)
}

type AlignmentIndex<'a> = (&'a HashMap<String, (String, String)>, &'a HashMap<String, (String, String)>);

#[allow(dead_code)]
fn best_hits_by(path: &Path, index: Option<AlignmentIndex<'_>>) -> Result<BestHits> {
    Ok(best_hits_by_audit(path, index)?.0)
}

fn best_hits_by_audit(
    path: &Path,
    index: Option<AlignmentIndex<'_>>,
) -> Result<(BestHits, HashMap<String, HitAudit>)> {
    let mut best = HashMap::<String, (f64, BTreeSet<String>, bool)>::new();
    let mut audit = HashMap::<String, HitAudit>::new();
    let mut rows = ReaderBuilder::new().delimiter(b'\t').has_headers(false)
        .from_reader(fasta::open_reader(path)?);
    for row in rows.records() {
        let row = row?;
        let raw_query = row.get(0).unwrap_or("").trim();
        let raw_subject = row.get(1).unwrap_or("").trim();
        if !matches!(row.len(), 8 | 12) || raw_query.is_empty() || raw_subject.is_empty() {
            bail!("expected eight or twelve BLAST columns in {}", path.display());
        }
        let (query, subject, compatible) = if let Some((queries, subjects)) = index {
            let (Some((query_id, query_family)), Some((subject_id, subject_family))) =
                (queries.get(raw_query), subjects.get(raw_subject)) else { continue };
            let compatible = family_matches(query_family, subject_family)
                || family_matches(subject_family, query_family);
            (query_id.as_str(), subject_id.as_str(), compatible)
        } else {
            (raw_query, raw_subject, true)
        };
        let audit_entry = audit.entry(query.to_owned()).or_default();
        audit_entry.all_alignments += 1;
        if !compatible {
            continue;
        }
        audit_entry.family_alignments += 1;
        let number = |i: usize| -> Result<f64> {
            let value = row.get(i).unwrap_or("").trim();
            let n: f64 = value
                .parse()
                .with_context(|| format!("parse BLAST column {i} ('{value}') in {}", path.display()))?;
            if !n.is_finite() || n < 0.0 {
                bail!("invalid number in {}", path.display());
            }
            Ok(n)
        };
        let score = number(2)?;
        let identity = number(3)?;
        let length = number(4)?;
        let qlen = number(5)?;
        let slen = number(6)?;
        let evalue = number(7)?;
        if qlen == 0.0 || slen == 0.0 || identity > 100.0 {
            bail!("invalid alignment lengths or identity in {}", path.display());
        }
        let (query_span, subject_span) = if row.len() == 12 {
            let span = |start, end, total| -> Result<f64> {
                let start = number(start)?;
                let end = number(end)?;
                if start < 1.0 || end < start || end > total {
                    bail!("invalid BLAST coordinates in {}", path.display());
                }
                Ok(end - start + 1.0)
            };
            (span(8, 9, qlen)?, span(10, 11, slen)?)
        } else {
            (length, length)
        };
        let (min_identity, min_coverage) = if index.is_some() { (35.0, 0.60) } else { (80.0, 0.80) };
        let eligible = identity >= min_identity && query_span / qlen >= min_coverage
            && subject_span / slen >= min_coverage && evalue <= 1e-5;
        if eligible {
            let slot = &mut audit_entry.best_eligible_score;
            *slot = Some(slot.map_or(score, |current| current.max(score)));
        }
        // Legacy strict semantics: rank by raw bitscore, then apply the
        // eligibility flag to that best row. A lower eligible row is not
        // substituted; Phase 0 diagnostics exposes this case explicitly.
        let item = best.entry(query.to_owned())
            .or_insert_with(|| (score, BTreeSet::new(), eligible));
        if score > item.0 {
            *item = (score, BTreeSet::new(), eligible);
        }
        if score == item.0 {
            item.1.insert(subject.to_owned());
            item.2 &= eligible;
        }
    }
    Ok((best, audit))
}

pub(super) fn run(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let inputs = [
        "--tf-genes", "--reference-map", "--forward", "--reverse", "--proteins",
        "--cds", "--reference-domains", "--dbd-models", "--reference-proteins",
    ];
    for key in opts.keys() {
        if !matches!(
            key.as_str(),
            "--out" | "--blastp" | "--diagnostics-out" | "--diagnostics-family-out"
        ) && !inputs.contains(&key.as_str()) {
            bail!("unknown option: {key}");
        }
    }
    let get = |key: &str| opts.get(key).map(String::as_str)
        .with_context(|| format!("missing {key}"));
    let out = get("--out")?;
    let diagnostics_out = opts.get("--diagnostics-out").map(Path::new);
    let diagnostics_family_out = opts.get("--diagnostics-family-out").map(Path::new);
    for input in inputs {
        if opts.get(input).is_some_and(|path| same_path(path, out)) {
            bail!("{input} and --out must be different files");
        }
    }
    let whole_protein = opts.contains_key("--reference-proteins");
    let external_alignment = opts.contains_key("--forward") || opts.contains_key("--reverse");
    let dbd_sequence_required = !whole_protein && !external_alignment;
    if whole_protein && ["--forward", "--reverse", "--reference-domains", "--dbd-models"]
        .iter().any(|key| opts.contains_key(*key)) {
        bail!("--reference-proteins requires automatic whole-protein alignment");
    }
    let genes = io::read_tf_genes_allow_empty(Path::new(get("--tf-genes")?))?;
    let reference = io::read_tf_map(Path::new(get("--reference-map")?))?;
    if reference.iter().any(|donor| donor.family.is_empty()) {
        bail!("reference map needs a nonempty family for every donor");
    }
    let (forward, reverse, local_info, forward_audit, _reverse_audit) =
        if opts.contains_key("--forward") || opts.contains_key("--reverse") {
            for option in ["--proteins", "--cds", "--reference-domains", "--dbd-models", "--reference-proteins", "--blastp"] {
                if opts.contains_key(option) {
                    bail!("{option} applies only to automatic alignment");
                }
            }
            let (forward, forward_audit) = best_hits_with_audit(Path::new(get("--forward")?))?;
            let (reverse, reverse_audit) = best_hits_with_audit(Path::new(get("--reverse")?))?;
            (forward, reverse, HashMap::new(), forward_audit, reverse_audit)
        } else {
            let (sequence_path, translate) = match (opts.get("--proteins"), opts.get("--cds")) {
                (Some(path), None) => (path, false),
                (None, Some(path)) => (path, true),
                _ => bail!("automatic alignment needs exactly one of --proteins or --cds"),
            };
            automatic_hits(
                Path::new(get("--tf-genes")?), Path::new(sequence_path), translate,
                Path::new(get(if whole_protein { "--reference-proteins" } else { "--reference-domains" })?),
                opts.get("--dbd-models").map(Path::new),
                whole_protein.then_some((&genes, reference.as_slice())),
                opts.get("--blastp").map_or("blastp", String::as_str), Path::new(out),
            )?
        };
    let mut mapped = BTreeSet::new();
    for donor in &reference {
        let Some((_, recipients, true)) = reverse.get(&donor.tf_id) else { continue };
        for recipient in recipients {
            let info = local_info.get(recipient);
            let gene_id = info.map_or(recipient.as_str(), |info| info.gene_id.as_str());
            let Some(gene) = genes.get(gene_id) else { continue };
            let Some((_, donors, true)) = forward.get(recipient) else { continue };
            if gene.is_identified()
                && family_matches(&gene.family, &donor.family)
                && donors.contains(&donor.tf_id)
            {
                let tied = recipients.len() > 1 || donors.len() > 1;
                let method = if whole_protein {
                    if tied { "reciprocal_protein_tied_candidate" } else { "reciprocal_protein_similarity" }
                } else if !tied {
                    "reciprocal_dbd_similarity"
                } else {
                    "reciprocal_dbd_tied_candidate"
                };
                mapped.insert(MappingRow {
                    tf_id: gene_id.to_owned(),
                    tf_symbol: gene.tf_symbol.clone(),
                    family: gene.family.clone(),
                    motif_id: donor.motif_id.clone(),
                    reference_tf_id: donor.tf_id.clone(),
                    mapping_method: method.to_owned(),
                    domain_mode: info.map_or(if whole_protein { "whole_protein" } else { "external_alignment" }, |info| info.mode).to_owned(),
                    extraction_family: info.map_or(String::new(), |info| info.family.clone()),
                    extraction_models: info.map_or(String::new(), |info| info.models.join("+")),
                    domain_count: info.map_or(0, |info| info.domain_count),
                });
            }
        }
    }
    let out = Path::new(out);
    let tmp = temp_output_path(out);
    let result = (|| -> Result<()> {
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&tmp)?;
        writer.write_record([
            "tf_id", "tf_symbol", "family", "motif_id", "reference_tf_id", "mapping_method",
            "domain_mode", "extraction_family", "extraction_models", "domain_count",
        ])?;
        for row in &mapped {
            writer.write_record([
                row.tf_id.as_str(),
                row.tf_symbol.as_str(),
                row.family.as_str(),
                row.motif_id.as_str(),
                row.reference_tf_id.as_str(),
                row.mapping_method.as_str(),
                row.domain_mode.as_str(),
                row.extraction_family.as_str(),
                row.extraction_models.as_str(),
                &row.domain_count.to_string(),
            ])?;
        }
        writer.flush()?;
        drop(writer);
        publish_temp(&tmp, out)
    })();
    if result.is_err() {
        cleanup_output_temp(&tmp);
    }
    result?;
    if let Some(path) = diagnostics_out {
        let (diagnostics, families) = build_mapping_diagnostics(
            &genes,
            &reference,
            &mapped,
            &local_info,
            &forward,
            &reverse,
            &forward_audit,
            whole_protein,
            dbd_sequence_required,
        );
        write_mapping_diagnostics(path, &diagnostics, diagnostics_family_out, &families)?;
    }
    println!("wrote {} inferred TF–motif mappings", mapped.len());
    Ok(())
}

#[derive(Clone, Debug)]
struct MappingDiagnostic {
    tf_id: String,
    tf_symbol: String,
    family: String,
    status: String,
    mapping_status: &'static str,
    primary_reason: &'static str,
    all_reasons: String,
    mapping_method: String,
    motif_id: String,
    reference_tf_id: String,
    forward_best_donor: String,
    forward_best_score: String,
    forward_best_eligible: String,
    reverse_best_recipients: String,
    family_donor_count: usize,
    family_motif_count: usize,
    notes: String,
}

#[derive(Clone, Debug, Default)]
struct FamilyDiagnostic {
    family: String,
    catalogue: usize,
    identified: usize,
    ambiguous: usize,
    rejected: usize,
    mapped_tfs: BTreeSet<String>,
    mapped_rows: usize,
    unmapped_identified: usize,
    donors: BTreeSet<String>,
    motifs: BTreeSet<String>,
    primary_reasons: BTreeMap<&'static str, usize>,
}

const REASON_ORDER: [&str; 13] = [
    "NOT_IDENTIFIED",
    "AMBIGUOUS_EXCLUDED",
    "REJECTED_EXCLUDED",
    "NO_DBD_SEQUENCE",
    "NO_REFERENCE_DONOR",
    "NO_BLAST_HIT",
    "FAMILY_MISMATCH",
    "BELOW_IDENTITY_OR_COVERAGE",
    "NOT_RECIPROCAL_BEST",
    "MULTIPLE_BEST_DONORS",
    "UNRESOLVED",
    "MAPPED_TIED",
    "MAPPED",
];

fn ordered_reasons(reasons: &[&'static str]) -> Vec<&'static str> {
    let mut out = reasons.to_vec();
    out.sort_by_key(|reason| {
        REASON_ORDER
            .iter()
            .position(|known| known == reason)
            .unwrap_or(REASON_ORDER.len())
    });
    out.dedup();
    out
}

fn joined_mapping_values(rows: &[&MappingRow], value: impl Fn(&MappingRow) -> String) -> String {
    let mut values = BTreeSet::new();
    for row in rows {
        let item = value(row);
        if !item.is_empty() {
            values.insert(item);
        }
    }
    values.into_iter().collect::<Vec<_>>().join("|")
}

fn build_mapping_diagnostics(
    genes: &HashMap<String, io::TfGene>,
    reference: &[io::TfMotif],
    mapped: &BTreeSet<MappingRow>,
    local_info: &HashMap<String, LocalDomainInfo>,
    forward: &BestHits,
    reverse: &BestHits,
    forward_audit: &HashMap<String, HitAudit>,
    whole_protein: bool,
    dbd_sequence_required: bool,
) -> (Vec<MappingDiagnostic>, Vec<FamilyDiagnostic>) {
    let mut mapped_by_tf = BTreeMap::<&str, Vec<&MappingRow>>::new();
    for row in mapped {
        mapped_by_tf.entry(row.tf_id.as_str()).or_default().push(row);
    }
    let mut diagnostics = Vec::with_capacity(genes.len());
    let mut family_stats = BTreeMap::<String, FamilyDiagnostic>::new();

    let mut gene_ids: Vec<&String> = genes.keys().collect();
    gene_ids.sort();
    for gene_id in gene_ids {
        let gene = &genes[gene_id];
        let rows = mapped_by_tf.get(gene.gene_id.as_str()).cloned().unwrap_or_default();
        let mut family_donors = BTreeSet::new();
        let mut family_motifs = BTreeSet::new();
        for donor in reference {
            if family_matches(&gene.family, &donor.family) {
                family_donors.insert(donor.tf_id.as_str());
                family_motifs.insert(donor.motif_id.as_str());
            }
        }
        let forward_entry = forward.get(gene.gene_id.as_str());
        let mut reverse_recipients = BTreeSet::new();
        if let Some((_, donors, _)) = forward_entry {
            for donor in donors {
                if let Some((_, recipients, _)) = reverse.get(donor) {
                    reverse_recipients.extend(recipients.iter().cloned());
                }
            }
        }

        let mut reasons: Vec<&'static str> = Vec::new();
        let mapping_status;
        if !rows.is_empty() {
            let tied = rows.iter().any(|row| row.mapping_method.contains("tied"));
            mapping_status = if tied { "MAPPED_TIED" } else { "MAPPED" };
            reasons.push(mapping_status);
        } else if gene.status == "ambiguous" {
            mapping_status = "UNMAPPED";
            reasons.push("AMBIGUOUS_EXCLUDED");
        } else if gene.status == "rejected" {
            mapping_status = "UNMAPPED";
            reasons.push("REJECTED_EXCLUDED");
        } else {
            mapping_status = "UNMAPPED";
            if dbd_sequence_required && !local_info.contains_key(gene.gene_id.as_str()) {
                reasons.push("NO_DBD_SEQUENCE");
            }
            if family_donors.is_empty() {
                reasons.push("NO_REFERENCE_DONOR");
            }
            if let Some((_, donors, eligible)) = forward_entry {
                if donors.is_empty() {
                    reasons.push("NO_BLAST_HIT");
                }
                if !*eligible {
                    reasons.push("BELOW_IDENTITY_OR_COVERAGE");
                }
                if donors.len() > 1 {
                    reasons.push("MULTIPLE_BEST_DONORS");
                }
                if !reverse_recipients.contains(gene.gene_id.as_str()) {
                    reasons.push("NOT_RECIPROCAL_BEST");
                }
                let compatible_donor = donors.iter().any(|donor_id| {
                    reference.iter().any(|donor| {
                        donor.tf_id == *donor_id && family_matches(&gene.family, &donor.family)
                    })
                });
                if !compatible_donor {
                    reasons.push("FAMILY_MISMATCH");
                }
            } else if let Some(audit) = forward_audit.get(gene.gene_id.as_str()) {
                if audit.all_alignments > 0 && audit.family_alignments == 0 {
                    reasons.push("FAMILY_MISMATCH");
                } else {
                    reasons.push("NO_BLAST_HIT");
                }
            } else {
                reasons.push("NO_BLAST_HIT");
            }
            if reasons.is_empty() {
                reasons.push("UNRESOLVED");
            }
        }
        let reasons = ordered_reasons(&reasons);
        let primary_reason = reasons.first().copied().unwrap_or("UNRESOLVED");
        let notes = if whole_protein {
            "sequence_mode=whole_protein; strict_reciprocal_protein_rbh"
        } else if dbd_sequence_required {
            "sequence_mode=dbd; legacy_reciprocal_dbd"
        } else {
            "sequence_mode=external_forward_reverse; legacy_reciprocal_dbd"
        };
        diagnostics.push(MappingDiagnostic {
            tf_id: gene.gene_id.clone(),
            tf_symbol: gene.tf_symbol.clone(),
            family: gene.family.clone(),
            status: gene.status.clone(),
            mapping_status,
            primary_reason,
            all_reasons: reasons.join("|"),
            mapping_method: joined_mapping_values(&rows, |row| row.mapping_method.clone()),
            motif_id: joined_mapping_values(&rows, |row| row.motif_id.clone()),
            reference_tf_id: joined_mapping_values(&rows, |row| row.reference_tf_id.clone()),
            forward_best_donor: forward_entry
                .map(|(_, donors, _)| donors.iter().cloned().collect::<Vec<_>>().join("|"))
                .unwrap_or_default(),
            forward_best_score: forward_entry
                .map(|(score, _, _)| format!("{score:.4}"))
                .unwrap_or_default(),
            forward_best_eligible: forward_entry
                .map(|(_, _, eligible)| eligible.to_string())
                .unwrap_or_default(),
            reverse_best_recipients: reverse_recipients.into_iter().collect::<Vec<_>>().join("|"),
            family_donor_count: family_donors.len(),
            family_motif_count: family_motifs.len(),
            notes: notes.to_owned(),
        });

        let family_key = if gene.family.is_empty() {
            "(unassigned)".to_owned()
        } else {
            gene.family.clone()
        };
        let stats = family_stats.entry(family_key.clone()).or_default();
        stats.family = family_key;
        stats.catalogue += 1;
        match gene.status.as_str() {
            "identified" => stats.identified += 1,
            "ambiguous" => stats.ambiguous += 1,
            "rejected" => stats.rejected += 1,
            _ => {}
        }
        stats.donors.extend(family_donors.iter().map(|value| (*value).to_owned()));
        stats.motifs.extend(family_motifs.iter().map(|value| (*value).to_owned()));
        if !rows.is_empty() {
            stats.mapped_tfs.insert(gene.gene_id.clone());
            stats.mapped_rows += rows.len();
        } else if gene.status == "identified" {
            stats.unmapped_identified += 1;
            *stats.primary_reasons.entry(primary_reason).or_insert(0) += 1;
        }
    }
    let mut families: Vec<FamilyDiagnostic> = family_stats.into_values().collect();
    families.sort_by(|a, b| a.family.cmp(&b.family));
    (diagnostics, families)
}

fn write_mapping_diagnostics(
    path: &Path,
    rows: &[MappingDiagnostic],
    family_path: Option<&Path>,
    families: &[FamilyDiagnostic],
) -> Result<()> {
    let tmp = temp_output_path(path);
    let result = (|| -> Result<()> {
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&tmp)?;
        writer.write_record([
            "tf_id", "tf_symbol", "family", "status", "mapping_status", "primary_reason",
            "all_reasons", "mapping_method", "motif_id", "reference_tf_id",
            "forward_best_donor", "forward_best_score", "forward_best_eligible",
            "reverse_best_recipients", "family_donor_count", "family_motif_count", "notes",
        ])?;
        for row in rows {
            writer.write_record([
                row.tf_id.as_str(), row.tf_symbol.as_str(), row.family.as_str(), row.status.as_str(),
                row.mapping_status, row.primary_reason, row.all_reasons.as_str(),
                row.mapping_method.as_str(), row.motif_id.as_str(), row.reference_tf_id.as_str(),
                row.forward_best_donor.as_str(), row.forward_best_score.as_str(),
                row.forward_best_eligible.as_str(), row.reverse_best_recipients.as_str(),
                &row.family_donor_count.to_string(), &row.family_motif_count.to_string(),
                row.notes.as_str(),
            ])?;
        }
        writer.flush()?;
        drop(writer);
        publish_temp(&tmp, path)
    })();
    if result.is_err() {
        cleanup_output_temp(&tmp);
    }
    result?;

    if let Some(family_path) = family_path {
        let tmp = temp_output_path(family_path);
        let result = (|| -> Result<()> {
            let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&tmp)?;
            writer.write_record([
                "family", "catalogue", "identified", "ambiguous", "rejected",
                "mapped_tfs", "mapped_rows", "unmapped_identified", "donor_count", "motif_count",
                "reason_no_reference_donor", "reason_no_dbd_sequence", "reason_no_blast_hit",
                "reason_family_mismatch", "reason_below_identity_or_coverage",
                "reason_not_reciprocal_best", "reason_multiple_best_donors", "reason_unresolved",
            ])?;
            for family in families {
                let reason = |code: &str| family.primary_reasons.get(code).copied().unwrap_or(0);
                writer.write_record([
                    family.family.as_str(),
                    &family.catalogue.to_string(),
                    &family.identified.to_string(),
                    &family.ambiguous.to_string(),
                    &family.rejected.to_string(),
                    &family.mapped_tfs.len().to_string(),
                    &family.mapped_rows.to_string(),
                    &family.unmapped_identified.to_string(),
                    &family.donors.len().to_string(),
                    &family.motifs.len().to_string(),
                    &reason("NO_REFERENCE_DONOR").to_string(),
                    &reason("NO_DBD_SEQUENCE").to_string(),
                    &reason("NO_BLAST_HIT").to_string(),
                    &reason("FAMILY_MISMATCH").to_string(),
                    &reason("BELOW_IDENTITY_OR_COVERAGE").to_string(),
                    &reason("NOT_RECIPROCAL_BEST").to_string(),
                    &reason("MULTIPLE_BEST_DONORS").to_string(),
                    &reason("UNRESOLVED").to_string(),
                ])?;
            }
            writer.flush()?;
            drop(writer);
            publish_temp(&tmp, family_path)
        })();
        if result.is_err() {
            cleanup_output_temp(&tmp);
        }
        result?;
    }
    Ok(())
}

fn automatic_hits(
    catalogue: &Path,
    sequences: &Path,
    translate: bool,
    reference: &Path,
    models: Option<&Path>,
    whole_protein: Option<(&HashMap<String, io::TfGene>, &[io::TfMotif])>,
    blastp: &str,
    out: &Path,
) -> Result<(BestHits, BestHits, HashMap<String, LocalDomainInfo>, HashMap<String, HitAudit>, HashMap<String, HitAudit>)> {
    let tmp = temp_output_path(out);
    fs::create_dir(&tmp)?;
    let paths = ["translated.faa", "local.faa", "reference.faa", "forward.tsv", "reverse.tsv"]
        .map(|name| tmp.join(name));
    let [translated, local, donor, forward, reverse] = &paths;
    let result = (|| {
        let input = if translate {
            plant_tfs::translate_cds(sequences, translated)?;
            translated.as_path()
        } else {
            sequences
        };
        let mut local_index = HashMap::new();
        let mut donor_index = HashMap::new();
        let info = if let Some((genes, reference_map)) = whole_protein {
            let proteins = fasta::read_fasta(input)?;
            let mut rows = ReaderBuilder::new().delimiter(b'\t').from_path(catalogue)?;
            let evidence_col = rows.headers()?.iter().position(|h| h == "evidence_json")
                .context("TF catalogue missing evidence_json")?;
            let mut output = BufWriter::new(File::create(local)?);
            for row in rows.records() {
                let row = row?;
                let Some(gene) = row.get(0).and_then(|id| genes.get(id)).filter(|g| g.is_identified())
                else { continue };
                let entries: Vec<EvidenceProtein> = serde_json::from_str(row.get(evidence_col).unwrap_or(""))?;
                for entry in entries {
                    let seq = proteins.seqs.get(&entry.protein_id)
                        .with_context(|| format!("TF protein {} missing from inputs", entry.protein_id))?;
                    local_index.insert(entry.protein_id.clone(), (gene.gene_id.clone(), gene.family.clone()));
                    writeln!(output, ">{}", entry.protein_id)?;
                    output.write_all(seq)?;
                    writeln!(output)?;
                }
            }
            output.flush()?;
            for donor in reference_map {
                if donor_index.insert(donor.tf_id.clone(), (donor.tf_id.clone(), donor.family.clone())).is_some() {
                    bail!("whole-protein reference map requires one representative matrix per donor: {}", donor.tf_id);
                }
            }
            HashMap::new()
        } else {
            write_local_domains(catalogue, input, models.context("missing --dbd-models")?, local)?
        };
        if info.is_empty() && local_index.is_empty() {
            return Ok((
                BestHits::new(),
                BestHits::new(),
                HashMap::new(),
                HashMap::new(),
                HashMap::new(),
            ));
        }
        // BLAST accepts plain FASTA. Stage a decompressed copy for either input format.
        std::io::copy(&mut fasta::open_reader(reference)?, &mut File::create(donor)?)?;
        let donors = fasta::read_fasta(donor)?;
        let locals = fasta::read_fasta(local)?;
        if whole_protein.is_some() && (donors.order.is_empty()
            || donor_index.keys().any(|id| !donors.seqs.contains_key(id))) {
            bail!("reference protein FASTA does not contain every representative motif donor");
        }
        for (query, subject, table, subject_count) in [
            (local, donor, forward, donors.order.len()),
            (donor, local, reverse, locals.order.len()),
        ] {
            let result = std::process::Command::new(blastp)
                .arg("-query").arg(query)
                .arg("-subject").arg(subject)
                .args(["-outfmt", "6 qseqid sseqid bitscore pident length qlen slen evalue qstart qend sstart send"])
                .args(["-max_target_seqs", &subject_count.to_string(), "-evalue", "1e-5"])
                .arg("-out").arg(table)
                .output().with_context(|| format!("start {blastp}; install BLAST+ or pass --blastp"))?;
            if !result.status.success() {
                bail!("blastp failed ({}): {}", result.status,
                    String::from_utf8_lossy(&result.stderr).trim());
            }
        }
        if whole_protein.is_some() {
            let (forward_hits, forward_audit) =
                best_hits_by_audit(forward, Some((&local_index, &donor_index)))?;
            let (reverse_hits, reverse_audit) =
                best_hits_by_audit(reverse, Some((&donor_index, &local_index)))?;
            Ok((forward_hits, reverse_hits, info, forward_audit, reverse_audit))
        } else {
            let (forward_hits, forward_audit) = best_hits_with_audit(forward)?;
            let (reverse_hits, reverse_audit) = best_hits_with_audit(reverse)?;
            Ok((forward_hits, reverse_hits, info, forward_audit, reverse_audit))
        }
    })();
    for path in paths {
        let _ = fs::remove_file(path);
    }
    let _ = fs::remove_dir(tmp);
    result
}

#[derive(serde::Deserialize)]
struct EvidenceProtein {
    protein_id: String,
    domains: Vec<EvidenceDomain>,
}

#[derive(serde::Deserialize)]
struct EvidenceDomain {
    model_name: String,
    protein_start: usize,
    protein_end: usize,
    #[serde(default)]
    domain_score: f64,
}

#[derive(Clone, Debug, serde::Deserialize)]
#[serde(untagged)]
enum DomainArchitecture {
    Single(String),
    Ordered(Vec<String>),
}

#[derive(Clone, Debug, Default)]
struct DbdModelSet {
    version: String,
    source: String,
    families: BTreeMap<String, DomainArchitecture>,
    abstained: BTreeMap<String, String>,
}

#[derive(serde::Deserialize)]
struct DetailedDbdFamily {
    models: DomainArchitecture,
    #[serde(default)]
    status: String,
    #[serde(default)]
    source: String,
    #[serde(default)]
    notes: String,
}

/// Read either the legacy `{family: model-or-array}` map or the versioned
/// `{version, source, families: {family: {models, status, source}}}` form.
/// Abstained families are retained for diagnostics but never extracted.
fn read_dbd_model_set(path: &Path) -> Result<DbdModelSet> {
    let value: serde_json::Value = serde_json::from_reader(fasta::open_reader(path)?)
        .with_context(|| format!("read DBD models {}", path.display()))?;
    let object = value
        .as_object()
        .context("DBD model file must be a JSON object")?;
    let mut set = DbdModelSet::default();
    if let Some(families) = object.get("families") {
        set.version = object
            .get("version")
            .and_then(|value| value.as_str())
            .unwrap_or("")
            .to_owned();
        set.source = object
            .get("source")
            .and_then(|value| value.as_str())
            .unwrap_or("")
            .to_owned();
        let families = families
            .as_object()
            .context("DBD model file families must be an object")?;
        for (family, entry) in families {
            if entry.is_string() || entry.is_array() {
                let architecture: DomainArchitecture = serde_json::from_value(entry.clone())
                    .with_context(|| format!("read DBD architecture for {family}"))?;
                set.families.insert(family.clone(), architecture);
                continue;
            }
            let detailed: DetailedDbdFamily = serde_json::from_value(entry.clone())
                .with_context(|| format!("read DBD family config for {family}"))?;
            let status = if detailed.status.is_empty() {
                "supported".to_owned()
            } else {
                detailed.status.clone()
            };
            match status.as_str() {
                "supported" | "calibrated" => {
                    set.families.insert(family.clone(), detailed.models);
                }
                "abstain" => {
                    let reason = if detailed.notes.is_empty() {
                        detailed.source
                    } else {
                        detailed.notes
                    };
                    set.abstained.insert(family.clone(), reason);
                }
                other => bail!("unknown DBD family status '{other}' for {family}"),
            }
        }
    } else {
        set.families = serde_json::from_value(value)?;
    }
    if set.families.is_empty() {
        bail!("DNA-binding domain model list is empty");
    }
    for (family, architecture) in &set.families {
        if architecture.models().is_empty()
            || architecture
                .models()
                .iter()
                .any(|name| name.trim().is_empty())
        {
            bail!("empty domain architecture for {family}");
        }
    }
    Ok(set)
}

impl DomainArchitecture {
    fn models(&self) -> &[String] {
        match self {
            Self::Single(model) => std::slice::from_ref(model),
            Self::Ordered(models) => models,
        }
    }
}

fn select_architecture_domains<'a>(
    matching: &[&'a EvidenceDomain],
    models: &[String],
) -> Option<Vec<&'a EvidenceDomain>> {
    let mut used = vec![false; matching.len()];
    let mut selected = Vec::with_capacity(models.len());
    for model in models {
        let (index, _) = matching
            .iter()
            .enumerate()
            .filter(|(index, hit)| !used[*index] && hit.model_name == *model)
            .max_by(|(_, a), (_, b)| {
                a.domain_score
                    .partial_cmp(&b.domain_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })?;
        used[index] = true;
        selected.push(matching[index]);
    }
    selected.sort_by_key(|hit| (hit.protein_start, hit.protein_end));
    if selected
        .windows(2)
        .any(|pair| pair[0].protein_end >= pair[1].protein_start)
    {
        return None;
    }
    Some(selected)
}

#[derive(Clone, Debug)]
struct ExtractedDbd {
    gene_id: String,
    family: String,
    models: Vec<String>,
    start: usize,
    end: usize,
    mode: &'static str,
    domain_count: usize,
    score: f64,
    sequence: Vec<u8>,
}

fn dbd_sequence_checksum(sequence: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(sequence);
    format!("sha256:{:x}", hasher.finalize())
}

fn extract_dbds_from_catalogue(
    catalogue: &Path,
    sequences: &Path,
    models: &Path,
    include_ambiguous: bool,
) -> Result<(Vec<ExtractedDbd>, DbdModelSet)> {
    let model_set = read_dbd_model_set(models)?;
    let proteins = fasta::read_fasta(sequences)?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(fasta::open_reader(catalogue)?);
    let headers = reader.headers()?.clone();
    let column = |name: &str| {
        headers
            .iter()
            .position(|col| col == name)
            .with_context(|| format!("TF catalogue needs {name} for DBD extraction"))
    };
    let gene = column("gene_id")?;
    let status = column("identification_status")?;
    let family = column("family")?;
    let evidence = column("evidence_json")?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let id = &row[gene];
        if !matches!(&row[status], "identified" | "ambiguous") {
            continue;
        }
        if !include_ambiguous && &row[status] != "identified" {
            continue;
        }
        let evidence: Vec<EvidenceProtein> = serde_json::from_str(&row[evidence])
            .with_context(|| format!("read evidence for {id}"))?;
        let Some(family_label) = choose_family_label(&row[family], &model_set.families, &evidence)
        else {
            continue;
        };
        let architecture = model_set
            .families
            .get(family_label)
            .expect("chosen family label has a configured architecture");
        let model_names = architecture.models().to_vec();
        let mut domains =
            BTreeMap::<Vec<u8>, (&'static str, usize, f64, usize, usize)>::new();
        for protein in evidence {
            let sequence = proteins
                .seqs
                .get(&protein.protein_id)
                .with_context(|| format!("{} is absent from input FASTA", protein.protein_id))?;
            let mut matching: Vec<_> = protein
                .domains
                .iter()
                .filter(|hit| model_names.contains(&hit.model_name))
                .collect();
            matching.sort_by_key(|hit| (hit.protein_start, hit.protein_end));
            let architecture = select_architecture_domains(&matching, &model_names);
            let (start, end, mode, domain_count, score) = if let Some(selected) = architecture {
                let first = selected.first().context("empty architecture")?;
                let last = selected.last().context("empty architecture")?;
                let score = selected
                    .iter()
                    .map(|hit| hit.domain_score)
                    .fold(f64::NEG_INFINITY, f64::max);
                (
                    first.protein_start,
                    last.protein_end,
                    "full_architecture",
                    model_names.len(),
                    score,
                )
            } else {
                // Multi-domain TFs do not satisfy every configured architecture.
                // Fall back to the best-scoring single domain so the TF still
                // reaches the alignment competition instead of being silently
                // dropped. The output labels this as single-domain evidence.
                let Some(best) = matching.iter().max_by(|a, b| {
                    a.domain_score
                        .partial_cmp(&b.domain_score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                }) else {
                    continue;
                };
                (
                    best.protein_start,
                    best.protein_end,
                    "single_domain_fallback",
                    1,
                    best.domain_score,
                )
            };
            let dbd_start = start;
            let start = start
                .checked_sub(1)
                .with_context(|| format!("invalid domain start for {id}"))?;
            let slice = sequence
                .get(start..end)
                .filter(|slice| !slice.is_empty())
                .with_context(|| format!("invalid domain coordinates for {id}"))?;
            let entry = domains
                .entry(slice.to_vec())
                .or_insert((mode, domain_count, score, dbd_start, end));
            if score > entry.2 {
                *entry = (mode, domain_count, score, dbd_start, end);
            }
        }
        // Keep ambiguous genes in the competition when they have a configured
        // label, but only the final mapping step can accept identified genes.
        // Isoforms can carry distinct but valid DBD sequences; keep the
        // best-scoring representative instead of dropping the whole gene.
        if let Some((sequence, (mode, domain_count, score, dbd_start, dbd_end))) = domains
            .into_iter()
            .max_by(|(seq_a, a), (seq_b, b)| {
                a.2.partial_cmp(&b.2)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| seq_a.len().cmp(&seq_b.len()))
            })
        {
            out.push(ExtractedDbd {
                gene_id: id.to_string(),
                family: family_label.to_string(),
                models: model_names,
                start: dbd_start,
                end: dbd_end,
                mode,
                domain_count,
                score,
                sequence,
            });
        }
    }
    Ok((out, model_set))
}

fn write_local_domains(
    catalogue: &Path,
    sequences: &Path,
    models: &Path,
    output: &Path,
) -> Result<HashMap<String, LocalDomainInfo>> {
    let (dbds, _) = extract_dbds_from_catalogue(catalogue, sequences, models, true)?;
    let mut writer = BufWriter::new(File::create(output)?);
    let mut info_by_gene = HashMap::new();
    for dbd in dbds {
        writeln!(writer, ">{}", dbd.gene_id)?;
        writer.write_all(&dbd.sequence)?;
        writeln!(writer)?;
        info_by_gene.insert(
            dbd.gene_id.clone(),
            LocalDomainInfo {
                gene_id: dbd.gene_id,
                family: dbd.family,
                models: dbd.models,
                mode: dbd.mode,
                domain_count: dbd.domain_count,
            },
        );
    }
    writer.flush()?;
    Ok(info_by_gene)
}

pub(super) fn run_dbd_extract(args: &[String]) -> Result<()> {
    let opts = parse_options_with_flags(args, &["--include-ambiguous"])?;
    for key in opts.keys() {
        if !matches!(
            key.as_str(),
            "--tf-genes"
                | "--proteins"
                | "--dbd-models"
                | "--out-faa"
                | "--out-tsv"
                | "--include-ambiguous"
        ) {
            bail!("unknown option: {key}");
        }
    }
    let get = |key: &str| {
        opts.get(key)
            .map(String::as_str)
            .with_context(|| format!("missing {key}"))
    };
    let include_ambiguous = opts
        .get("--include-ambiguous")
        .is_none_or(|value| !matches!(value.to_ascii_lowercase().as_str(), "false" | "0"));
    let (dbds, model_set) = extract_dbds_from_catalogue(
        Path::new(get("--tf-genes")?),
        Path::new(get("--proteins")?),
        Path::new(get("--dbd-models")?),
        include_ambiguous,
    )?;
    let out_faa = Path::new(get("--out-faa")?);
    let out_tsv = Path::new(get("--out-tsv")?);
    let mut fasta_writer = BufWriter::new(File::create(out_faa)?);
    let mut tsv = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(out_tsv)?;
    tsv.write_record([
        "sequence_id",
        "gene_id",
        "family",
        "models",
        "dbd_start",
        "dbd_end",
        "domain_count",
        "mode",
        "score",
        "sequence_sha256",
        "sequence_length",
    ])?;
    for dbd in &dbds {
        writeln!(fasta_writer, ">{}", dbd.gene_id)?;
        fasta_writer.write_all(&dbd.sequence)?;
        writeln!(fasta_writer)?;
        tsv.write_record([
            dbd.gene_id.as_str(),
            dbd.gene_id.as_str(),
            dbd.family.as_str(),
            &dbd.models.join("+"),
            &dbd.start.to_string(),
            &dbd.end.to_string(),
            &dbd.domain_count.to_string(),
            dbd.mode,
            &format!("{:.4}", dbd.score),
            &dbd_sequence_checksum(&dbd.sequence),
            &dbd.sequence.len().to_string(),
        ])?;
    }
    fasta_writer.flush()?;
    tsv.flush()?;
    println!(
        "extracted {} DBD sequence(s) with {} ({} supported families, {} abstained)",
        dbds.len(),
        if model_set.version.is_empty() {
            "legacy DBD models"
        } else {
            &model_set.version
        },
        model_set.families.len(),
        model_set.abstained.len()
    );
    Ok(())
}

#[derive(Clone, Debug, serde::Deserialize)]
struct DbdThresholdValues {
    min_identity: f64,
    min_coverage: f64,
    max_evalue: f64,
    #[serde(default)]
    source: String,
}

#[derive(Clone, Debug, serde::Deserialize)]
struct DbdThresholdConfig {
    version: String,
    #[serde(default)]
    calibration_status: String,
    #[serde(default)]
    source: String,
    #[serde(default)]
    defaults: Option<DbdThresholdValues>,
    #[serde(default)]
    families: BTreeMap<String, DbdThresholdValues>,
}

impl DbdThresholdConfig {
    fn validate_value(&self, label: &str, value: &DbdThresholdValues) -> Result<()> {
        if !value.min_identity.is_finite()
            || value.min_identity < 0.0
            || value.min_identity > 100.0
        {
            bail!("{label}: min_identity must be within 0..=100");
        }
        if !value.min_coverage.is_finite()
            || value.min_coverage <= 0.0
            || value.min_coverage > 1.0
        {
            bail!("{label}: min_coverage must be within (0, 1]");
        }
        if !value.max_evalue.is_finite() || value.max_evalue <= 0.0 {
            bail!("{label}: max_evalue must be positive");
        }
        Ok(())
    }

    fn validate(&self) -> Result<()> {
        if self.version.trim().is_empty() {
            bail!("DBD threshold config needs a version");
        }
        let status = self.calibration_status.trim().to_ascii_lowercase();
        if !matches!(status.as_str(), "calibrated" | "uncalibrated") {
            bail!("DBD threshold config calibration_status must be calibrated or uncalibrated");
        }
        if self.defaults.is_none() && self.families.is_empty() {
            bail!("DBD threshold config has neither defaults nor family thresholds");
        }
        if let Some(defaults) = &self.defaults {
            self.validate_value("defaults", defaults)?;
        }
        for (family, values) in &self.families {
            self.validate_value(family, values)?;
        }
        Ok(())
    }

    fn values_for(&self, family_field: &str) -> Option<&DbdThresholdValues> {
        for family in family_field.split(';').map(str::trim) {
            if let Some(values) = self.families.get(family) {
                return Some(values);
            }
        }
        self.defaults.as_ref()
    }

    fn max_evalue(&self) -> f64 {
        self.families
            .values()
            .chain(self.defaults.iter())
            .map(|values| values.max_evalue)
            .fold(0.0_f64, f64::max)
    }

    fn is_uncalibrated(&self) -> bool {
        self.calibration_status.trim().eq_ignore_ascii_case("uncalibrated")
    }
}

fn read_dbd_thresholds(path: &Path) -> Result<DbdThresholdConfig> {
    let config: DbdThresholdConfig = serde_json::from_reader(fasta::open_reader(path)?)
        .with_context(|| format!("read DBD threshold config {}", path.display()))?;
    config.validate()?;
    Ok(config)
}

#[derive(Clone, Debug)]
struct DonorDbd {
    donor_id: String,
    family: String,
    models: String,
    start: usize,
    mode: String,
    sequence: Vec<u8>,
}

fn read_donor_dbds(tsv: &Path, fasta_path: &Path) -> Result<BTreeMap<String, DonorDbd>> {
    let sequences = fasta::read_fasta(fasta_path)?;
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv)?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| {
        headers
            .iter()
            .position(|header| header.trim().eq_ignore_ascii_case(name))
    };
    let sequence_id = col("sequence_id").context("donor DBD table needs sequence_id")?;
    let family = col("family").context("donor DBD table needs family")?;
    let start = col("dbd_start").context("donor DBD table needs dbd_start")?;
    let end = col("dbd_end").context("donor DBD table needs dbd_end")?;
    let checksum = col("sequence_sha256").context("donor DBD table needs sequence_sha256")?;
    let models = col("models");
    let mode = col("mode");
    let sequence_length = col("sequence_length");
    let cell = |record: &csv::StringRecord, index: Option<usize>| {
        index
            .and_then(|i| record.get(i))
            .map(str::trim)
            .unwrap_or("")
            .to_owned()
    };
    let mut out = BTreeMap::new();
    for record in reader.records() {
        let record = record?;
        let donor_id = record.get(sequence_id).unwrap_or("").trim().to_owned();
        if donor_id.is_empty() {
            bail!("donor DBD table contains an empty sequence_id");
        }
        let sequence = sequences
            .seqs
            .get(&donor_id)
            .with_context(|| format!("donor DBD {donor_id} is absent from the DBD FASTA"))?;
        let start_value = record
            .get(start)
            .unwrap_or("")
            .trim()
            .parse::<usize>()
            .with_context(|| format!("donor DBD {donor_id}: invalid dbd_start"))?;
        let end_value = record
            .get(end)
            .unwrap_or("")
            .trim()
            .parse::<usize>()
            .with_context(|| format!("donor DBD {donor_id}: invalid dbd_end"))?;
        if start_value == 0 || end_value < start_value {
            bail!(
                "donor DBD {donor_id}: invalid coordinates {start_value}-{end_value}"
            );
        }
        if end_value - start_value + 1 != sequence.len() {
            bail!(
                "donor DBD {donor_id}: FASTA length {} does not match full-protein coordinates {start_value}-{end_value}",
                sequence.len()
            );
        }
        let expected = record.get(checksum).unwrap_or("").trim();
        let actual = dbd_sequence_checksum(sequence);
        if expected != actual {
            bail!(
                "donor DBD {donor_id}: checksum mismatch (table {expected}, FASTA {actual})"
            );
        }
        let declared_length = cell(&record, sequence_length);
        if !declared_length.is_empty()
            && declared_length.parse::<usize>().ok() != Some(sequence.len())
        {
            bail!("donor DBD {donor_id}: declared sequence_length does not match FASTA");
        }
        let entry = DonorDbd {
            donor_id: donor_id.clone(),
            family: cell(&record, Some(family)),
            models: cell(&record, models),
            start: start_value,
            mode: cell(&record, mode),
            sequence: sequence.clone(),
        };
        if entry.family.is_empty() {
            bail!("donor DBD {donor_id}: empty family");
        }
        if out.insert(donor_id.clone(), entry).is_some() {
            bail!("duplicate donor DBD record: {donor_id}");
        }
    }
    Ok(out)
}

#[derive(Clone, Debug)]
struct DbdAlignment {
    query: String,
    subject: String,
    bitscore: f64,
    identity: f64,
    alignment_length: f64,
    evalue: f64,
    query_start: f64,
    query_end: f64,
    subject_start: f64,
    subject_end: f64,
    query_coverage: f64,
    subject_coverage: f64,
}

fn read_dbd_alignments(path: &Path) -> Result<Vec<DbdAlignment>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(fasta::open_reader(path)?);
    let mut out = Vec::new();
    for record in reader.records() {
        let record = record?;
        if record.len() != 12 {
            bail!(
                "expected twelve BLAST columns in {}, found {}",
                path.display(),
                record.len()
            );
        }
        let number = |index: usize| -> Result<f64> {
            let value = record.get(index).unwrap_or("").trim();
            let parsed = value.parse::<f64>().with_context(|| {
                format!("parse BLAST column {index} ('{value}') in {}", path.display())
            })?;
            if !parsed.is_finite() || parsed < 0.0 {
                bail!("invalid number in {}", path.display());
            }
            Ok(parsed)
        };
        let query = record.get(0).unwrap_or("").trim().to_owned();
        let subject = record.get(1).unwrap_or("").trim().to_owned();
        if query.is_empty() || subject.is_empty() {
            bail!("empty query or subject in {}", path.display());
        }
        let bitscore = number(2)?;
        let identity = number(3)?;
        let alignment_length = number(4)?;
        let query_length = number(5)?;
        let subject_length = number(6)?;
        let evalue = number(7)?;
        if query_length == 0.0
            || subject_length == 0.0
            || identity > 100.0
            || alignment_length == 0.0
        {
            bail!("invalid alignment metrics in {}", path.display());
        }
        let span = |start_index: usize, end_index: usize, total: f64| -> Result<f64> {
            let start = number(start_index)?;
            let end = number(end_index)?;
            if start < 1.0 || end < start || end > total {
                bail!("invalid BLAST coordinates in {}", path.display());
            }
            Ok(end - start + 1.0)
        };
        let query_span = span(8, 9, query_length)?;
        let subject_span = span(10, 11, subject_length)?;
        out.push(DbdAlignment {
            query,
            subject,
            bitscore,
            identity,
            alignment_length,
            evalue,
            query_start: number(8)?,
            query_end: number(9)?,
            subject_start: number(10)?,
            subject_end: number(11)?,
            query_coverage: query_span / query_length,
            subject_coverage: subject_span / subject_length,
        });
    }
    Ok(out)
}

fn read_specificity_groups(path: &Path) -> Result<HashMap<String, String>> {
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    let headers = reader.headers()?.clone();
    let motif = headers
        .iter()
        .position(|header| header.trim().eq_ignore_ascii_case("motif_id"))
        .context("specificity group table needs motif_id")?;
    let group = headers
        .iter()
        .position(|header| header.trim().eq_ignore_ascii_case("group_090"))
        .context("specificity group table needs group_090")?;
    let mut out = HashMap::new();
    for record in reader.records() {
        let record = record?;
        let motif_id = record.get(motif).unwrap_or("").trim();
        let group_id = record.get(group).unwrap_or("").trim();
        if motif_id.is_empty() || group_id.is_empty() {
            continue;
        }
        out.insert(motif_id.to_owned(), group_id.to_owned());
    }
    Ok(out)
}

fn read_reverse_dbd_best(path: &Path) -> Result<HashMap<String, String>> {
    let alignments = read_dbd_alignments(path)?;
    let mut best = HashMap::<String, (f64, String)>::new();
    for alignment in alignments {
        let donor = alignment.query;
        let local = alignment.subject;
        let replace = best
            .get(&donor)
            .is_none_or(|(score, _)| alignment.bitscore > *score);
        if replace {
            best.insert(donor, (alignment.bitscore, local));
        }
    }
    Ok(best
        .into_iter()
        .map(|(donor, (_, local))| (donor, local))
        .collect())
}

#[derive(Clone, Debug)]
struct DbdCandidate {
    donor_id: String,
    bitscore: f64,
    identity: f64,
    query_coverage: f64,
    subject_coverage: f64,
    alignment_length: f64,
    evalue: f64,
    query_start: f64,
    query_end: f64,
    subject_start: f64,
    subject_end: f64,
}

fn select_dbd_candidates(
    gene_family: &str,
    donor_dbds: &BTreeMap<String, DonorDbd>,
    reference_family: &HashMap<String, String>,
    alignments: &[DbdAlignment],
    thresholds: &DbdThresholdValues,
) -> (Vec<DbdCandidate>, Vec<&'static str>) {
    let mut any_alignment = false;
    let mut any_family = false;
    let mut any_below = false;
    let mut by_donor = BTreeMap::<String, DbdCandidate>::new();
    for alignment in alignments {
        any_alignment = true;
        let Some(donor) = donor_dbds.get(&alignment.subject) else {
            continue;
        };
        let donor_family = reference_family
            .get(&alignment.subject)
            .map(String::as_str)
            .unwrap_or(donor.family.as_str());
        if !family_matches(gene_family, donor_family) {
            continue;
        }
        any_family = true;
        if alignment.identity < thresholds.min_identity
            || alignment.query_coverage < thresholds.min_coverage
            || alignment.subject_coverage < thresholds.min_coverage
            || alignment.evalue > thresholds.max_evalue
        {
            any_below = true;
            continue;
        }
        let candidate = DbdCandidate {
            donor_id: alignment.subject.clone(),
            bitscore: alignment.bitscore,
            identity: alignment.identity,
            query_coverage: alignment.query_coverage,
            subject_coverage: alignment.subject_coverage,
            alignment_length: alignment.alignment_length,
            evalue: alignment.evalue,
            query_start: alignment.query_start,
            query_end: alignment.query_end,
            subject_start: alignment.subject_start,
            subject_end: alignment.subject_end,
        };
        let replace = by_donor
            .get(&candidate.donor_id)
            .is_none_or(|current| {
                candidate
                    .bitscore
                    .partial_cmp(&current.bitscore)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then(
                        candidate
                            .identity
                            .partial_cmp(&current.identity)
                            .unwrap_or(std::cmp::Ordering::Equal),
                    )
                    .is_gt()
            });
        if replace {
            by_donor.insert(candidate.donor_id.clone(), candidate);
        }
    }
    let mut candidates = by_donor.into_values().collect::<Vec<_>>();
    candidates.sort_by(|a, b| {
        b.bitscore
            .partial_cmp(&a.bitscore)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(
                b.identity
                    .partial_cmp(&a.identity)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
            .then_with(|| a.donor_id.cmp(&b.donor_id))
    });
    let mut reasons = Vec::new();
    if candidates.is_empty() {
        if !any_alignment {
            reasons.push("NO_BLAST_HIT");
        } else if !any_family {
            reasons.push("FAMILY_MISMATCH");
        } else if any_below {
            reasons.push("BELOW_DBD_THRESHOLD");
        } else {
            reasons.push("NO_BLAST_HIT");
        }
    }
    (candidates, reasons)
}

#[derive(Clone, Debug)]
struct DbdMappingRow {
    tf_id: String,
    tf_symbol: String,
    family: String,
    motif_id: String,
    reference_tf_id: String,
    mapping_method: String,
    domain_mode: String,
    extraction_family: String,
    extraction_models: String,
    domain_count: usize,
    mapping_confidence: String,
    mapping_status: String,
    reciprocal_dbd: String,
    specificity_groups: String,
    specificity_group_count: usize,
    dbd_identity: f64,
    dbd_bitscore: f64,
    dbd_evalue: f64,
    query_dbd_coverage: f64,
    donor_dbd_coverage: f64,
    alignment_length: f64,
    query_dbd_start: f64,
    query_dbd_end: f64,
    donor_dbd_start: f64,
    donor_dbd_end: f64,
    query_full_dbd_start: f64,
    query_full_dbd_end: f64,
    donor_full_dbd_start: f64,
    donor_full_dbd_end: f64,
    threshold_config_version: String,
    threshold_source: String,
    architecture_status: String,
    ambiguity_status: String,
    reason_code: String,
    reference_pack_version: String,
}

#[derive(Clone, Debug)]
struct DbdDiagnostic {
    tf_id: String,
    tf_symbol: String,
    family: String,
    status: String,
    mapping_status: String,
    primary_reason: &'static str,
    all_reasons: String,
    local_dbd_start: String,
    local_dbd_end: String,
    local_dbd_mode: String,
    query_models: String,
    family_donor_count: usize,
    eligible_donor_count: usize,
    best_donor: String,
    best_identity: String,
    best_query_coverage: String,
    best_donor_coverage: String,
    best_evalue: String,
    threshold_config_version: String,
    ambiguity_status: String,
    notes: String,
}

pub(super) fn run_map_dbd_motifs(args: &[String]) -> Result<()> {
    let opts = parse_options_with_flags(
        args,
        &["--allow-uncalibrated", "--require-reciprocal-dbd"],
    )?;
    for key in opts.keys() {
        if !matches!(
            key.as_str(),
            "--tf-genes"
                | "--proteins"
                | "--dbd-models"
                | "--reference-map"
                | "--donor-dbd"
                | "--donor-dbd-faa"
                | "--threshold-config"
                | "--out"
                | "--blastp"
                | "--dbd-alignments"
                | "--reverse-dbd-alignments"
                | "--specificity-groups"
                | "--diagnostics-out"
                | "--reference-pack-version"
                | "--allow-uncalibrated"
                | "--require-reciprocal-dbd"
        ) {
            bail!("unknown option: {key}");
        }
    }
    let get = |key: &str| {
        opts.get(key)
            .map(String::as_str)
            .with_context(|| format!("missing {key}"))
    };
    let out_value = get("--out")?;
    let out = Path::new(out_value);
    for input in [
        "--tf-genes",
        "--proteins",
        "--dbd-models",
        "--reference-map",
        "--donor-dbd",
        "--donor-dbd-faa",
        "--threshold-config",
        "--dbd-alignments",
        "--reverse-dbd-alignments",
        "--specificity-groups",
        "--diagnostics-out",
    ] {
        if opts.get(input).is_some_and(|path| same_path(path, out_value)) {
            bail!("{input} and --out must be different files");
        }
    }
    let allow_uncalibrated = opts
        .get("--allow-uncalibrated")
        .is_some_and(|value| !matches!(value.to_ascii_lowercase().as_str(), "false" | "0"));
    let thresholds = read_dbd_thresholds(Path::new(get("--threshold-config")?))?;
    if thresholds.is_uncalibrated() && !allow_uncalibrated {
        bail!(
            "threshold config {} is uncalibrated; pass --allow-uncalibrated true to run an explicitly exploratory mapping",
            thresholds.version
        );
    }
    let require_reciprocal = opts
        .get("--require-reciprocal-dbd")
        .is_some_and(|value| !matches!(value.to_ascii_lowercase().as_str(), "false" | "0"));
    if !thresholds.is_uncalibrated() && opts.get("--specificity-groups").is_none() {
        bail!("calibrated mapping requires --specificity-groups");
    }
    if !thresholds.is_uncalibrated() && !require_reciprocal {
        bail!("calibrated mapping requires --require-reciprocal-dbd true");
    }
    if require_reciprocal
        && opts.get("--dbd-alignments").is_some()
        && opts.get("--reverse-dbd-alignments").is_none()
    {
        bail!("--require-reciprocal-dbd needs --reverse-dbd-alignments with --dbd-alignments");
    }
    let specificity_groups = opts
        .get("--specificity-groups")
        .map(|path| read_specificity_groups(Path::new(path)))
        .transpose()?
        .unwrap_or_default();
    let genes = io::read_tf_genes_allow_empty(Path::new(get("--tf-genes")?))?;
    let reference = io::read_tf_map(Path::new(get("--reference-map")?))?;
    let model_set = read_dbd_model_set(Path::new(get("--dbd-models")?))?;
    let (local_dbds, _) = extract_dbds_from_catalogue(
        Path::new(get("--tf-genes")?),
        Path::new(get("--proteins")?),
        Path::new(get("--dbd-models")?),
        true,
    )?;
    let local_by_gene: HashMap<String, ExtractedDbd> = local_dbds
        .into_iter()
        .map(|dbd| (dbd.gene_id.clone(), dbd))
        .collect();
    let donor_dbds = read_donor_dbds(
        Path::new(get("--donor-dbd")?),
        Path::new(get("--donor-dbd-faa")?),
    )?;
    let reference_family: HashMap<String, String> = reference
        .iter()
        .map(|donor| (donor.tf_id.clone(), donor.family.clone()))
        .collect();
    for donor in donor_dbds.values() {
        if let Some(expected) = reference_family.get(&donor.donor_id) {
            if expected != &donor.family {
                bail!(
                    "donor DBD {} family {} does not match reference map family {}",
                    donor.donor_id,
                    donor.family,
                    expected
                );
            }
            let expected_models = model_set
                .families
                .get(&donor.family)
                .map(|architecture| architecture.models().join("+"))
                .unwrap_or_default();
            if !donor.models.is_empty() && !expected_models.is_empty() && donor.models != expected_models {
                bail!(
                    "donor DBD {} models {} do not match config models {}",
                    donor.donor_id,
                    donor.models,
                    expected_models
                );
            }
        }
    }
    let reference_donor_ids: BTreeSet<&str> =
        reference.iter().map(|donor| donor.tf_id.as_str()).collect();
    let donor_dbd_missing = reference_donor_ids
        .iter()
        .filter(|donor| !donor_dbds.contains_key(**donor))
        .count();
    let donor_dbd_extra = donor_dbds
        .keys()
        .filter(|donor| !reference_donor_ids.contains(donor.as_str()))
        .count();
    let (alignments, reverse_best) = if let Some(path) = opts.get("--dbd-alignments") {
        let forward = read_dbd_alignments(Path::new(path))?;
        let reverse = opts
            .get("--reverse-dbd-alignments")
            .map(|path| read_reverse_dbd_best(Path::new(path)))
            .transpose()?
            .unwrap_or_default();
        (forward, reverse)
    } else {
        let tmp = temp_output_path(out);
        fs::create_dir_all(&tmp)?;
        let local_faa = tmp.join("local-dbd.faa");
        let donor_faa = tmp.join("donor-dbd.faa");
        let table = tmp.join("dbd-alignments.tsv");
        let reverse_table = tmp.join("reverse-dbd-alignments.tsv");
        let result = (|| -> Result<(Vec<DbdAlignment>, HashMap<String, String>)> {
            let mut local_writer = BufWriter::new(File::create(&local_faa)?);
            for dbd in local_by_gene.values() {
                if genes
                    .get(&dbd.gene_id)
                    .is_some_and(|gene| gene.is_identified())
                {
                    writeln!(local_writer, ">{}", dbd.gene_id)?;
                    local_writer.write_all(&dbd.sequence)?;
                    writeln!(local_writer)?;
                }
            }
            local_writer.flush()?;
            let mut donor_writer = BufWriter::new(File::create(&donor_faa)?);
            for donor in donor_dbds.values() {
                writeln!(donor_writer, ">{}", donor.donor_id)?;
                donor_writer.write_all(&donor.sequence)?;
                writeln!(donor_writer)?;
            }
            donor_writer.flush()?;
            let blastp = opts.get("--blastp").map_or("blastp", String::as_str);
            let evalue = format!("{:.3e}", thresholds.max_evalue().max(1e-300));
            let local_query_count = local_by_gene
                .values()
                .filter(|dbd| {
                    genes
                        .get(&dbd.gene_id)
                        .is_some_and(|gene| gene.is_identified())
                })
                .count()
                .max(1);
            let run = |query: &Path, subject: &Path, table: &Path, max_targets: usize| -> Result<()> {
                let search = std::process::Command::new(blastp)
                    .arg("-query")
                    .arg(query)
                    .arg("-subject")
                    .arg(subject)
                    .args([
                        "-outfmt",
                        "6 qseqid sseqid bitscore pident length qlen slen evalue qstart qend sstart send",
                    ])
                    .args([
                        "-max_target_seqs",
                        &max_targets.to_string(),
                        "-evalue",
                        &evalue,
                    ])
                    .arg("-out")
                    .arg(table)
                    .output()
                    .with_context(|| format!("start {blastp}; install BLAST+ or pass --blastp"))?;
                if !search.status.success() {
                    bail!(
                        "blastp failed ({}): {}",
                        search.status,
                        String::from_utf8_lossy(&search.stderr).trim()
                    );
                }
                Ok(())
            };
            run(&local_faa, &donor_faa, &table, donor_dbds.len().max(1))?;
            let forward = read_dbd_alignments(&table)?;
            let reverse = if require_reciprocal {
                run(&donor_faa, &local_faa, &reverse_table, local_query_count)?;
                read_reverse_dbd_best(&reverse_table)?
            } else {
                HashMap::new()
            };
            Ok((forward, reverse))
        })();
        for path in [&local_faa, &donor_faa, &table, &reverse_table] {
            let _ = fs::remove_file(path);
        }
        let _ = fs::remove_dir(&tmp);
        result?
    };
    let mut alignments_by_query = BTreeMap::<&str, Vec<&DbdAlignment>>::new();
    for alignment in &alignments {
        alignments_by_query
            .entry(alignment.query.as_str())
            .or_default()
            .push(alignment);
    }
    let mut rows = Vec::new();
    let mut diagnostics = Vec::new();
    let mut gene_ids = genes.keys().collect::<Vec<_>>();
    gene_ids.sort();
    for gene_id in gene_ids {
        let gene = &genes[gene_id];
        if !gene.is_identified() {
            diagnostics.push(dbd_diagnostic(
                gene, "NOT_IDENTIFIED", &[], &local_by_gene, &donor_dbds,
                &thresholds, "no_supported_motif", &[],
            ));
            continue;
        }
        let Some(local) = local_by_gene.get(gene.gene_id.as_str()) else {
            diagnostics.push(dbd_diagnostic(
                gene, "NO_DBD_SEQUENCE", &["NO_DBD_SEQUENCE"], &local_by_gene, &donor_dbds,
                &thresholds, "no_supported_motif", &[],
            ));
            continue;
        };
        let Some(threshold_values) = thresholds.values_for(&gene.family) else {
            diagnostics.push(dbd_diagnostic(
                gene, "NO_CALIBRATED_THRESHOLD", &["NO_CALIBRATED_THRESHOLD"],
                &local_by_gene, &donor_dbds, &thresholds, "no_supported_motif", &[],
            ));
            continue;
        };
        let family_donor_count = donor_dbds
            .values()
            .filter(|donor| family_matches(&gene.family, &donor.family))
            .count();
        if family_donor_count == 0 {
            diagnostics.push(dbd_diagnostic(
                gene, "NO_DONOR_DBD", &["NO_DONOR_DBD"], &local_by_gene, &donor_dbds,
                &thresholds, "no_supported_motif", &[],
            ));
            continue;
        }
        let empty = Vec::new();
        let query_alignments = alignments_by_query
            .get(gene.gene_id.as_str())
            .unwrap_or(&empty);
        let (candidates, mut reasons) = select_dbd_candidates(
            &gene.family,
            &donor_dbds,
            &reference_family,
            &query_alignments.iter().map(|a| (*a).clone()).collect::<Vec<_>>(),
            threshold_values,
        );
        if candidates.is_empty() {
            if reasons.is_empty() {
                reasons.push("UNRESOLVED");
            }
            let primary = reasons[0];
            diagnostics.push(dbd_diagnostic(
                gene, primary, &reasons, &local_by_gene, &donor_dbds,
                &thresholds, "no_supported_motif", &[],
            ));
            continue;
        }
        let mut candidates = candidates;
        if require_reciprocal {
            candidates.retain(|candidate| {
                reverse_best
                    .get(&candidate.donor_id)
                    .is_some_and(|local| local == &gene.gene_id)
            });
            if candidates.is_empty() {
                diagnostics.push(dbd_diagnostic(
                    gene,
                    "NOT_RECIPROCAL_DBD_BEST",
                    &["NOT_RECIPROCAL_DBD_BEST"],
                    &local_by_gene,
                    &donor_dbds,
                    &thresholds,
                    "no_supported_motif",
                    &[],
                ));
                continue;
            }
        }
        let donor_ids = candidates
            .iter()
            .map(|candidate| candidate.donor_id.as_str())
            .collect::<BTreeSet<_>>();
        let motif_ids = candidates
            .iter()
            .flat_map(|candidate| {
                reference
                    .iter()
                    .filter(move |donor| donor.tf_id == candidate.donor_id)
                    .map(|donor| donor.motif_id.as_str())
            })
            .collect::<BTreeSet<_>>();
        let mut specificity_group_ids = BTreeSet::new();
        let mut missing_specificity_group = false;
        for candidate in &candidates {
            let mut found = false;
            for donor in reference
                .iter()
                .filter(|donor| donor.tf_id == candidate.donor_id)
            {
                if let Some(group) = specificity_groups.get(&donor.motif_id) {
                    specificity_group_ids.insert(group.clone());
                    found = true;
                }
            }
            if !found {
                missing_specificity_group = true;
            }
        }
        let ambiguity_status = match (donor_ids.len(), motif_ids.len()) {
            (1, 1) => "unique",
            (1, _) => "multiple_motifs",
            (_, 1) => "multiple_donors",
            _ => "multiple_donors_and_motifs",
        };
        let mapping_status = if thresholds.is_uncalibrated() {
            "dbd_exploratory"
        } else if specificity_group_ids.len() == 1 && !missing_specificity_group {
            "dbd_calibrated"
        } else {
            "dbd_ambiguous"
        };
        let reason_code = match mapping_status {
            "dbd_calibrated" => "DBD_CALIBRATED",
            "dbd_ambiguous" => "DBD_AMBIGUOUS",
            _ => {
                if ambiguity_status == "unique" {
                    "DBD_TRANSFER"
                } else {
                    "DBD_TRANSFER_AMBIGUOUS"
                }
            }
        };
        for candidate in &candidates {
            for donor in reference
                .iter()
                .filter(|donor| donor.tf_id == candidate.donor_id)
            {
                rows.push(DbdMappingRow {
                    tf_id: gene.gene_id.clone(),
                    tf_symbol: gene.tf_symbol.clone(),
                    family: gene.family.clone(),
                    motif_id: donor.motif_id.clone(),
                    reference_tf_id: donor.tf_id.clone(),
                    mapping_method: if thresholds.is_uncalibrated() {
                        "dbd_transfer_exploratory".to_owned()
                    } else {
                        "dbd_transfer_calibrated".to_owned()
                    },
                    domain_mode: local.mode.to_owned(),
                    extraction_family: local.family.clone(),
                    extraction_models: local.models.join("+"),
                    domain_count: local.domain_count,
                    mapping_confidence: if thresholds.is_uncalibrated() {
                        "exploratory_uncalibrated".to_owned()
                    } else {
                        "calibrated".to_owned()
                    },
                    mapping_status: mapping_status.to_owned(),
                    reciprocal_dbd: if require_reciprocal {
                        "yes".to_owned()
                    } else {
                        "not_required".to_owned()
                    },
                    specificity_groups: specificity_group_ids
                        .iter()
                        .cloned()
                        .collect::<Vec<_>>()
                        .join("|"),
                    specificity_group_count: specificity_group_ids.len(),
                    dbd_identity: candidate.identity,
                    dbd_bitscore: candidate.bitscore,
                    dbd_evalue: candidate.evalue,
                    query_dbd_coverage: candidate.query_coverage,
                    donor_dbd_coverage: candidate.subject_coverage,
                    alignment_length: candidate.alignment_length,
                    query_dbd_start: candidate.query_start,
                    query_dbd_end: candidate.query_end,
                    donor_dbd_start: candidate.subject_start,
                    donor_dbd_end: candidate.subject_end,
                    query_full_dbd_start: local.start as f64 + candidate.query_start - 1.0,
                    query_full_dbd_end: local.start as f64 + candidate.query_end - 1.0,
                    donor_full_dbd_start: donor_dbds
                        .get(&candidate.donor_id)
                        .map(|donor| donor.start as f64 + candidate.subject_start - 1.0)
                        .unwrap_or(0.0),
                    donor_full_dbd_end: donor_dbds
                        .get(&candidate.donor_id)
                        .map(|donor| donor.start as f64 + candidate.subject_end - 1.0)
                        .unwrap_or(0.0),
                    threshold_config_version: thresholds.version.clone(),
                    threshold_source: if threshold_values.source.is_empty() {
                        thresholds.source.clone()
                    } else {
                        threshold_values.source.clone()
                    },
                    architecture_status: format!(
                        "query={};donor={}",
                        local.mode,
                        donor_dbds
                            .get(&candidate.donor_id)
                            .map(|donor| donor.mode.as_str())
                            .unwrap_or("unknown")
                    ),
                    ambiguity_status: ambiguity_status.to_owned(),
                    reason_code: reason_code.to_owned(),
                    reference_pack_version: opts
                        .get("--reference-pack-version")
                        .cloned()
                        .unwrap_or_default(),
                });
            }
        }
        let primary_reason = match mapping_status {
            "dbd_calibrated" => "DBD_CALIBRATED",
            "dbd_ambiguous" => "DBD_AMBIGUOUS",
            _ => "DBD_EXPLORATORY",
        };
        diagnostics.push(dbd_diagnostic(
            gene,
            primary_reason,
            &[],
            &local_by_gene,
            &donor_dbds,
            &thresholds,
            mapping_status,
            &candidates,
        ));
    }
    rows.sort_by(|a, b| {
        a.tf_id
            .cmp(&b.tf_id)
            .then_with(|| a.reference_tf_id.cmp(&b.reference_tf_id))
            .then_with(|| a.motif_id.cmp(&b.motif_id))
    });
    let tmp = temp_output_path(out);
    let result = (|| -> Result<()> {
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&tmp)?;
        writer.write_record([
            "tf_id",
            "tf_symbol",
            "family",
            "motif_id",
            "reference_tf_id",
            "mapping_method",
            "domain_mode",
            "extraction_family",
            "extraction_models",
            "domain_count",
            "mapping_confidence",
            "mapping_status",
            "reciprocal_dbd",
            "specificity_groups",
            "specificity_group_count",
            "dbd_identity",
            "dbd_bitscore",
            "dbd_evalue",
            "query_dbd_coverage",
            "donor_dbd_coverage",
            "alignment_length",
            "query_dbd_start",
            "query_dbd_end",
            "donor_dbd_start",
            "donor_dbd_end",
            "query_full_dbd_start",
            "query_full_dbd_end",
            "donor_full_dbd_start",
            "donor_full_dbd_end",
            "threshold_config_version",
            "threshold_source",
            "architecture_status",
            "ambiguity_status",
            "reason_code",
            "reference_pack_version",
        ])?;
        for row in &rows {
            writer.write_record([
                row.tf_id.as_str(),
                row.tf_symbol.as_str(),
                row.family.as_str(),
                row.motif_id.as_str(),
                row.reference_tf_id.as_str(),
                row.mapping_method.as_str(),
                row.domain_mode.as_str(),
                row.extraction_family.as_str(),
                row.extraction_models.as_str(),
                &row.domain_count.to_string(),
                row.mapping_confidence.as_str(),
                row.mapping_status.as_str(),
                row.reciprocal_dbd.as_str(),
                row.specificity_groups.as_str(),
                &row.specificity_group_count.to_string(),
                &format!("{:.4}", row.dbd_identity),
                &format!("{:.1}", row.dbd_bitscore),
                &format!("{:.3e}", row.dbd_evalue),
                &format!("{:.4}", row.query_dbd_coverage),
                &format!("{:.4}", row.donor_dbd_coverage),
                &format!("{:.1}", row.alignment_length),
                &format!("{:.1}", row.query_dbd_start),
                &format!("{:.1}", row.query_dbd_end),
                &format!("{:.1}", row.donor_dbd_start),
                &format!("{:.1}", row.donor_dbd_end),
                &format!("{:.1}", row.query_full_dbd_start),
                &format!("{:.1}", row.query_full_dbd_end),
                &format!("{:.1}", row.donor_full_dbd_start),
                &format!("{:.1}", row.donor_full_dbd_end),
                row.threshold_config_version.as_str(),
                row.threshold_source.as_str(),
                row.architecture_status.as_str(),
                row.ambiguity_status.as_str(),
                row.reason_code.as_str(),
                row.reference_pack_version.as_str(),
            ])?;
        }
        writer.flush()?;
        drop(writer);
        publish_temp(&tmp, out)
    })();
    if result.is_err() {
        cleanup_output_temp(&tmp);
    }
    result?;
    if let Some(path) = opts.get("--diagnostics-out") {
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
        writer.write_record([
            "tf_id",
            "tf_symbol",
            "family",
            "status",
            "mapping_status",
            "primary_reason",
            "all_reasons",
            "local_dbd_start",
            "local_dbd_end",
            "local_dbd_mode",
            "query_models",
            "family_donor_count",
            "eligible_donor_count",
            "best_donor",
            "best_identity",
            "best_query_coverage",
            "best_donor_coverage",
            "best_evalue",
            "threshold_config_version",
            "ambiguity_status",
            "notes",
        ])?;
        for diagnostic in &diagnostics {
            writer.write_record([
                diagnostic.tf_id.as_str(),
                diagnostic.tf_symbol.as_str(),
                diagnostic.family.as_str(),
                diagnostic.status.as_str(),
                diagnostic.mapping_status.as_str(),
                diagnostic.primary_reason,
                diagnostic.all_reasons.as_str(),
                diagnostic.local_dbd_start.as_str(),
                diagnostic.local_dbd_end.as_str(),
                diagnostic.local_dbd_mode.as_str(),
                diagnostic.query_models.as_str(),
                &diagnostic.family_donor_count.to_string(),
                &diagnostic.eligible_donor_count.to_string(),
                diagnostic.best_donor.as_str(),
                diagnostic.best_identity.as_str(),
                diagnostic.best_query_coverage.as_str(),
                diagnostic.best_donor_coverage.as_str(),
                diagnostic.best_evalue.as_str(),
                diagnostic.threshold_config_version.as_str(),
                diagnostic.ambiguity_status.as_str(),
                diagnostic.notes.as_str(),
            ])?;
        }
        writer.flush()?;
    }
    println!(
        "wrote {} {} DBD-based TF–motif mapping row(s); {} donor DBDs missing, {} extra; threshold config {} ({})",
        rows.len(),
        if thresholds.is_uncalibrated() {
            "exploratory"
        } else {
            "calibrated"
        },
        donor_dbd_missing,
        donor_dbd_extra,
        thresholds.version,
        thresholds.calibration_status
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn dbd_diagnostic(
    gene: &io::TfGene,
    primary_reason: &'static str,
    reasons: &[&'static str],
    local_by_gene: &HashMap<String, ExtractedDbd>,
    donor_dbds: &BTreeMap<String, DonorDbd>,
    thresholds: &DbdThresholdConfig,
    mapping_status: &str,
    candidates: &[DbdCandidate],
) -> DbdDiagnostic {
    let local = local_by_gene.get(gene.gene_id.as_str());
    let family_donor_count = donor_dbds
        .values()
        .filter(|donor| family_matches(&gene.family, &donor.family))
        .count();
    let all_reasons = if reasons.is_empty() {
        primary_reason.to_owned()
    } else {
        reasons.join("|")
    };
    let best = candidates.first();
    DbdDiagnostic {
        tf_id: gene.gene_id.clone(),
        tf_symbol: gene.tf_symbol.clone(),
        family: gene.family.clone(),
        status: gene.status.clone(),
        mapping_status: mapping_status.to_owned(),
        primary_reason,
        all_reasons,
        local_dbd_start: local.map(|dbd| dbd.start.to_string()).unwrap_or_default(),
        local_dbd_end: local.map(|dbd| dbd.end.to_string()).unwrap_or_default(),
        local_dbd_mode: local.map(|dbd| dbd.mode.to_owned()).unwrap_or_default(),
        query_models: local
            .map(|dbd| dbd.models.join("+"))
            .unwrap_or_default(),
        family_donor_count,
        eligible_donor_count: candidates.len(),
        best_donor: best.map(|candidate| candidate.donor_id.clone()).unwrap_or_default(),
        best_identity: best
            .map(|candidate| format!("{:.4}", candidate.identity))
            .unwrap_or_default(),
        best_query_coverage: best
            .map(|candidate| format!("{:.4}", candidate.query_coverage))
            .unwrap_or_default(),
        best_donor_coverage: best
            .map(|candidate| format!("{:.4}", candidate.subject_coverage))
            .unwrap_or_default(),
        best_evalue: best
            .map(|candidate| format!("{:.3e}", candidate.evalue))
            .unwrap_or_default(),
        threshold_config_version: thresholds.version.clone(),
        ambiguity_status: mapping_status.to_owned(),
        notes: String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::{
        best_hits, choose_family_label, family_matches, select_architecture_domains,
        select_dbd_candidates, DbdAlignment, DbdThresholdValues, DomainArchitecture, DonorDbd,
        EvidenceDomain, EvidenceProtein,
    };
    use std::collections::{BTreeMap, HashMap};

    #[test]
    fn matches_any_family_label() {
        assert!(family_matches("Hap4;bZIP", "Hap4"));
        assert!(family_matches("Hap4;bZIP", "bZIP"));
        assert!(family_matches("Hap4; bZIP", "bZIP"));
        assert!(family_matches("Fungal_trans;zf-C2H2", "zf-C2H2"));
    }

    #[test]
    fn rejects_mismatched_or_empty_labels() {
        assert!(!family_matches("Hap4;bZIP", "C2H2"));
        assert!(!family_matches("", "bZIP"));
        assert!(!family_matches("Hap4;bZIP", ""));
    }

    #[test]
    fn chooses_family_by_domain_score_not_input_order() {
        let mut models = BTreeMap::new();
        models.insert(
            "Hap4".to_owned(),
            DomainArchitecture::Single("Hap4_Hap_bind".to_owned()),
        );
        models.insert(
            "bZIP".to_owned(),
            DomainArchitecture::Single("bZIP_1".to_owned()),
        );
        let evidence = vec![EvidenceProtein {
            protein_id: "p1".to_owned(),
            domains: vec![
                EvidenceDomain {
                    model_name: "Hap4_Hap_bind".to_owned(),
                    protein_start: 1,
                    protein_end: 17,
                    domain_score: 10.0,
                },
                EvidenceDomain {
                    model_name: "bZIP_1".to_owned(),
                    protein_start: 50,
                    protein_end: 80,
                    domain_score: 30.0,
                },
            ],
        }];
        assert_eq!(
            choose_family_label("Hap4;bZIP", &models, &evidence),
            Some("bZIP")
        );
        assert_eq!(
            choose_family_label("bZIP;Hap4", &models, &evidence),
            Some("bZIP")
        );
    }

    #[test]
    fn architecture_selection_accepts_non_input_order() {
        let pou = EvidenceDomain {
            model_name: "Pou".to_owned(),
            protein_start: 10,
            protein_end: 80,
            domain_score: 90.0,
        };
        let homeo = EvidenceDomain {
            model_name: "Homeodomain".to_owned(),
            protein_start: 100,
            protein_end: 160,
            domain_score: 60.0,
        };
        let hits = [&homeo, &pou];
        let selected = select_architecture_domains(
            &hits,
            &["Homeodomain".to_owned(), "Pou".to_owned()],
        )
        .unwrap();
        assert_eq!(selected.len(), 2);
        assert_eq!(selected[0].model_name, "Pou");
        assert_eq!(selected[1].model_name, "Homeodomain");
    }

    #[test]
    fn blast_parser_accepts_padded_numeric_fields() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_blast_padded_{}.tsv",
            std::process::id()
        ));
        std::fs::write(
            &path,
            "q1\ts1\t 9.2\t25.000\t4\t57\t83\t8.2\t52\t55\t72\t75\n",
        )
        .unwrap();
        let hits = best_hits(&path).unwrap();
        assert_eq!(hits.get("q1").unwrap().0, 9.2);
        let _ = std::fs::remove_file(path);
    }

    fn donor(id: &str, family: &str) -> DonorDbd {
        DonorDbd {
            donor_id: id.to_owned(),
            family: family.to_owned(),
            models: "PF00000".to_owned(),
            start: 1,
            mode: "full_architecture".to_owned(),
            sequence: b"AAAA".to_vec(),
        }
    }

    fn alignment(
        query: &str,
        subject: &str,
        bitscore: f64,
        identity: f64,
        query_coverage: f64,
        subject_coverage: f64,
    ) -> DbdAlignment {
        DbdAlignment {
            query: query.to_owned(),
            subject: subject.to_owned(),
            bitscore,
            identity,
            alignment_length: 30.0,
            evalue: 1e-10,
            query_start: 1.0,
            query_end: 30.0,
            subject_start: 1.0,
            subject_end: 30.0,
            query_coverage,
            subject_coverage,
        }
    }

    #[test]
    fn calibrated_dbd_model_status_is_accepted() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_calibrated_models_{}.json",
            std::process::id()
        ));
        std::fs::write(
            &path,
            "{\"version\":\"test\",\"families\":{\"A\":{\"models\":[\"PF00000\"],\"status\":\"calibrated\"}}}",
        )
        .unwrap();
        let models = super::read_dbd_model_set(&path).unwrap();
        assert!(models.families.contains_key("A"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn dbd_candidates_apply_family_and_threshold_filters() {
        let donors = BTreeMap::from([
            ("D1".to_owned(), donor("D1", "A")),
            ("D2".to_owned(), donor("D2", "A")),
            ("D3".to_owned(), donor("D3", "B")),
        ]);
        let families = HashMap::from([
            ("D1".to_owned(), "A".to_owned()),
            ("D2".to_owned(), "A".to_owned()),
            ("D3".to_owned(), "B".to_owned()),
        ]);
        let thresholds = DbdThresholdValues {
            min_identity: 60.0,
            min_coverage: 0.7,
            max_evalue: 1e-5,
            source: "test".to_owned(),
        };
        let hits = vec![
            alignment("q1", "D1", 100.0, 80.0, 0.95, 0.95),
            alignment("q1", "D2", 90.0, 70.0, 0.90, 0.90),
            alignment("q1", "D3", 120.0, 95.0, 0.95, 0.95),
        ];
        let (candidates, reasons) =
            select_dbd_candidates("A", &donors, &families, &hits, &thresholds);
        assert!(reasons.is_empty());
        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[0].donor_id, "D1");
        assert_eq!(candidates[1].donor_id, "D2");

        let below = vec![alignment("q2", "D1", 100.0, 40.0, 0.95, 0.95)];
        let (candidates, reasons) =
            select_dbd_candidates("A", &donors, &families, &below, &thresholds);
        assert!(candidates.is_empty());
        assert_eq!(reasons, ["BELOW_DBD_THRESHOLD"]);

        let mismatch = vec![alignment("q3", "D3", 100.0, 90.0, 0.95, 0.95)];
        let (candidates, reasons) =
            select_dbd_candidates("A", &donors, &families, &mismatch, &thresholds);
        assert!(candidates.is_empty());
        assert_eq!(reasons, ["FAMILY_MISMATCH"]);

        let (candidates, reasons) =
            select_dbd_candidates("A", &donors, &families, &[], &thresholds);
        assert!(candidates.is_empty());
        assert_eq!(reasons, ["NO_BLAST_HIT"]);
    }

    #[test]
    fn dbd_candidates_keep_best_alignment_per_donor() {
        let donors = BTreeMap::from([("D1".to_owned(), donor("D1", "A"))]);
        let families = HashMap::from([("D1".to_owned(), "A".to_owned())]);
        let thresholds = DbdThresholdValues {
            min_identity: 50.0,
            min_coverage: 0.7,
            max_evalue: 1e-5,
            source: "test".to_owned(),
        };
        let hits = vec![
            alignment("q1", "D1", 80.0, 60.0, 0.8, 0.8),
            alignment("q1", "D1", 120.0, 75.0, 0.9, 0.9),
        ];
        let (candidates, _) =
            select_dbd_candidates("A", &donors, &families, &hits, &thresholds);
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].bitscore, 120.0);
        assert_eq!(candidates[0].identity, 75.0);
    }
}
