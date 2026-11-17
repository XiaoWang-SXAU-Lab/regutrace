//! Traceable TF → motif → target → site outputs for the normal CLI flow.
//!
//! The files written here are deliberately small derived views. They do not
//! replace `evidence.tsv`; they make the provenance and ambiguity already
//! present in `tf_genes.tsv`, `tf_map.tsv` and the motif library visible
//! without manual joins.

use super::*;
use regutrace_core::pwm::Pwm;

pub(super) const TF_MOTIF_TRACE_FILE: &str = "tf_motif_trace.tsv";
pub(super) const MOTIF_AMBIGUITY_FILE: &str = "motif_ambiguity.tsv";
pub(super) const TARGET_CANDIDATE_TIES_FILE: &str = "target_candidate_ties.tsv";

pub(super) struct TraceInput<'a> {
    pub genes: &'a HashMap<String, io::TfGene>,
    pub mappings: &'a [io::TfMotif],
    /// `(tf_id, motif_id)` pairs that were actually scanned. Mappings outside
    /// this set are retained in the trace with an explicit scan_status.
    pub scanned_pairs: &'a BTreeSet<(String, String)>,
    pub reference_map: &'a [io::TfMotif],
    pub motifs: &'a [Pwm],
    pub edges: &'a [io::EdgeRecord],
    pub include_fallback: bool,
    pub pack_version: &'a str,
    pub motif_library: &'a str,
}

#[derive(Clone, Debug, Default)]
pub(super) struct TraceSummary {
    pub tf_motif_rows: usize,
    pub ambiguous_tf_motif_rows: usize,
    pub motifs_with_multiple_tfs: usize,
    pub tied_mapping_rows: usize,
    pub site_tie_groups: usize,
    pub site_tie_rows: usize,
}

#[derive(Clone, Debug, Default)]
struct DonorInfo {
    tf_symbol: String,
    family: String,
}

#[derive(Default)]
struct MappingAgg {
    rows: usize,
    scanned_rows: usize,
    excluded_fallback_rows: usize,
    tf_symbols: BTreeSet<String>,
    families: BTreeSet<String>,
    donors: BTreeMap<String, DonorInfo>,
    methods: BTreeSet<String>,
    domain_modes: BTreeSet<String>,
    extraction_families: BTreeSet<String>,
    extraction_models: BTreeSet<String>,
    domain_counts: BTreeSet<String>,
}

#[derive(Default, Clone)]
struct EvidenceAgg {
    rows: usize,
    targets: BTreeSet<String>,
    best_p: Option<f64>,
    best_score: Option<f64>,
    example_target_id: String,
    example_target_symbol: String,
    example_promoter: String,
    example_site: String,
}

impl EvidenceAgg {
    fn add(&mut self, edge: &io::EdgeRecord) {
        self.rows += 1;
        self.targets.insert(edge.target_id.clone());
        if let Some(score) = edge.score {
            self.best_score = Some(self.best_score.map_or(score, |best| best.max(score)));
        }
        let better_p = match (edge.p_value, self.best_p) {
            (Some(p), Some(best)) => p < best,
            (Some(_), None) => true,
            _ => false,
        };
        if better_p || self.example_target_id.is_empty() {
            if better_p {
                self.best_p = edge.p_value;
            }
            self.example_target_id = edge.target_id.clone();
            self.example_target_symbol = edge.target_symbol.clone();
            self.example_promoter = promoter_label(edge);
            self.example_site = site_label(edge);
        }
    }

    fn merge(&mut self, other: &EvidenceAgg) {
        self.rows += other.rows;
        self.targets.extend(other.targets.iter().cloned());
        if let Some(score) = other.best_score {
            self.best_score = Some(self.best_score.map_or(score, |best| best.max(score)));
        }
        let better_p = match (other.best_p, self.best_p) {
            (Some(p), Some(best)) => p < best,
            (Some(_), None) => true,
            _ => false,
        };
        if better_p || self.example_target_id.is_empty() {
            if better_p {
                self.best_p = other.best_p;
            }
            self.example_target_id = other.example_target_id.clone();
            self.example_target_symbol = other.example_target_symbol.clone();
            self.example_promoter = other.example_promoter.clone();
            self.example_site = other.example_site.clone();
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SiteKey {
    target_id: String,
    motif_id: String,
    promoter_start: Option<u64>,
    promoter_end: Option<u64>,
    hit_start: u64,
    hit_end: u64,
    strand: String,
}

#[derive(Default)]
struct SiteAgg {
    target_symbol: String,
    tf_ids: BTreeSet<String>,
    tf_symbols: BTreeSet<String>,
    promoter_start: Option<u64>,
    promoter_end: Option<u64>,
    score: Option<f64>,
    p_value: Option<f64>,
    evidence_level: String,
}

fn promoter_label(edge: &io::EdgeRecord) -> String {
    match (edge.promoter_start, edge.promoter_end) {
        (Some(start), Some(end)) => format!("{start}-{end}"),
        _ => String::new(),
    }
}

fn site_label(edge: &io::EdgeRecord) -> String {
    match (edge.hit_start, edge.hit_end) {
        (Some(start), Some(end)) => format!("{start}-{end}({})", edge.strand),
        _ => String::new(),
    }
}

fn join_set(set: &BTreeSet<String>) -> String {
    set.iter().map(String::as_str).collect::<Vec<_>>().join(";")
}

fn donor_ids(donors: &BTreeMap<String, DonorInfo>) -> String {
    donors.keys().map(String::as_str).collect::<Vec<_>>().join(";")
}

fn donor_symbols(donors: &BTreeMap<String, DonorInfo>) -> String {
    donors
        .iter()
        .map(|(id, donor)| {
            if donor.tf_symbol.is_empty() {
                id.clone()
            } else {
                donor.tf_symbol.clone()
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn donor_families(donors: &BTreeMap<String, DonorInfo>) -> String {
    donors
        .values()
        .map(|donor| donor.family.clone())
        .collect::<Vec<_>>()
        .join(";")
}

fn donor_json(donors: &BTreeMap<String, DonorInfo>) -> String {
    let records: Vec<serde_json::Value> = donors
        .iter()
        .map(|(id, donor)| {
            serde_json::json!({
                "tf_id": id,
                "tf_symbol": donor.tf_symbol.clone(),
                "tf_family": donor.family.clone(),
            })
        })
        .collect();
    serde_json::Value::Array(records).to_string()
}

fn symbols_for(ids: &BTreeSet<String>, genes: &HashMap<String, io::TfGene>) -> String {
    ids.iter()
        .map(|id| {
            genes
                .get(id)
                .map(|gene| gene.tf_symbol.clone())
                .filter(|symbol| !symbol.is_empty())
                .unwrap_or_else(|| id.clone())
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn opt_f64(value: Option<f64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn opt_u64(value: Option<u64>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn tsv_writer(path: &Path, headers: &[&str]) -> Result<csv::Writer<File>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("create {}", path.display()))?;
    writer
        .write_record(headers)
        .with_context(|| format!("write header {}", path.display()))?;
    Ok(writer)
}

fn write_fields(writer: &mut csv::Writer<File>, fields: &[String]) -> Result<()> {
    writer
        .write_record(fields.iter().map(String::as_str))
        .context("write trace TSV row")?;
    Ok(())
}

pub(super) fn write_trace_outputs(
    out_dir: &Path,
    input: &TraceInput<'_>,
) -> Result<TraceSummary> {
    let motif_names: HashMap<&str, &str> = input
        .motifs
        .iter()
        .map(|motif| (motif.id.as_str(), motif.name.as_str()))
        .collect();
    let mut reference_by_key: HashMap<(&str, &str), &io::TfMotif> = HashMap::new();
    let mut reference_by_id: HashMap<&str, &io::TfMotif> = HashMap::new();
    for reference in input.reference_map {
        reference_by_key.insert(
            (reference.tf_id.as_str(), reference.motif_id.as_str()),
            reference,
        );
        reference_by_id
            .entry(reference.tf_id.as_str())
            .or_insert(reference);
    }
    let reference_for = |mapping: &io::TfMotif| -> Option<&io::TfMotif> {
        reference_by_key
            .get(&(
                mapping.reference_tf_id.as_str(),
                mapping.motif_id.as_str(),
            ))
            .copied()
            .or_else(|| reference_by_id.get(mapping.reference_tf_id.as_str()).copied())
    };

    let mut mapping_aggs: BTreeMap<(String, String), MappingAgg> = BTreeMap::new();
    let mut motif_tfs: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut motif_scanned_tfs: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut motif_reference_ids: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut motif_methods: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut motif_tied_rows: BTreeMap<String, usize> = BTreeMap::new();
    let mut tied_mapping_rows = 0usize;

    for mapping in input.mappings {
        let key = (mapping.tf_id.clone(), mapping.motif_id.clone());
        let fallback = mapping
            .domain_mode
            .trim()
            .eq_ignore_ascii_case("single_domain_fallback");
        let motif_present = motif_names.contains_key(mapping.motif_id.as_str());
        let scanned = motif_present
            && if fallback {
                input.include_fallback
            } else {
                input.scanned_pairs.contains(&key)
            };
        let agg = mapping_aggs.entry(key).or_default();
        agg.rows += 1;
        if scanned {
            agg.scanned_rows += 1;
        } else if fallback {
            agg.excluded_fallback_rows += 1;
        }
        if !mapping.tf_symbol.is_empty() {
            agg.tf_symbols.insert(mapping.tf_symbol.clone());
        }
        if !mapping.family.is_empty() {
            agg.families.insert(mapping.family.clone());
        }
        if !mapping.reference_tf_id.is_empty() {
            let donor = agg
                .donors
                .entry(mapping.reference_tf_id.clone())
                .or_default();
            if let Some(reference) = reference_for(mapping) {
                if donor.tf_symbol.is_empty() && !reference.tf_symbol.is_empty() {
                    donor.tf_symbol = reference.tf_symbol.clone();
                }
                if donor.family.is_empty() && !reference.family.is_empty() {
                    donor.family = reference.family.clone();
                }
            }
        }
        if !mapping.mapping_method.is_empty() {
            agg.methods.insert(mapping.mapping_method.clone());
        }
        if !mapping.domain_mode.is_empty() {
            agg.domain_modes.insert(mapping.domain_mode.clone());
        }
        if !mapping.extraction_family.is_empty() {
            agg.extraction_families
                .insert(mapping.extraction_family.clone());
        }
        if !mapping.extraction_models.is_empty() {
            agg.extraction_models
                .insert(mapping.extraction_models.clone());
        }
        if !mapping.domain_count.is_empty() {
            agg.domain_counts.insert(mapping.domain_count.clone());
        }

        motif_tfs
            .entry(mapping.motif_id.clone())
            .or_default()
            .insert(mapping.tf_id.clone());
        if scanned {
            motif_scanned_tfs
                .entry(mapping.motif_id.clone())
                .or_default()
                .insert(mapping.tf_id.clone());
        }
        if !mapping.reference_tf_id.is_empty() {
            motif_reference_ids
                .entry(mapping.motif_id.clone())
                .or_default()
                .insert(mapping.reference_tf_id.clone());
        }
        if !mapping.mapping_method.is_empty() {
            motif_methods
                .entry(mapping.motif_id.clone())
                .or_default()
                .insert(mapping.mapping_method.clone());
        }
        if mapping
            .mapping_method
            .to_ascii_lowercase()
            .contains("tied")
        {
            *motif_tied_rows.entry(mapping.motif_id.clone()).or_insert(0) += 1;
            tied_mapping_rows += 1;
        }
    }

    let tie_motifs: BTreeSet<&str> = motif_scanned_tfs
        .iter()
        .filter(|(_, tfs)| tfs.len() > 1)
        .map(|(motif_id, _)| motif_id.as_str())
        .collect();
    let mut evidence_aggs: BTreeMap<(String, String), EvidenceAgg> = BTreeMap::new();
    let mut site_aggs: BTreeMap<SiteKey, SiteAgg> = BTreeMap::new();
    for edge in input.edges {
        if edge.motif_id.is_empty() {
            continue;
        }
        evidence_aggs
            .entry((edge.tf_id.clone(), edge.motif_id.clone()))
            .or_default()
            .add(edge);
        if !tie_motifs.contains(edge.motif_id.as_str()) {
            continue;
        }
        let (Some(hit_start), Some(hit_end)) = (edge.hit_start, edge.hit_end) else {
            continue;
        };
        let key = SiteKey {
            target_id: edge.target_id.clone(),
            motif_id: edge.motif_id.clone(),
            promoter_start: edge.promoter_start,
            promoter_end: edge.promoter_end,
            hit_start,
            hit_end,
            strand: edge.strand.clone(),
        };
        let site = site_aggs.entry(key).or_default();
        if site.target_symbol.is_empty() {
            site.target_symbol = if edge.target_symbol.is_empty() {
                edge.target_id.clone()
            } else {
                edge.target_symbol.clone()
            };
        }
        site.tf_ids.insert(edge.tf_id.clone());
        site.tf_symbols.insert(if edge.tf_symbol.is_empty() {
            edge.tf_id.clone()
        } else {
            edge.tf_symbol.clone()
        });
        if site.promoter_start.is_none() {
            site.promoter_start = edge.promoter_start;
        }
        if site.promoter_end.is_none() {
            site.promoter_end = edge.promoter_end;
        }
        if site.score.is_none() {
            site.score = edge.score;
        }
        if site.p_value.is_none() {
            site.p_value = edge.p_value;
        }
        if site.evidence_level.is_empty() {
            site.evidence_level = edge.evidence_level.clone();
        }
    }

    let mut motif_evidence: BTreeMap<String, EvidenceAgg> = BTreeMap::new();
    for ((_, motif_id), agg) in &evidence_aggs {
        motif_evidence
            .entry(motif_id.clone())
            .or_default()
            .merge(agg);
    }

    let mut summary = TraceSummary {
        tf_motif_rows: mapping_aggs.len(),
        tied_mapping_rows,
        ..TraceSummary::default()
    };

    let trace_path = out_dir.join(TF_MOTIF_TRACE_FILE);
    let mut trace_writer = tsv_writer(
        &trace_path,
        &[
            "tf_id",
            "tf_symbol",
            "tf_family",
            "identification_status",
            "tf_source",
            "motif_id",
            "motif_name",
            "motif_library",
            "reference_pack",
            "reference_tf_id",
            "reference_tf_symbol",
            "reference_tf_family",
            "reference_donors",
            "mapping_methods",
            "mapping_ambiguous",
            "mapping_rows",
            "excluded_fallback_rows",
            "scan_status",
            "domain_modes",
            "extraction_families",
            "extraction_models",
            "domain_counts",
            "shared_motif_tf_count",
            "shared_motif_scanned_tf_count",
            "shared_motif_tfs",
            "evidence_rows",
            "target_genes",
            "best_p_value",
            "best_score",
            "example_target_id",
            "example_target_symbol",
            "example_promoter",
            "example_site",
        ],
    )?;
    let empty_evidence = EvidenceAgg::default();
    for ((tf_id, motif_id), agg) in &mapping_aggs {
        let gene = input.genes.get(tf_id);
        let tf_symbol = gene
            .map(|gene| gene.tf_symbol.clone())
            .filter(|symbol| !symbol.is_empty())
            .or_else(|| agg.tf_symbols.iter().next().cloned())
            .unwrap_or_else(|| tf_id.clone());
        let tf_family = gene
            .map(|gene| gene.family.clone())
            .filter(|family| !family.is_empty())
            .unwrap_or_else(|| join_set(&agg.families));
        let identification_status = gene
            .map(|gene| gene.status.clone())
            .unwrap_or_else(|| "unknown".to_owned());
        let tf_source = gene.map(|gene| gene.source.clone()).unwrap_or_default();
        let mapping_ambiguous = agg.rows > 1
            || agg
                .methods
                .iter()
                .any(|method| method.to_ascii_lowercase().contains("tied"));
        if mapping_ambiguous {
            summary.ambiguous_tf_motif_rows += 1;
        }
        let scan_status = if agg.scanned_rows > 0 {
            "scanned"
        } else if !motif_names.contains_key(motif_id.as_str()) {
            "missing_motif"
        } else if agg.excluded_fallback_rows > 0 {
            "excluded_single_domain_fallback"
        } else {
            "not_scanned"
        };
        let shared_tfs = motif_tfs.get(motif_id).cloned().unwrap_or_default();
        let shared_symbols = symbols_for(&shared_tfs, input.genes);
        let evidence = evidence_aggs
            .get(&(tf_id.clone(), motif_id.clone()))
            .unwrap_or(&empty_evidence);
        let row = vec![
            tf_id.clone(),
            tf_symbol,
            tf_family,
            identification_status,
            tf_source,
            motif_id.clone(),
            motif_names.get(motif_id.as_str()).copied().unwrap_or("").to_owned(),
            input.motif_library.to_owned(),
            input.pack_version.to_owned(),
            donor_ids(&agg.donors),
            donor_symbols(&agg.donors),
            donor_families(&agg.donors),
            donor_json(&agg.donors),
            join_set(&agg.methods),
            mapping_ambiguous.to_string(),
            agg.rows.to_string(),
            agg.excluded_fallback_rows.to_string(),
            scan_status.to_owned(),
            join_set(&agg.domain_modes),
            join_set(&agg.extraction_families),
            join_set(&agg.extraction_models),
            join_set(&agg.domain_counts),
            shared_tfs.len().to_string(),
            motif_scanned_tfs
                .get(motif_id)
                .map(|tfs| tfs.len())
                .unwrap_or(0)
                .to_string(),
            shared_symbols,
            evidence.rows.to_string(),
            evidence.targets.len().to_string(),
            opt_f64(evidence.best_p),
            opt_f64(evidence.best_score),
            evidence.example_target_id.clone(),
            evidence.example_target_symbol.clone(),
            evidence.example_promoter.clone(),
            evidence.example_site.clone(),
        ];
        write_fields(&mut trace_writer, &row)?;
    }
    trace_writer
        .flush()
        .with_context(|| format!("finish {}", trace_path.display()))?;

    let ambiguity_path = out_dir.join(MOTIF_AMBIGUITY_FILE);
    let mut ambiguity_writer = tsv_writer(
        &ambiguity_path,
        &[
            "motif_id",
            "motif_name",
            "motif_library",
            "reference_pack",
            "tf_count",
            "scanned_tf_count",
            "reference_tf_count",
            "tied_mapping_rows",
            "tf_candidates",
            "scanned_tf_candidates",
            "reference_tf_candidates",
            "reference_donors",
            "mapping_methods",
            "evidence_rows",
            "target_genes",
            "example_target_id",
            "example_target_symbol",
            "example_site",
        ],
    )?;
    for (motif_id, tfs) in &motif_tfs {
        let scanned_tfs = motif_scanned_tfs.get(motif_id).cloned().unwrap_or_default();
        let reference_ids = motif_reference_ids
            .get(motif_id)
            .cloned()
            .unwrap_or_default();
        let tied_rows = motif_tied_rows.get(motif_id).copied().unwrap_or(0);
        let multiple_tfs = tfs.len() > 1;
        if !multiple_tfs && reference_ids.len() <= 1 && tied_rows == 0 {
            continue;
        }
        if multiple_tfs {
            summary.motifs_with_multiple_tfs += 1;
        }
        let donors: BTreeMap<String, DonorInfo> = reference_ids
            .iter()
            .map(|id| {
                let info = reference_by_id
                    .get(id.as_str())
                    .map(|reference| DonorInfo {
                        tf_symbol: reference.tf_symbol.clone(),
                        family: reference.family.clone(),
                    })
                    .unwrap_or_default();
                (id.clone(), info)
            })
            .collect();
        let evidence = motif_evidence.get(motif_id).unwrap_or(&empty_evidence);
        let row = vec![
            motif_id.clone(),
            motif_names.get(motif_id.as_str()).copied().unwrap_or("").to_owned(),
            input.motif_library.to_owned(),
            input.pack_version.to_owned(),
            tfs.len().to_string(),
            scanned_tfs.len().to_string(),
            reference_ids.len().to_string(),
            tied_rows.to_string(),
            symbols_for(tfs, input.genes),
            symbols_for(&scanned_tfs, input.genes),
            donor_symbols(&donors),
            donor_json(&donors),
            motif_methods
                .get(motif_id)
                .cloned()
                .map(|methods| join_set(&methods))
                .unwrap_or_default(),
            evidence.rows.to_string(),
            evidence.targets.len().to_string(),
            evidence.example_target_id.clone(),
            evidence.example_target_symbol.clone(),
            evidence.example_site.clone(),
        ];
        write_fields(&mut ambiguity_writer, &row)?;
    }
    ambiguity_writer
        .flush()
        .with_context(|| format!("finish {}", ambiguity_path.display()))?;

    let ties_path = out_dir.join(TARGET_CANDIDATE_TIES_FILE);
    let mut ties_writer = tsv_writer(
        &ties_path,
        &[
            "target_id",
            "target_symbol",
            "motif_id",
            "motif_name",
            "promoter_start",
            "promoter_end",
            "hit_start",
            "hit_end",
            "strand",
            "tf_count",
            "tf_ids",
            "tf_candidates",
            "score",
            "p_value",
            "evidence_level",
        ],
    )?;
    for (site, agg) in &site_aggs {
        if agg.tf_ids.len() < 2 {
            continue;
        }
        summary.site_tie_groups += 1;
        summary.site_tie_rows += agg.tf_ids.len();
        let row = vec![
            site.target_id.clone(),
            agg.target_symbol.clone(),
            site.motif_id.clone(),
            motif_names.get(site.motif_id.as_str()).copied().unwrap_or("").to_owned(),
            opt_u64(agg.promoter_start),
            opt_u64(agg.promoter_end),
            site.hit_start.to_string(),
            site.hit_end.to_string(),
            site.strand.clone(),
            agg.tf_ids.len().to_string(),
            join_set(&agg.tf_ids),
            join_set(&agg.tf_symbols),
            opt_f64(agg.score),
            opt_f64(agg.p_value),
            agg.evidence_level.clone(),
        ];
        write_fields(&mut ties_writer, &row)?;
    }
    ties_writer
        .flush()
        .with_context(|| format!("finish {}", ties_path.display()))?;

    Ok(summary)
}

/// Standalone entry point for tracing an existing scan/pipeline output.
pub(super) fn run(args: &[String]) -> Result<()> {
    let opts = parse_options_with_flags(args, &["--include-single-domain-fallback"])?;
    for key in opts.keys() {
        if !matches!(
            key.as_str(),
            "--tf-genes"
                | "--tf-map"
                | "--evidence"
                | "--motifs"
                | "--reference-map"
                | "--pack-version"
                | "--out-dir"
                | "--include-single-domain-fallback"
        ) {
            bail!("unknown option: {key}");
        }
    }
    let get = |key: &str| {
        opts.get(key)
            .map(String::as_str)
            .with_context(|| format!("missing {key}"))
    };
    let tf_genes_path = PathBuf::from(get("--tf-genes")?);
    let tf_map_path = PathBuf::from(get("--tf-map")?);
    let evidence_path = PathBuf::from(get("--evidence")?);
    let motifs_path = PathBuf::from(get("--motifs")?);
    let out_dir = PathBuf::from(get("--out-dir")?);
    let genes = io::read_tf_genes(&tf_genes_path)
        .with_context(|| format!("read {}", tf_genes_path.display()))?;
    let mappings = io::read_tf_map(&tf_map_path)
        .with_context(|| format!("read {}", tf_map_path.display()))?;
    let reference_map = match opts.get("--reference-map") {
        Some(path) => io::read_tf_map(Path::new(path))
            .with_context(|| format!("read {}", path))?,
        None => Vec::new(),
    };
    let motifs = pwm::read_jaspar(&motifs_path)
        .with_context(|| format!("read {}", motifs_path.display()))?;
    let edges = io::read_edges(&evidence_path)
        .with_context(|| format!("read {}", evidence_path.display()))?;
    let include_fallback = opts
        .get("--include-single-domain-fallback")
        .is_some_and(|value| matches!(value.to_ascii_lowercase().as_str(), "true" | "1"));
    let scanned_pairs: BTreeSet<(String, String)> = mappings
        .iter()
        .filter(|mapping| {
            include_fallback
                || !mapping
                    .domain_mode
                    .trim()
                    .eq_ignore_ascii_case("single_domain_fallback")
        })
        .map(|mapping| (mapping.tf_id.clone(), mapping.motif_id.clone()))
        .collect();
    let pack_version = opts.get("--pack-version").cloned().unwrap_or_default();
    let motif_library = motifs_path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    fs::create_dir_all(&out_dir)
        .with_context(|| format!("create output directory {}", out_dir.display()))?;
    let summary = write_trace_outputs(
        &out_dir,
        &TraceInput {
            genes: &genes,
            mappings: &mappings,
            scanned_pairs: &scanned_pairs,
            reference_map: &reference_map,
            motifs: &motifs,
            edges: &edges,
            include_fallback,
            pack_version: &pack_version,
            motif_library: &motif_library,
        },
    )?;
    println!(
        "trace: {} TF–motif rows; {} ambiguous TF–motif row(s); {} motif(s) shared by >1 TF; {} tied mapping row(s); {} site tie group(s) / {} candidate row(s)",
        summary.tf_motif_rows,
        summary.ambiguous_tf_motif_rows,
        summary.motifs_with_multiple_tfs,
        summary.tied_mapping_rows,
        summary.site_tie_groups,
        summary.site_tie_rows
    );
    println!("outputs: {}", out_dir.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gene(id: &str, symbol: &str, family: &str, status: &str) -> io::TfGene {
        io::TfGene {
            gene_id: id.to_owned(),
            tf_symbol: symbol.to_owned(),
            family: family.to_owned(),
            status: status.to_owned(),
            source: "test".to_owned(),
        }
    }

    fn mapping(
        tf_id: &str,
        tf_symbol: &str,
        motif_id: &str,
        reference_tf_id: &str,
        method: &str,
        domain_mode: &str,
    ) -> io::TfMotif {
        io::TfMotif {
            tf_id: tf_id.to_owned(),
            tf_symbol: tf_symbol.to_owned(),
            family: "F".to_owned(),
            motif_id: motif_id.to_owned(),
            reference_tf_id: reference_tf_id.to_owned(),
            mapping_method: method.to_owned(),
            domain_mode: domain_mode.to_owned(),
            extraction_family: "F".to_owned(),
            extraction_models: "MODEL".to_owned(),
            domain_count: "1".to_owned(),
        }
    }

    fn reference_with_family(
        tf_id: &str,
        symbol: &str,
        family: &str,
        motif_id: &str,
    ) -> io::TfMotif {
        io::TfMotif {
            tf_id: tf_id.to_owned(),
            tf_symbol: symbol.to_owned(),
            family: family.to_owned(),
            motif_id: motif_id.to_owned(),
            ..io::TfMotif::default()
        }
    }

    fn reference(tf_id: &str, symbol: &str, motif_id: &str) -> io::TfMotif {
        reference_with_family(tf_id, symbol, "F", motif_id)
    }

    fn motif(id: &str, name: &str) -> Pwm {
        Pwm {
            id: id.to_owned(),
            name: name.to_owned(),
            matrix: vec![[0.25; 4]; 2],
        }
    }

    fn edge_window(
        tf_id: &str,
        target_id: &str,
        motif_id: &str,
        hit_start: u64,
        hit_end: u64,
        p_value: f64,
        score: f64,
        promoter_start: u64,
        promoter_end: u64,
    ) -> io::EdgeRecord {
        io::EdgeRecord {
            tf_id: tf_id.to_owned(),
            tf_symbol: tf_id.to_owned(),
            target_id: target_id.to_owned(),
            target_symbol: target_id.to_owned(),
            evidence_type: "prediction_motif".to_owned(),
            evidence_level: "D".to_owned(),
            method: "jaspar_pwm".to_owned(),
            source: "promoter_scan".to_owned(),
            promoter_start: Some(promoter_start),
            promoter_end: Some(promoter_end),
            motif_id: motif_id.to_owned(),
            hit_start: Some(hit_start),
            hit_end: Some(hit_end),
            strand: "+".to_owned(),
            score: Some(score),
            p_value: Some(p_value),
            ..io::EdgeRecord::default()
        }
    }

    fn edge(
        tf_id: &str,
        target_id: &str,
        motif_id: &str,
        hit_start: u64,
        hit_end: u64,
        p_value: f64,
        score: f64,
    ) -> io::EdgeRecord {
        edge_window(
            tf_id,
            target_id,
            motif_id,
            hit_start,
            hit_end,
            p_value,
            score,
            100,
            200,
        )
    }

    fn read_tsv(path: &Path) -> Vec<HashMap<String, String>> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .unwrap();
        let headers = reader.headers().unwrap().clone();
        reader
            .records()
            .map(|record| {
                let record = record.unwrap();
                headers
                    .iter()
                    .zip(record.iter())
                    .map(|(header, value)| (header.to_owned(), value.to_owned()))
                    .collect()
            })
            .collect()
    }

    #[test]
    fn trace_tables_expose_provenance_and_ties() {
        let mut genes = HashMap::new();
        genes.insert("TF1".to_owned(), gene("TF1", "TF1", "F", "identified"));
        genes.insert("TF2".to_owned(), gene("TF2", "TF2", "F", "identified"));
        genes.insert("TF3".to_owned(), gene("TF3", "TF3", "F", "ambiguous"));

        let mappings = vec![
            mapping("TF1", "TF1", "MA1", "R1", "reciprocal_dbd_similarity", "full_architecture"),
            mapping("TF1", "TF1", "MA1", "R2", "reciprocal_dbd_similarity", "full_architecture"),
            mapping("TF2", "TF2", "MA1", "R1", "reciprocal_dbd_tied_candidate", "full_architecture"),
            mapping("TF3", "TF3", "MA1", "R1", "reciprocal_dbd_similarity", "single_domain_fallback"),
        ];
        let scanned_pairs = BTreeSet::from([
            ("TF1".to_owned(), "MA1".to_owned()),
            ("TF2".to_owned(), "MA1".to_owned()),
        ]);
        let reference_map = vec![reference("R1", "REF1", "MA1"), reference("R2", "REF2", "MA1")];
        let motifs = vec![motif("MA1", "MOTIF ONE")];
        let edges = vec![
            edge("TF1", "G1", "MA1", 10, 20, 1e-5, 10.0),
            edge("TF2", "G1", "MA1", 10, 20, 1e-5, 10.0),
            edge("TF1", "G2", "MA1", 30, 40, 1e-4, 8.0),
            edge("TF2", "G1", "MA1", 50, 60, 1e-6, 12.0),
            // Same target, motif, relative hit and strand, but a different
            // promoter window: this must remain a separate tied site.
            edge_window("TF1", "G1", "MA1", 10, 20, 1e-5, 10.0, 300, 400),
            edge_window("TF2", "G1", "MA1", 10, 20, 1e-5, 10.0, 300, 400),
        ];
        let dir = std::env::temp_dir().join(format!(
            "regutrace-trace-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        let summary = write_trace_outputs(
            &dir,
            &TraceInput {
                genes: &genes,
                mappings: &mappings,
                scanned_pairs: &scanned_pairs,
                reference_map: &reference_map,
                motifs: &motifs,
                edges: &edges,
                include_fallback: false,
                pack_version: "test-pack",
                motif_library: "motifs.txt",
            },
        )
        .unwrap();
        assert_eq!(summary.tf_motif_rows, 3);
        assert_eq!(summary.ambiguous_tf_motif_rows, 2);
        assert_eq!(summary.motifs_with_multiple_tfs, 1);
        assert_eq!(summary.tied_mapping_rows, 1);
        assert_eq!(summary.site_tie_groups, 2);
        assert_eq!(summary.site_tie_rows, 4);

        let trace = read_tsv(&dir.join(TF_MOTIF_TRACE_FILE));
        let tf1 = trace.iter().find(|row| row["tf_id"] == "TF1").unwrap();
        assert_eq!(tf1["motif_name"], "MOTIF ONE");
        assert_eq!(tf1["reference_tf_symbol"], "REF1;REF2");
        assert_eq!(tf1["mapping_ambiguous"], "true");
        assert_eq!(tf1["shared_motif_tf_count"], "3");
        assert_eq!(tf1["shared_motif_scanned_tf_count"], "2");
        assert_eq!(tf1["evidence_rows"], "3");
        assert_eq!(tf1["target_genes"], "2");
        assert_eq!(tf1["example_site"], "10-20(+)");
        let tf1_donors: serde_json::Value =
            serde_json::from_str(&tf1["reference_donors"]).unwrap();
        let tf1_donors = tf1_donors.as_array().unwrap();
        assert_eq!(tf1_donors.len(), 2);
        let donor = |id: &str| {
            tf1_donors
                .iter()
                .find(|donor| donor["tf_id"] == id)
                .unwrap()
        };
        assert_eq!(donor("R1")["tf_symbol"], "REF1");
        assert_eq!(donor("R1")["tf_family"], "F");
        assert_eq!(donor("R2")["tf_symbol"], "REF2");
        assert_eq!(donor("R2")["tf_family"], "F");
        let tf3 = trace.iter().find(|row| row["tf_id"] == "TF3").unwrap();
        assert_eq!(tf3["scan_status"], "excluded_single_domain_fallback");
        assert_eq!(tf3["evidence_rows"], "0");

        let ambiguity = read_tsv(&dir.join(MOTIF_AMBIGUITY_FILE));
        assert_eq!(ambiguity.len(), 1);
        assert_eq!(ambiguity[0]["motif_id"], "MA1");
        assert_eq!(ambiguity[0]["tf_count"], "3");
        assert_eq!(ambiguity[0]["scanned_tf_count"], "2");
        assert_eq!(ambiguity[0]["reference_tf_count"], "2");
        assert_eq!(ambiguity[0]["tied_mapping_rows"], "1");
        assert_eq!(ambiguity[0]["tf_candidates"], "TF1;TF2;TF3");
        assert_eq!(ambiguity[0]["scanned_tf_candidates"], "TF1;TF2");
        assert_eq!(ambiguity[0]["reference_tf_candidates"], "REF1;REF2");
        let ambiguity_donors: serde_json::Value =
            serde_json::from_str(&ambiguity[0]["reference_donors"]).unwrap();
        assert_eq!(ambiguity_donors.as_array().unwrap().len(), 2);

        let ties = read_tsv(&dir.join(TARGET_CANDIDATE_TIES_FILE));
        assert_eq!(ties.len(), 2);
        for tie in &ties {
            assert_eq!(tie["target_id"], "G1");
            assert_eq!(tie["hit_start"], "10");
            assert_eq!(tie["hit_end"], "20");
            assert_eq!(tie["tf_count"], "2");
            assert_eq!(tie["tf_candidates"], "TF1;TF2");
        }
        assert_eq!(ties[0]["promoter_start"], "100");
        assert_eq!(ties[0]["promoter_end"], "200");
        assert_eq!(ties[1]["promoter_start"], "300");
        assert_eq!(ties[1]["promoter_end"], "400");
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn donor_records_keep_id_symbol_family_pairs() {
        let mut genes = HashMap::new();
        genes.insert("TF1".to_owned(), gene("TF1", "TF1", "F", "identified"));
        let mappings = vec![
            mapping("TF1", "TF1", "MA1", "D2", "reciprocal_dbd_similarity", "full_architecture"),
            mapping("TF1", "TF1", "MA1", "D1", "reciprocal_dbd_similarity", "full_architecture"),
            mapping("TF1", "TF1", "MA1", "D3", "reciprocal_dbd_similarity", "full_architecture"),
        ];
        let scanned_pairs = BTreeSet::from([("TF1".to_owned(), "MA1".to_owned())]);
        let reference_map = vec![
            reference_with_family("D1", "BETA", "F2", "MA1"),
            reference_with_family("D2", "ALPHA", "F1", "MA1"),
            reference_with_family("D3", "", "", "MA1"),
        ];
        let motifs = vec![motif("MA1", "MOTIF ONE")];
        let dir = std::env::temp_dir().join(format!(
            "regutrace-trace-donors-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        let summary = write_trace_outputs(
            &dir,
            &TraceInput {
                genes: &genes,
                mappings: &mappings,
                scanned_pairs: &scanned_pairs,
                reference_map: &reference_map,
                motifs: &motifs,
                edges: &[],
                include_fallback: false,
                pack_version: "test-pack",
                motif_library: "motifs.txt",
            },
        )
        .unwrap();
        assert_eq!(summary.tf_motif_rows, 1);
        assert_eq!(summary.ambiguous_tf_motif_rows, 1);

        let trace = read_tsv(&dir.join(TF_MOTIF_TRACE_FILE));
        let tf = trace.iter().find(|row| row["tf_id"] == "TF1").unwrap();
        assert_eq!(tf["reference_tf_id"], "D1;D2;D3");
        assert_eq!(tf["reference_tf_symbol"], "BETA;ALPHA;D3");
        assert_eq!(tf["reference_tf_family"], "F2;F1;");
        let donors: serde_json::Value = serde_json::from_str(&tf["reference_donors"]).unwrap();
        let donors = donors.as_array().unwrap();
        assert_eq!(donors.len(), 3);
        let donor = |id: &str| donors.iter().find(|donor| donor["tf_id"] == id).unwrap();
        assert_eq!(donor("D1")["tf_symbol"], "BETA");
        assert_eq!(donor("D1")["tf_family"], "F2");
        assert_eq!(donor("D2")["tf_symbol"], "ALPHA");
        assert_eq!(donor("D2")["tf_family"], "F1");
        assert_eq!(donor("D3")["tf_symbol"], "");
        assert_eq!(donor("D3")["tf_family"], "");

        let ambiguity = read_tsv(&dir.join(MOTIF_AMBIGUITY_FILE));
        assert_eq!(ambiguity.len(), 1);
        assert_eq!(ambiguity[0]["reference_tf_count"], "3");
        assert_eq!(ambiguity[0]["reference_tf_candidates"], "BETA;ALPHA;D3");
        let ambiguity_donors: serde_json::Value =
            serde_json::from_str(&ambiguity[0]["reference_donors"]).unwrap();
        assert_eq!(ambiguity_donors.as_array().unwrap().len(), 3);
        let _ = fs::remove_dir_all(&dir);
    }
}
