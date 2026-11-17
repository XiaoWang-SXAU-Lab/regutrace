use super::*;

// BLAST outfmt 6: qseqid sseqid bitscore pident length qlen slen evalue.
type BestHits = HashMap<String, (f64, BTreeSet<String>, bool)>;

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

fn best_hits(path: &Path) -> Result<BestHits> {
    best_hits_by(path, None)
}

type AlignmentIndex<'a> = (&'a HashMap<String, (String, String)>, &'a HashMap<String, (String, String)>);

fn best_hits_by(path: &Path, index: Option<AlignmentIndex<'_>>) -> Result<BestHits> {
    let mut best = HashMap::<String, (f64, BTreeSet<String>, bool)>::new();
    let mut rows = ReaderBuilder::new().delimiter(b'\t').has_headers(false)
        .from_reader(fasta::open_reader(path)?);
    for row in rows.records() {
        let row = row?;
        let query = row.get(0).unwrap_or("").trim();
        let subject = row.get(1).unwrap_or("").trim();
        if !matches!(row.len(), 8 | 12) || query.is_empty() || subject.is_empty() {
            bail!("expected eight or twelve BLAST columns in {}", path.display());
        }
        let (query, subject) = if let Some((queries, subjects)) = index {
            let (Some((query_id, query_family)), Some((subject_id, subject_family))) =
                (queries.get(query), subjects.get(subject)) else { continue };
            if !family_matches(query_family, subject_family)
                && !family_matches(subject_family, query_family) {
                continue;
            }
            (query_id.as_str(), subject_id.as_str())
        } else {
            (query, subject)
        };
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
    Ok(best)
}

pub(super) fn run(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let inputs = [
        "--tf-genes", "--reference-map", "--forward", "--reverse", "--proteins",
        "--cds", "--reference-domains", "--dbd-models", "--reference-proteins",
    ];
    for key in opts.keys() {
        if !matches!(key.as_str(), "--out" | "--blastp") && !inputs.contains(&key.as_str()) {
            bail!("unknown option: {key}");
        }
    }
    let get = |key: &str| opts.get(key).map(String::as_str)
        .with_context(|| format!("missing {key}"));
    let out = get("--out")?;
    for input in inputs {
        if opts.get(input).is_some_and(|path| same_path(path, out)) {
            bail!("{input} and --out must be different files");
        }
    }
    let whole_protein = opts.contains_key("--reference-proteins");
    if whole_protein && ["--forward", "--reverse", "--reference-domains", "--dbd-models"]
        .iter().any(|key| opts.contains_key(*key)) {
        bail!("--reference-proteins requires automatic whole-protein alignment");
    }
    let genes = io::read_tf_genes_allow_empty(Path::new(get("--tf-genes")?))?;
    let reference = io::read_tf_map(Path::new(get("--reference-map")?))?;
    if reference.iter().any(|donor| donor.family.is_empty()) {
        bail!("reference map needs a nonempty family for every donor");
    }
    let (forward, reverse, local_info) =
        if opts.contains_key("--forward") || opts.contains_key("--reverse") {
            for option in ["--proteins", "--cds", "--reference-domains", "--dbd-models", "--reference-proteins", "--blastp"] {
                if opts.contains_key(option) {
                    bail!("{option} applies only to automatic alignment");
                }
            }
            (
                best_hits(Path::new(get("--forward")?))?,
                best_hits(Path::new(get("--reverse")?))?,
                HashMap::new(),
            )
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
    println!("wrote {} inferred TF–motif mappings", mapped.len());
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
) -> Result<(BestHits, BestHits, HashMap<String, LocalDomainInfo>)> {
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
            return Ok((BestHits::new(), BestHits::new(), HashMap::new()));
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
            Ok((
                best_hits_by(forward, Some((&local_index, &donor_index)))?,
                best_hits_by(reverse, Some((&donor_index, &local_index)))?,
                info,
            ))
        } else {
            Ok((best_hits(forward)?, best_hits(reverse)?, info))
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

#[derive(serde::Deserialize)]
#[serde(untagged)]
enum DomainArchitecture {
    Single(String),
    Ordered(Vec<String>),
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

fn write_local_domains(
    catalogue: &Path,
    sequences: &Path,
    models: &Path,
    output: &Path,
) -> Result<HashMap<String, LocalDomainInfo>> {
    let model_by_family: BTreeMap<String, DomainArchitecture> =
        serde_json::from_reader(fasta::open_reader(models)?)?;
    if model_by_family.is_empty() {
        bail!("DNA-binding domain model list is empty");
    }
    for (family, architecture) in &model_by_family {
        if architecture.models().is_empty()
            || architecture.models().iter().any(|name| name.trim().is_empty())
        {
            bail!("empty domain architecture for {family}");
        }
    }
    let proteins = fasta::read_fasta(sequences)?;
    let mut reader = ReaderBuilder::new().delimiter(b'\t')
        .from_reader(fasta::open_reader(catalogue)?);
    let headers = reader.headers()?.clone();
    let column = |name: &str| headers.iter().position(|col| col == name)
        .with_context(|| format!("TF catalogue needs {name} for automatic alignment"));
    let gene = column("gene_id")?;
    let status = column("identification_status")?;
    let family = column("family")?;
    let evidence = column("evidence_json")?;
    let mut writer = BufWriter::new(File::create(output)?);
    let mut info_by_gene = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let id = &row[gene];
        let evidence: Vec<EvidenceProtein> = serde_json::from_str(&row[evidence])
            .with_context(|| format!("read evidence for {id}"))?;
        let Some(family_label) = choose_family_label(&row[family], &model_by_family, &evidence)
        else { continue };
        let architecture = model_by_family
            .get(family_label)
            .expect("chosen family label has a configured architecture");
        let model_names = architecture.models().to_vec();
        let mut domains = BTreeMap::<Vec<u8>, (&'static str, usize, f64)>::new();
        for protein in evidence {
            let sequence = proteins.seqs.get(&protein.protein_id)
                .with_context(|| format!("{} is absent from input FASTA", protein.protein_id))?;
            let mut matching: Vec<_> = protein.domains.iter()
                .filter(|hit| model_names.contains(&hit.model_name))
                .collect();
            matching.sort_by_key(|hit| (hit.protein_start, hit.protein_end));
            let architecture = select_architecture_domains(&matching, &model_names);
            let (start, end, mode, domain_count, score) = if let Some(selected) = architecture {
                let first = selected.first().context("empty architecture")?;
                let last = selected.last().context("empty architecture")?;
                let score = selected.iter().map(|hit| hit.domain_score)
                    .fold(f64::NEG_INFINITY, f64::max);
                (
                    first.protein_start,
                    last.protein_end,
                    "full_architecture",
                    model_names.len(),
                    score,
                )
            } else {
                // Multi-domain TFs (C2H2 arrays, Myb repeats, two Zn_clus
                // domains, ...) do not satisfy every configured architecture.
                // Fall back to the best-scoring single domain so the TF still
                // reaches the BLAST competition instead of being silently
                // dropped. The output labels this as single-domain evidence.
                let Some(best) = matching.iter().max_by(|a, b| {
                    a.domain_score
                        .partial_cmp(&b.domain_score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                }) else {
                    continue;
                };
                (best.protein_start, best.protein_end, "single_domain_fallback", 1, best.domain_score)
            };
            let start = start.checked_sub(1)
                .with_context(|| format!("invalid domain start for {id}"))?;
            let slice = sequence.get(start..end)
                .filter(|slice| !slice.is_empty())
                .with_context(|| format!("invalid domain coordinates for {id}"))?;
            let entry = domains.entry(slice.to_vec())
                .or_insert((mode, domain_count, score));
            if score > entry.2 {
                *entry = (mode, domain_count, score);
            }
        }
        // Keep ambiguous genes in the competition when they have a configured
        // label, but only the final mapping step can accept identified genes.
        // Isoforms can carry distinct but valid DBD sequences; keep the
        // best-scoring representative instead of dropping the whole gene.
        if matches!(&row[status], "identified" | "ambiguous") {
            if let Some((sequence, (mode, domain_count, _))) = domains
                .into_iter()
                .max_by(|(seq_a, a), (seq_b, b)| {
                    a.2.partial_cmp(&b.2)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then_with(|| seq_a.len().cmp(&seq_b.len()))
                })
            {
                writeln!(writer, ">{id}")?;
                writer.write_all(&sequence)?;
                writeln!(writer)?;
                info_by_gene.insert(id.to_string(), LocalDomainInfo {
                    gene_id: id.to_string(),
                    family: family_label.to_string(),
                    models: model_names.clone(),
                    mode,
                    domain_count,
                });
            }
        }
    }
    writer.flush()?;
    Ok(info_by_gene)
}

#[cfg(test)]
mod tests {
    use super::{
        best_hits, choose_family_label, family_matches, select_architecture_domains,
        DomainArchitecture, EvidenceDomain, EvidenceProtein,
    };
    use std::collections::BTreeMap;

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
}
