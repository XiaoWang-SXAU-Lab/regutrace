use super::*;
use serde::Serialize;
use sha2::{Digest, Sha256};

/// One provenance record in a versioned motif-reference registry. This is the
/// Phase 1 schema; it intentionally keeps unknown fields empty instead of
/// inventing source, method or license metadata.
#[derive(Clone, Debug, Serialize)]
pub struct RegistryRecord {
    pub motif_id: String,
    pub donor_id: String,
    pub family: String,
    pub source: String,
    pub source_version: String,
    pub accession: String,
    pub evidence_type: String,
    pub method: String,
    pub species: String,
    pub tax_id: String,
    pub license: String,
    pub donor_sequence_id: String,
    pub dbd_start: Option<usize>,
    pub dbd_end: Option<usize>,
    pub dbd_models: String,
    pub checksum: String,
    pub qc_status: String,
}

#[derive(Clone, Debug, Serialize)]
struct SourceCount {
    source: String,
    source_version: String,
    records: usize,
    unique_motifs: usize,
    unique_donors: usize,
}

#[derive(Clone, Debug, Serialize)]
struct FamilyCount {
    family: String,
    records: usize,
    unique_motifs: usize,
    unique_donors: usize,
    missing_dbd: usize,
}

#[derive(Clone, Debug, Serialize)]
struct RegistrySummary {
    records: usize,
    unique_motifs: usize,
    unique_donors: usize,
    missing_dbd: usize,
    qc_status: BTreeMap<String, usize>,
    licenses: BTreeMap<String, usize>,
    evidence_types: BTreeMap<String, usize>,
    donor_species: BTreeMap<String, usize>,
    sources: Vec<SourceCount>,
    families: Vec<FamilyCount>,
    #[serde(skip_serializing_if = "Option::is_none")]
    library: Option<LibrarySummary>,
}

#[derive(Clone, Debug, Serialize)]
struct DuplicateGroup {
    checksum: String,
    length: usize,
    motif_ids: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
struct LibrarySummary {
    path: String,
    motifs: usize,
    registry_motifs: usize,
    registry_missing_from_library: usize,
    registry_missing_examples: Vec<String>,
    library_unmapped: usize,
    library_unmapped_examples: Vec<String>,
    duplicate_group_count: usize,
    duplicate_motif_count: usize,
    duplicate_groups: Vec<DuplicateGroup>,
}

#[derive(serde::Deserialize)]
struct PackManifest {
    #[serde(default)]
    pack_version: String,
    #[serde(default)]
    species: String,
    #[serde(default)]
    files: BTreeMap<String, String>,
    #[serde(default)]
    checksums: BTreeMap<String, String>,
}

/// SHA-256 over the uppercase residue string as returned by `fasta::read_fasta`.
/// Formatting (line wrapping, case) is therefore not part of the checksum.
fn sequence_checksum(sequence: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(sequence);
    format!("sha256:{:x}", hasher.finalize())
}

fn required(record: &csv::StringRecord, index: usize, label: &str, row: usize) -> Result<String> {
    let value = record.get(index).unwrap_or("").trim();
    if value.is_empty() {
        bail!("registry row {row}: empty {label}");
    }
    Ok(value.to_owned())
}

pub fn read_registry(path: &Path) -> Result<Vec<RegistryRecord>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("open registry {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| {
        headers
            .iter()
            .position(|header| header.trim().eq_ignore_ascii_case(name))
    };
    let motif = col("motif_id").context("registry needs motif_id")?;
    let donor = col("donor_id").context("registry needs donor_id")?;
    let family = col("family");
    let source = col("source");
    let source_version = col("source_version");
    let accession = col("accession");
    let evidence_type = col("evidence_type");
    let method = col("method");
    let species = col("species");
    let tax_id = col("tax_id");
    let license = col("license");
    let donor_sequence = col("donor_sequence_id");
    let dbd_start = col("dbd_start");
    let dbd_end = col("dbd_end");
    let dbd_models = col("dbd_models");
    let checksum = col("checksum");
    let qc_status = col("qc_status");
    let cell = |record: &csv::StringRecord, index: Option<usize>| {
        index
            .and_then(|i| record.get(i))
            .map(str::trim)
            .unwrap_or("")
            .to_owned()
    };
    let parse_pos = |value: &str, label: &str, row: usize| -> Result<Option<usize>> {
        if value.is_empty() {
            return Ok(None);
        }
        let parsed = value
            .parse::<usize>()
            .with_context(|| format!("registry row {row}: invalid {label} '{value}'"))?;
        Ok(Some(parsed))
    };
    let mut records = Vec::new();
    for (index, record) in reader.records().enumerate() {
        let record = record?;
        let row = index + 2;
        let motif_id = required(&record, motif, "motif_id", row)?;
        let donor_id = required(&record, donor, "donor_id", row)?;
        records.push(RegistryRecord {
            motif_id,
            donor_id,
            family: cell(&record, family),
            source: cell(&record, source),
            source_version: cell(&record, source_version),
            accession: cell(&record, accession),
            evidence_type: cell(&record, evidence_type),
            method: cell(&record, method),
            species: cell(&record, species),
            tax_id: cell(&record, tax_id),
            license: cell(&record, license),
            donor_sequence_id: cell(&record, donor_sequence),
            dbd_start: parse_pos(&cell(&record, dbd_start), "dbd_start", row)?,
            dbd_end: parse_pos(&cell(&record, dbd_end), "dbd_end", row)?,
            dbd_models: cell(&record, dbd_models),
            checksum: cell(&record, checksum),
            qc_status: cell(&record, qc_status),
        });
    }
    Ok(records)
}

#[derive(Clone, Debug, Default)]
struct RowProvenance {
    species: String,
    tax_id: String,
}

/// Read optional per-donor provenance columns from the legacy reference map.
/// Absent columns stay empty; nothing is inferred from the manifest except the
/// pack-wide fallback species.
fn read_legacy_provenance(
    path: &Path,
) -> Result<HashMap<(String, String), RowProvenance>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("open {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| {
        headers
            .iter()
            .position(|header| header.trim().eq_ignore_ascii_case(name))
    };
    let tf = col("tf_id").context("reference map needs tf_id")?;
    let motif = col("motif_id").context("reference map needs motif_id")?;
    let species = col("reference_species");
    let tax_id = col("tax_id");
    let mut out = HashMap::new();
    for record in reader.records() {
        let record = record?;
        let tf_id = record.get(tf).unwrap_or("").trim().to_owned();
        let motif_id = record.get(motif).unwrap_or("").trim().to_owned();
        if tf_id.is_empty() || motif_id.is_empty() {
            continue;
        }
        out.insert(
            (tf_id, motif_id),
            RowProvenance {
                species: species
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .trim()
                    .to_owned(),
                tax_id: tax_id
                    .and_then(|i| record.get(i))
                    .unwrap_or("")
                    .trim()
                    .to_owned(),
            },
        );
    }
    Ok(out)
}

/// Verify every file recorded in `manifest.checksums` before trusting the pack.
/// This is deliberately strict: a missing or altered asset makes the pack
/// unreproducible and must stop registry construction.
fn verify_pack_checksums(pack_dir: &Path, manifest: &PackManifest) -> Result<()> {
    for (relative, expected) in &manifest.checksums {
        let path = pack_dir.join(relative);
        let bytes = fs::read(&path)
            .with_context(|| format!("pack file {} is missing", path.display()))?;
        let actual = format!("{:x}", Sha256::digest(&bytes));
        if &actual != expected {
            bail!(
                "pack checksum mismatch for {}: manifest {}, file {}",
                relative,
                expected,
                actual
            );
        }
    }
    for (key, relative) in &manifest.files {
        let path = pack_dir.join(relative);
        if !path.is_file() {
            bail!("pack file {} ({}) is missing", path.display(), key);
        }
    }
    Ok(())
}

pub fn read_legacy_pack_registry(pack_dir: &Path) -> Result<Vec<RegistryRecord>> {
    let manifest_path = pack_dir.join("manifest.json");
    let manifest: PackManifest = serde_json::from_reader(
        File::open(&manifest_path)
            .with_context(|| format!("open {}", manifest_path.display()))?,
    )
    .with_context(|| format!("parse {}", manifest_path.display()))?;
    if manifest.pack_version.trim().is_empty() {
        bail!("legacy pack manifest has no pack_version");
    }
    verify_pack_checksums(pack_dir, &manifest)?;
    if let Some(relative) = manifest.files.get("registry") {
        let path = pack_dir.join(relative);
        let records = read_registry(&path)
            .with_context(|| format!("read pack registry {}", path.display()))?;
        if records.is_empty() {
            bail!("pack registry {} has no records", path.display());
        }
        return Ok(records);
    }
    let reference_map_rel = manifest
        .files
        .get("reference_map")
        .context("legacy pack manifest has no files.reference_map")?;
    let reference_map = pack_dir.join(reference_map_rel);
    let map = io::read_tf_map(&reference_map)
        .with_context(|| format!("read legacy reference map {}", reference_map.display()))?;
    let provenance = read_legacy_provenance(&reference_map)?;
    let donor_sequences = manifest
        .files
        .get("reference_proteins")
        .map(|rel| pack_dir.join(rel))
        .filter(|path| path.is_file())
        .map(|path| {
            fasta::read_fasta(&path)
                .with_context(|| format!("read legacy reference proteins {}", path.display()))
                .map(|fasta| fasta.seqs)
        })
        .transpose()?
        .unwrap_or_default();
    let mut records = Vec::with_capacity(map.len());
    for donor in map {
        let sequence_id = donor.tf_id.clone();
        if !donor_sequences.is_empty() && !donor_sequences.contains_key(&sequence_id) {
            bail!(
                "legacy pack reference-map donor {} is absent from reference proteins",
                donor.tf_id
            );
        }
        let row = provenance
            .get(&(donor.tf_id.clone(), donor.motif_id.clone()))
            .cloned()
            .unwrap_or_default();
        let species = if row.species.is_empty() {
            manifest.species.clone()
        } else {
            row.species
        };
        let checksum = donor_sequences
            .get(&sequence_id)
            .map(|sequence| sequence_checksum(sequence))
            .unwrap_or_default();
        records.push(RegistryRecord {
            motif_id: donor.motif_id,
            donor_id: donor.tf_id.clone(),
            family: donor.family,
            source: "legacy-pack".to_owned(),
            source_version: manifest.pack_version.clone(),
            accession: donor.reference_tf_id.clone(),
            evidence_type: "legacy_reference_map".to_owned(),
            method: donor.mapping_method,
            species,
            tax_id: row.tax_id,
            license: "unknown".to_owned(),
            donor_sequence_id: sequence_id,
            dbd_start: None,
            dbd_end: None,
            dbd_models: String::new(),
            checksum,
            qc_status: "legacy_unverified".to_owned(),
        });
    }
    Ok(records)
}

pub fn validate_registry(
    records: &[RegistryRecord],
    donor_sequences: &HashMap<String, Vec<u8>>,
) -> Result<()> {
    let mut exact_keys = BTreeSet::new();
    for record in records {
        if record.motif_id.trim().is_empty() || record.donor_id.trim().is_empty() {
            bail!("registry contains an empty motif_id or donor_id");
        }
        let key = format!(
            "{}\x00{}\x00{}\x00{}",
            record.motif_id, record.donor_id, record.source, record.source_version
        );
        if !exact_keys.insert(key) {
            bail!(
                "registry duplicate record: motif {} donor {} source {}",
                record.motif_id,
                record.donor_id,
                record.source
            );
        }
        if !record.donor_sequence_id.is_empty() && !donor_sequences.is_empty() {
            if !donor_sequences.contains_key(&record.donor_sequence_id) {
                bail!(
                    "registry donor sequence {} is absent from the supplied protein FASTA",
                    record.donor_sequence_id
                );
            }
        }
        if !record.checksum.is_empty() {
            if !record.checksum.starts_with("sha256:") {
                bail!(
                    "registry motif {} has unsupported checksum '{}' (expected sha256:<hex>)",
                    record.motif_id,
                    record.checksum
                );
            }
            if let Some(sequence) = donor_sequences.get(&record.donor_sequence_id) {
                let expected = sequence_checksum(sequence);
                if record.checksum != expected {
                    bail!(
                        "registry motif {} donor {} checksum mismatch: registry {}, FASTA {}",
                        record.motif_id,
                        record.donor_id,
                        record.checksum,
                        expected
                    );
                }
            }
        }
        match (record.dbd_start, record.dbd_end) {
            (Some(start), Some(end)) => {
                if start == 0 || end < start {
                    bail!(
                        "registry motif {} has invalid DBD coordinates {}-{}",
                        record.motif_id,
                        start,
                        end
                    );
                }
                if let Some(sequence) = donor_sequences.get(&record.donor_sequence_id) {
                    if end > sequence.len() {
                        bail!(
                            "registry motif {} DBD end {} exceeds donor sequence length {}",
                            record.motif_id,
                            end,
                            sequence.len()
                        );
                    }
                }
            }
            (None, None) => {}
            _ => bail!(
                "registry motif {} has only one DBD coordinate (start/end must be paired)",
                record.motif_id
            ),
        }
    }
    Ok(())
}

/// A registry motif that is absent from the packaged PWM library cannot be
/// scanned, so this is a hard pack-integrity failure rather than a warning.
fn validate_registry_library(records: &[RegistryRecord], library: &[pwm::Pwm]) -> Result<()> {
    let library_ids: BTreeSet<&str> = library.iter().map(|motif| motif.id.as_str()).collect();
    let missing = records
        .iter()
        .map(|record| record.motif_id.as_str())
        .filter(|motif| !library_ids.contains(motif))
        .collect::<BTreeSet<_>>();
    if !missing.is_empty() {
        let examples = missing.iter().take(10).copied().collect::<Vec<_>>();
        bail!(
            "registry references {} motif(s) absent from the PWM library; examples: {}",
            missing.len(),
            examples.join(", ")
        );
    }
    Ok(())
}

fn summarize_library(path: &Path, records: &[RegistryRecord]) -> Result<LibrarySummary> {
    let library = pwm::read_jaspar(path)
        .with_context(|| format!("read motif library {}", path.display()))?;
    validate_registry_library(records, &library)?;
    let library_ids: BTreeSet<&str> = library.iter().map(|motif| motif.id.as_str()).collect();
    let registry_ids: BTreeSet<&str> =
        records.iter().map(|record| record.motif_id.as_str()).collect();
    let registry_missing_examples = registry_ids
        .iter()
        .filter(|id| !library_ids.contains(*id))
        .take(20)
        .map(|id| (*id).to_owned())
        .collect::<Vec<_>>();
    let library_unmapped_examples = library_ids
        .iter()
        .filter(|id| !registry_ids.contains(*id))
        .take(20)
        .map(|id| (*id).to_owned())
        .collect::<Vec<_>>();
    let mut groups = BTreeMap::<String, (usize, Vec<String>)>::new();
    for motif in &library {
        let mut hasher = Sha256::new();
        hasher.update(b"regutrace-normalized-matrix-v1\0");
        hasher.update((motif.matrix.len() as u64).to_be_bytes());
        for column in &motif.matrix {
            for value in column {
                hasher.update(value.to_bits().to_be_bytes());
            }
        }
        let digest = format!("{:x}", hasher.finalize());
        let entry = groups.entry(digest).or_insert((motif.matrix.len(), Vec::new()));
        entry.1.push(motif.id.clone());
    }
    let mut duplicate_groups = groups
        .into_iter()
        .filter(|(_, (_, ids))| ids.len() > 1)
        .map(|(checksum, (length, mut motif_ids))| {
            motif_ids.sort();
            DuplicateGroup {
                checksum,
                length,
                motif_ids,
            }
        })
        .collect::<Vec<_>>();
    duplicate_groups.sort_by(|a, b| a.motif_ids.cmp(&b.motif_ids));
    Ok(LibrarySummary {
        path: path.display().to_string(),
        motifs: library.len(),
        registry_motifs: registry_ids.len(),
        registry_missing_from_library: registry_ids
            .iter()
            .filter(|id| !library_ids.contains(*id))
            .count(),
        registry_missing_examples,
        library_unmapped: library_ids
            .iter()
            .filter(|id| !registry_ids.contains(*id))
            .count(),
        library_unmapped_examples,
        duplicate_group_count: duplicate_groups.len(),
        duplicate_motif_count: duplicate_groups.iter().map(|g| g.motif_ids.len()).sum(),
        duplicate_groups,
    })
}

fn summarize(records: &[RegistryRecord]) -> RegistrySummary {
    let mut sources = BTreeMap::<(String, String), (usize, BTreeSet<&str>, BTreeSet<&str>)>::new();
    let mut families = BTreeMap::<String, (usize, BTreeSet<&str>, BTreeSet<&str>, usize)>::new();
    let mut qc = BTreeMap::<String, usize>::new();
    let mut licenses = BTreeMap::<String, usize>::new();
    let mut evidence_types = BTreeMap::<String, usize>::new();
    let mut species = BTreeMap::<String, usize>::new();
    let mut motifs = BTreeSet::new();
    let mut donors = BTreeSet::new();
    let mut missing_dbd = 0usize;
    for record in records {
        motifs.insert(record.motif_id.as_str());
        donors.insert(record.donor_id.as_str());
        *qc.entry(record.qc_status.clone()).or_insert(0) += 1;
        *licenses.entry(record.license.clone()).or_insert(0) += 1;
        *evidence_types
            .entry(record.evidence_type.clone())
            .or_insert(0) += 1;
        *species.entry(record.species.clone()).or_insert(0) += 1;
        if record.dbd_start.is_none() {
            missing_dbd += 1;
        }
        let source = sources
            .entry((record.source.clone(), record.source_version.clone()))
            .or_default();
        source.0 += 1;
        source.1.insert(record.motif_id.as_str());
        source.2.insert(record.donor_id.as_str());
        let family = families
            .entry(record.family.clone())
            .or_default();
        family.0 += 1;
        family.1.insert(record.motif_id.as_str());
        family.2.insert(record.donor_id.as_str());
        if record.dbd_start.is_none() {
            family.3 += 1;
        }
    }
    let source_rows = sources
        .into_iter()
        .map(|((source, source_version), (records, motifs, donors))| SourceCount {
            source,
            source_version,
            records,
            unique_motifs: motifs.len(),
            unique_donors: donors.len(),
        })
        .collect();
    let family_rows = families
        .into_iter()
        .map(|(family, (records, motifs, donors, missing_dbd))| FamilyCount {
            family,
            records,
            unique_motifs: motifs.len(),
            unique_donors: donors.len(),
            missing_dbd,
        })
        .collect();
    RegistrySummary {
        records: records.len(),
        unique_motifs: motifs.len(),
        unique_donors: donors.len(),
        missing_dbd,
        qc_status: qc,
        licenses,
        evidence_types,
        donor_species: species,
        sources: source_rows,
        families: family_rows,
        library: None,
    }
}

fn write_registry(path: &Path, records: &[RegistryRecord]) -> Result<()> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record([
        "motif_id", "donor_id", "family", "source", "source_version", "accession",
        "evidence_type", "method", "species", "tax_id", "license", "donor_sequence_id",
        "dbd_start", "dbd_end", "dbd_models", "checksum", "qc_status",
    ])?;
    for record in records {
        writer.write_record([
            record.motif_id.as_str(),
            record.donor_id.as_str(),
            record.family.as_str(),
            record.source.as_str(),
            record.source_version.as_str(),
            record.accession.as_str(),
            record.evidence_type.as_str(),
            record.method.as_str(),
            record.species.as_str(),
            record.tax_id.as_str(),
            record.license.as_str(),
            record.donor_sequence_id.as_str(),
            &record.dbd_start.map(|v| v.to_string()).unwrap_or_default(),
            &record.dbd_end.map(|v| v.to_string()).unwrap_or_default(),
            record.dbd_models.as_str(),
            record.checksum.as_str(),
            record.qc_status.as_str(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

pub(super) fn run(args: &[String]) -> Result<()> {
    let args = match args.first().map(String::as_str) {
        Some("inspect") => &args[1..],
        Some(other) => bail!("unknown registry command: {other} (expected inspect)"),
        None => bail!("registry requires a subcommand: inspect"),
    };
    let opts = parse_options(args)?;
    for key in opts.keys() {
        if !matches!(
            key.as_str(),
            "--reference-pack"
                | "--registry"
                | "--out"
                | "--summary-json"
                | "--motif-library"
                | "--donor-proteins"
        ) {
            bail!("unknown option: {key}");
        }
    }
    let out = opts.get("--out").context("missing --out")?;
    if opts.get("--reference-pack").is_none() && opts.get("--registry").is_none() {
        bail!("supply --reference-pack or --registry");
    }
    let records = if let Some(path) = opts.get("--registry") {
        read_registry(Path::new(path))?
    } else {
        let pack = opts.get("--reference-pack").expect("checked above");
        read_legacy_pack_registry(Path::new(pack))?
    };
    let donor_sequences = if let Some(path) = opts.get("--donor-proteins") {
        fasta::read_fasta(Path::new(path))
            .with_context(|| format!("read donor proteins {path}"))?
            .seqs
    } else {
        let manifest_path = opts
            .get("--reference-pack")
            .map(|pack| Path::new(pack).join("manifest.json"));
        manifest_path
            .as_deref()
            .and_then(|path| path.parent())
            .map(|pack| -> Result<HashMap<String, Vec<u8>>> {
                let data = fs::read(pack.join("manifest.json"))?;
                let manifest: PackManifest = serde_json::from_slice(&data)?;
                let path = manifest
                    .files
                    .get("reference_proteins")
                    .map(|rel| pack.join(rel));
                if let Some(path) = path.filter(|p| p.is_file()) {
                    Ok(fasta::read_fasta(&path)?.seqs)
                } else {
                    Ok(HashMap::new())
                }
            })
            .transpose()?
            .unwrap_or_default()
    };
    validate_registry(&records, &donor_sequences)?;
    let mut summary = summarize(&records);
    if let Some(path) = opts.get("--motif-library") {
        let library = summarize_library(Path::new(path), &records)?;
        summary.library = Some(library);
    }
    write_registry(Path::new(out), &records)?;
    if let Some(path) = opts.get("--summary-json") {
        fs::write(path, serde_json::to_vec_pretty(&summary)?)
            .with_context(|| format!("write registry summary {path}"))?;
    }
    println!(
        "registry: {} records, {} motifs, {} donors; {} missing DBD coordinates",
        summary.records, summary.unique_motifs, summary.unique_donors, summary.missing_dbd
    );
    if let Some(library) = &summary.library {
        println!(
            "motif library: {} motifs, {} unused by registry, {} missing from library, {} duplicate normalized-matrix groups",
            library.motifs,
            library.library_unmapped,
            library.registry_missing_from_library,
            library.duplicate_group_count
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_reads_and_validates_provenance() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-registry-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("registry.tsv");
        fs::write(
            &path,
            "motif_id\tdonor_id\tfamily\tsource\tsource_version\tlicense\tdonor_sequence_id\tdbd_start\tdbd_end\tqc_status\n\
MA1\tD1\tMYB\tJASPAR\t2026\tCC-BY-4.0\tD1\t10\t40\tcurated\n\
MA2\tD1\tMYB\tJASPAR\t2026\tCC-BY-4.0\tD1\t\t\treview\n",
        )
        .unwrap();
        let records = read_registry(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].dbd_start, Some(10));
        assert_eq!(records[1].dbd_end, None);
        let sequences = HashMap::from([("D1".to_owned(), b"A".repeat(50).to_vec())]);
        validate_registry(&records, &sequences).unwrap();
        let summary = summarize(&records);
        assert_eq!(summary.unique_motifs, 2);
        assert_eq!(summary.unique_donors, 1);
        assert_eq!(summary.missing_dbd, 1);
        assert_eq!(summary.sources.len(), 1);
        assert_eq!(summary.families[0].family, "MYB");
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn registry_rejects_unpaired_or_invalid_dbd_coordinates() {
        let base = RegistryRecord {
            motif_id: "MA1".to_owned(),
            donor_id: "D1".to_owned(),
            family: "MYB".to_owned(),
            source: "test".to_owned(),
            source_version: "1".to_owned(),
            accession: String::new(),
            evidence_type: "test".to_owned(),
            method: "test".to_owned(),
            species: "test".to_owned(),
            tax_id: "0".to_owned(),
            license: "test".to_owned(),
            donor_sequence_id: "D1".to_owned(),
            dbd_start: Some(10),
            dbd_end: None,
            dbd_models: String::new(),
            checksum: String::new(),
            qc_status: "test".to_owned(),
        };
        let sequences = HashMap::from([("D1".to_owned(), b"A".repeat(20).to_vec())]);
        assert!(validate_registry(&[base.clone()], &sequences).is_err());
        let mut bad = base;
        bad.dbd_end = Some(21);
        assert!(validate_registry(&[bad], &sequences).is_err());
    }

    #[test]
    fn registry_rejects_duplicate_records_and_missing_donor_sequences() {
        let base = RegistryRecord {
            motif_id: "MA1".to_owned(),
            donor_id: "D1".to_owned(),
            family: "MYB".to_owned(),
            source: "test".to_owned(),
            source_version: "1".to_owned(),
            accession: String::new(),
            evidence_type: "test".to_owned(),
            method: "test".to_owned(),
            species: "test".to_owned(),
            tax_id: String::new(),
            license: "test".to_owned(),
            donor_sequence_id: "D1".to_owned(),
            dbd_start: None,
            dbd_end: None,
            dbd_models: String::new(),
            checksum: String::new(),
            qc_status: "test".to_owned(),
        };
        let sequences = HashMap::from([("D1".to_owned(), b"AAAA".to_vec())]);
        validate_registry(&[base.clone()], &sequences).unwrap();
        let duplicate = [base.clone(), base.clone()];
        let error = validate_registry(&duplicate, &sequences).unwrap_err().to_string();
        assert!(error.contains("duplicate record"), "{error}");
        let mut missing = base;
        missing.donor_sequence_id = "D_MISSING".to_owned();
        let error = validate_registry(&[missing], &sequences).unwrap_err().to_string();
        assert!(error.contains("absent from the supplied protein FASTA"), "{error}");
    }

    #[test]
    fn legacy_pack_adapter_reads_per_donor_species() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-legacy-provenance-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            dir.join("manifest.json"),
            "{\"pack_version\":\"p1\",\"species\":\"Plants\",\
             \"files\":{\"reference_map\":\"reference-map.tsv\",\
             \"reference_proteins\":\"donor.faa\"}}",
        )
        .unwrap();
        fs::write(dir.join("donor.faa"), ">D1\nAAAA\n>D2\nCCCC\n").unwrap();
        fs::write(
            dir.join("reference-map.tsv"),
            "tf_id\ttf_symbol\tfamily\tmotif_id\treference_species\ttax_id\n\
D1\tS1\tMYB\tMA1\tArabidopsis thaliana\t3702\n\
D2\tS2\tMYB\tMA2\tZea mays\t4577\n",
        )
        .unwrap();
        let records = read_legacy_pack_registry(&dir).unwrap();
        assert_eq!(records[0].species, "Arabidopsis thaliana");
        assert_eq!(records[0].tax_id, "3702");
        assert_eq!(records[1].species, "Zea mays");
        assert_eq!(records[1].tax_id, "4577");
        assert_eq!(records[0].source, "legacy-pack");
        assert_eq!(records[0].source_version, "p1");
        assert_eq!(records[0].qc_status, "legacy_unverified");
        assert_eq!(records[0].checksum, sequence_checksum(b"AAAA"));
        assert_eq!(records[1].checksum, sequence_checksum(b"CCCC"));
        let sequences = HashMap::from([
            ("D1".to_owned(), b"AAAA".to_vec()),
            ("D2".to_owned(), b"CCCC".to_vec()),
        ]);
        validate_registry(&records, &sequences).unwrap();
        let mut bad = records.clone();
        bad[0].checksum = sequence_checksum(b"GGGG");
        let error = validate_registry(&bad, &sequences).unwrap_err().to_string();
        assert!(error.contains("checksum mismatch"), "{error}");
        fs::write(
            dir.join("manifest.json"),
            "{\"pack_version\":\"p1\",\"species\":\"Plants\",\
             \"files\":{\"reference_map\":\"reference-map.tsv\",\
             \"reference_proteins\":\"donor.faa\"},\
             \"checksums\":{\"reference-map.tsv\":\"deadbeef\"}}",
        )
        .unwrap();
        let error = read_legacy_pack_registry(&dir).unwrap_err().to_string();
        assert!(error.contains("pack checksum mismatch"), "{error}");
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn legacy_pack_uses_embedded_registry_when_present() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-embedded-registry-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            dir.join("manifest.json"),
            "{\"pack_version\":\"p2026\",\"files\":{\
             \"registry\":\"registry.tsv\",\
             \"reference_proteins\":\"proteins.faa\"}}",
        )
        .unwrap();
        fs::write(dir.join("proteins.faa"), ">D1\nAAAA\n").unwrap();
        fs::write(
            dir.join("registry.tsv"),
            "motif_id\tdonor_id\tfamily\tsource\tsource_version\tlicense\tdonor_sequence_id\tchecksum\tqc_status\n\
MA1\tD1\tMYB\tJASPAR\t2026\tCC-BY-4.0\tD1\tsha256:63c1dd951ffedf6f7fd968ad4efa39b8ed584f162f46e715114ee184f8de9201\tjaspar_metadata_verified\n",
        )
        .unwrap();
        let records = read_legacy_pack_registry(&dir).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].source, "JASPAR");
        assert_eq!(records[0].license, "CC-BY-4.0");
        assert_eq!(records[0].qc_status, "jaspar_metadata_verified");
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn motif_library_reports_unused_and_duplicate_matrices() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-library-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        let library = dir.join("motifs.txt");
        fs::write(
            &library,
            ">MA1 a\nA [ 1 0 ]\nC [ 0 1 ]\nG [ 0 0 ]\nT [ 0 0 ]\n\
>MA2 b\nA [ 1 0 ]\nC [ 0 1 ]\nG [ 0 0 ]\nT [ 0 0 ]\n\
>MA3 c\nA [ 0 1 ]\nC [ 0 1 ]\nG [ 1 0 ]\nT [ 1 0 ]\n",
        )
        .unwrap();
        let base = RegistryRecord {
            motif_id: "MA1".to_owned(),
            donor_id: "D1".to_owned(),
            family: "MYB".to_owned(),
            source: "test".to_owned(),
            source_version: "1".to_owned(),
            accession: String::new(),
            evidence_type: "test".to_owned(),
            method: "test".to_owned(),
            species: "test".to_owned(),
            tax_id: String::new(),
            license: "test".to_owned(),
            donor_sequence_id: "D1".to_owned(),
            dbd_start: None,
            dbd_end: None,
            dbd_models: String::new(),
            checksum: String::new(),
            qc_status: "test".to_owned(),
        };
        let records = vec![
            base.clone(),
            RegistryRecord {
                motif_id: "MA2".to_owned(),
                donor_id: "D2".to_owned(),
                ..base.clone()
            },
            RegistryRecord {
                motif_id: "MA3".to_owned(),
                donor_id: "D3".to_owned(),
                ..base.clone()
            },
        ];
        let summary = summarize_library(&library, &records).unwrap();
        assert_eq!(summary.motifs, 3);
        assert_eq!(summary.library_unmapped, 0);
        assert_eq!(summary.duplicate_group_count, 1);
        assert_eq!(summary.duplicate_motif_count, 2);
        assert_eq!(
            summary.duplicate_groups[0].motif_ids,
            vec!["MA1".to_owned(), "MA2".to_owned()]
        );
        let partial = &records[..2];
        let summary = summarize_library(&library, partial).unwrap();
        assert_eq!(summary.library_unmapped, 1);
        assert_eq!(
            summary.library_unmapped_examples,
            vec!["MA3".to_owned()]
        );
        let mut missing = records.clone();
        missing[0].motif_id = "MA_NOT_IN_LIBRARY".to_owned();
        assert!(summarize_library(&library, &missing).is_err());
        let _ = fs::remove_dir_all(dir);
    }
}
