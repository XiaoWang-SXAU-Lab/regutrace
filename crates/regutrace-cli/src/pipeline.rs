use super::*;
use super::trace;
use sha2::{Digest, Sha256};
use std::io::Read;

const REQUIRED_FILES: [&str; 6] = [
    "hmm",
    "rules",
    "reference_map",
    "reference_domains",
    "dbd_models",
    "motifs",
];

const OPTIONAL_MAPPING_FILES: [&str; 6] = [
    "donor_dbd",
    "donor_dbd_faa",
    "exploratory_thresholds",
    "calibrated_dbd_models",
    "calibrated_thresholds",
    "specificity_groups",
];

#[derive(Debug, serde::Deserialize)]
struct Manifest {
    pack_version: String,
    #[serde(default)]
    species: String,
    #[serde(default)]
    assembly: String,
    #[serde(default)]
    motif_scope: String,
    files: BTreeMap<String, String>,
    #[serde(default)]
    checksums: BTreeMap<String, String>,
    #[serde(default)]
    tools: BTreeMap<String, String>,
    #[serde(default)]
    defaults: Defaults,
}

#[derive(Clone, Debug, Default, serde::Deserialize)]
struct Defaults {
    promoter_len: Option<usize>,
    pvalue: Option<f64>,
    gene_biotype: Option<String>,
}

#[derive(Debug)]
struct ReferencePack {
    dir: PathBuf,
    manifest: Manifest,
    paths: BTreeMap<String, PathBuf>,
}

impl ReferencePack {
    fn path(&self, key: &str) -> &Path {
        self.paths
            .get(key)
            .map(PathBuf::as_path)
            .unwrap_or_else(|| Path::new(""))
    }

    fn tool(&self, key: &str) -> String {
        self.manifest
            .tools
            .get(key)
            .cloned()
            .unwrap_or_else(|| key.to_owned())
    }
}

fn resolve_pack_path(pack_dir: &Path, relative: &str, label: &str) -> Result<PathBuf> {
    let relative_path = Path::new(relative);
    if relative_path.is_absolute()
        || relative_path.components().any(|component| {
            matches!(
                component,
                std::path::Component::ParentDir
                    | std::path::Component::RootDir
                    | std::path::Component::Prefix(_)
            )
        })
    {
        bail!("reference pack {label} must be a relative path inside the pack: {relative}");
    }
    let joined = pack_dir.join(relative_path);
    let root = pack_dir
        .canonicalize()
        .with_context(|| format!("open reference pack directory {}", pack_dir.display()))?;
    let canonical = joined
        .canonicalize()
        .with_context(|| format!("open reference pack {label}: {}", joined.display()))?;
    if !canonical.starts_with(&root) {
        bail!(
            "reference pack {label} resolves outside the pack directory: {}",
            joined.display()
        );
    }
    Ok(joined)
}

fn sha256_file(path: &Path) -> Result<String> {
    let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let read = file
            .read(&mut buffer)
            .with_context(|| format!("hash {}", path.display()))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn load_pack(dir: &Path) -> Result<ReferencePack> {
    let manifest_path = dir.join("manifest.json");
    let manifest: Manifest = serde_json::from_reader(
        File::open(&manifest_path)
            .with_context(|| format!("open reference pack manifest {}", manifest_path.display()))?,
    )
    .with_context(|| format!("parse reference pack manifest {}", manifest_path.display()))?;
    if manifest.pack_version.trim().is_empty() {
        bail!("reference pack manifest has an empty pack_version");
    }
    if manifest.checksums.is_empty() {
        bail!("reference pack manifest has no checksums; refusing an unverified pack");
    }

    let mut paths = BTreeMap::new();
    for key in REQUIRED_FILES.into_iter().chain(manifest.files.contains_key("reference_proteins").then_some("reference_proteins")) {
        let relative = manifest
            .files
            .get(key)
            .with_context(|| format!("reference pack manifest is missing files.{key}"))?;
        let expected = manifest
            .checksums
            .get(relative)
            .with_context(|| format!("reference pack manifest is missing a checksum for {relative}"))?
            .trim()
            .trim_start_matches('*')
            .to_ascii_lowercase();
        let path = resolve_pack_path(dir, relative, key)?;
        let actual = sha256_file(&path)?;
        if actual != expected {
            bail!(
                "reference pack checksum mismatch for {relative}: expected {expected}, got {actual}"
            );
        }
        paths.insert(key.to_owned(), path);
    }
    for key in OPTIONAL_MAPPING_FILES {
        let Some(relative) = manifest.files.get(key) else {
            continue;
        };
        let expected = manifest
            .checksums
            .get(relative)
            .with_context(|| format!("reference pack manifest is missing a checksum for {relative}"))?
            .trim()
            .trim_start_matches('*')
            .to_ascii_lowercase();
        let path = resolve_pack_path(dir, relative, key)?;
        let actual = sha256_file(&path)?;
        if actual != expected {
            bail!(
                "reference pack checksum mismatch for {relative}: expected {expected}, got {actual}"
            );
        }
        paths.insert(key.to_owned(), path);
    }
    let pack = ReferencePack {
        dir: dir.to_owned(),
        manifest,
        paths,
    };
    // Parse the motif library while loading the pack so duplicate or malformed
    // records fail before TF identification. The scanner and the report must
    // both refer to exactly one matrix per motif ID.
    pwm::read_jaspar(pack.path("motifs")).with_context(|| {
        format!(
            "validate reference pack motifs {}",
            pack.path("motifs").display()
        )
    })?;
    Ok(pack)
}

fn gene_biotype_status(
    requested: &Option<String>,
    from_pack: bool,
    scan_report: &serde_json::Value,
) -> serde_json::Value {
    let diagnostics = scan_report
        .get("diagnostics")
        .and_then(serde_json::Value::as_array)
        .map(Vec::as_slice)
        .unwrap_or(&[]);
    let has_code = |code: &str| {
        diagnostics.iter().any(|diagnostic| {
            diagnostic.get("code").and_then(serde_json::Value::as_str) == Some(code)
        })
    };
    let applied = requested.is_some() && has_code("GENE_BIOTYPE_FILTER");
    let skipped = has_code("GENE_BIOTYPE_FILTER_SKIPPED");
    let status = if requested.is_none() {
        "not_requested"
    } else if applied {
        "applied"
    } else if skipped {
        "skipped"
    } else {
        "unknown"
    };
    serde_json::json!({
        "requested": requested.as_deref(),
        "source": if requested.is_none() {
            "none"
        } else if from_pack {
            "pack_default"
        } else {
            "user"
        },
        "applied": applied,
        "skipped": skipped,
        "status": status,
    })
}

fn bool_option(opts: &HashMap<String, String>, key: &str) -> bool {
    opts.get(key).is_some_and(|value| {
        matches!(value.to_ascii_lowercase().as_str(), "true" | "1")
    })
}

fn read_tsv_rows(path: &Path) -> Result<Vec<BTreeMap<String, String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("open {}", path.display()))?;
    let headers = reader.headers()?.clone();
    let mut rows = Vec::new();
    for record in reader.records() {
        let record = record?;
        let mut row = BTreeMap::new();
        for (index, header) in headers.iter().enumerate() {
            row.insert(header.to_owned(), record.get(index).unwrap_or("").to_owned());
        }
        rows.push(row);
    }
    Ok(rows)
}

fn scan_summary_context(report: &serde_json::Value) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    if let Some(diagnostics) = report.get("diagnostics").and_then(|value| value.as_array()) {
        for diagnostic in diagnostics {
            if diagnostic.get("code").and_then(|value| value.as_str()) == Some("SCAN_SUMMARY") {
                if let Some(context) = diagnostic.get("context").and_then(|value| value.as_object())
                {
                    for (key, value) in context {
                        if let Some(value) = value.as_str() {
                            out.insert(key.clone(), value.to_string());
                        }
                    }
                }
            }
        }
    }
    out
}

fn parse_background_context(value: &str) -> Option<[f64; 4]> {
    let values: Vec<f64> = value
        .split(',')
        .map(|part| part.trim().parse::<f64>().ok())
        .collect::<Option<Vec<_>>>()?;
    if values.len() != 4
        || values
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return None;
    }
    Some([values[0], values[1], values[2], values[3]])
}

fn output_paths(out_dir: &Path) -> BTreeMap<&'static str, PathBuf> {
    [
        ("tf_genes", "tf_genes.tsv"),
        ("tf_map", "tf_map.tsv"),
        ("mapping_diagnostics", "mapping-diagnostics.tsv"),
        ("mapping_diagnostics_family", "mapping-diagnostics.family.tsv"),
        ("evidence", "evidence.tsv"),
        ("promoters", "promoters.fa"),
        ("progress", "progress.json"),
        ("tf_motif_trace", trace::TF_MOTIF_TRACE_FILE),
        ("motif_ambiguity", trace::MOTIF_AMBIGUITY_FILE),
        ("target_candidate_ties", trace::TARGET_CANDIDATE_TIES_FILE),
        ("database", "evidence.db"),
        ("scan_report", "scan-report.json"),
        ("sites", "motif-site-hits.tsv"),
        ("motif_testability", "motif-testability.tsv"),
        ("promoter_summary", "promoter-summary.tsv"),
        ("scan_manifest", "scan-manifest.json"),
        ("pipeline_report", "pipeline-report.json"),
    ]
    .into_iter()
    .map(|(key, name)| (key, out_dir.join(name)))
    .collect()
}


// Stage boundaries report actual work; elapsed-time percentages are not estimated.
fn report_pipeline_stage(path: &Path, stage: &str, message: &str) -> Result<()> {
    let progress = serde_json::json!({"stage": stage, "progress": 0.0, "message": message});
    fs::write(path, serde_json::to_vec(&progress)?).context("write pipeline stage")
}

fn optional_pack_path<'a>(pack: &'a ReferencePack, key: &str) -> Result<&'a Path> {
    let path = pack.path(key);
    if path.as_os_str().is_empty() {
        bail!(
            "reference pack {} does not provide {key}; extended mapping is unavailable for this pack",
            pack.manifest.pack_version
        );
    }
    Ok(path)
}

fn run_dbd_mapping(
    pack: &ReferencePack,
    genes: &Path,
    sequence: &Path,
    is_cds: bool,
    blastp: &Path,
    models: &Path,
    threshold: &Path,
    out: &Path,
    allow_uncalibrated: bool,
    reciprocal_required: bool,
) -> Result<()> {
    let donor = optional_pack_path(pack, "donor_dbd")?;
    let donor_faa = optional_pack_path(pack, "donor_dbd_faa")?;
    let mut args = vec![
        "--tf-genes".to_owned(),
        genes.to_string_lossy().into_owned(),
        "--dbd-models".to_owned(),
        models.to_string_lossy().into_owned(),
        "--reference-map".to_owned(),
        pack.path("reference_map").to_string_lossy().into_owned(),
        "--donor-dbd".to_owned(),
        donor.to_string_lossy().into_owned(),
        "--donor-dbd-faa".to_owned(),
        donor_faa.to_string_lossy().into_owned(),
        "--threshold-config".to_owned(),
        threshold.to_string_lossy().into_owned(),
        "--blastp".to_owned(),
        blastp.to_string_lossy().into_owned(),
        "--reference-pack-version".to_owned(),
        pack.manifest.pack_version.clone(),
        "--out".to_owned(),
        out.to_string_lossy().into_owned(),
    ];
    if is_cds {
        args.extend(["--cds".to_owned(), sequence.to_string_lossy().into_owned()]);
    } else {
        args.extend(["--proteins".to_owned(), sequence.to_string_lossy().into_owned()]);
    }
    if let Some(specificity) = pack.paths.get("specificity_groups") {
        args.extend([
            "--specificity-groups".to_owned(),
            specificity.to_string_lossy().into_owned(),
        ]);
    }
    if allow_uncalibrated {
        args.extend(["--allow-uncalibrated".to_owned(), "true".to_owned()]);
    }
    if reciprocal_required {
        args.extend(["--require-reciprocal-dbd".to_owned(), "true".to_owned()]);
    }
    motif_transfer::run_map_dbd_motifs(&args).context("DBD motif mapping failed")
}

fn merge_tf_maps(strict: &Path, extra: &[PathBuf], out: &Path) -> Result<usize> {
    let mut headers: Vec<String> = Vec::new();
    let mut records: Vec<(String, BTreeMap<String, String>)> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for (layer_index, path) in std::iter::once(strict)
        .chain(extra.iter().map(PathBuf::as_path))
        .enumerate()
    {
        if !path.exists() {
            continue;
        }
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .flexible(true)
            .from_path(path)
            .with_context(|| format!("open {}", path.display()))?;
        let file_headers = reader.headers()?.clone();
        for header in file_headers.iter() {
            if !headers.iter().any(|item| item == header) {
                headers.push(header.to_owned());
            }
        }
        for record in reader.records() {
            let record = record?;
            let mut row = BTreeMap::new();
            for (index, header) in file_headers.iter().enumerate() {
                row.insert(header.to_owned(), record.get(index).unwrap_or("").to_owned());
            }
            if layer_index == 0
                && row
                    .get("mapping_status")
                    .is_none_or(|value| value.trim().is_empty())
            {
                row.insert("mapping_status".to_owned(), "strict_rbh".to_owned());
            }
            let key = format!(
                "{}|{}|{}",
                row.get("tf_id").cloned().unwrap_or_default(),
                row.get("motif_id").cloned().unwrap_or_default(),
                row.get("reference_tf_id").cloned().unwrap_or_default()
            );
            if seen.insert(key.clone()) {
                records.push((key, row));
            }
        }
    }
    if headers.is_empty() {
        bail!("TF–motif mapping produced no header");
    }
    let tmp = temp_output_path(out);
    {
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&tmp)
            .with_context(|| format!("write {}", tmp.display()))?;
        writer.write_record(&headers)?;
        for (_, row) in &records {
            writer.write_record(headers.iter().map(|header| {
                row.get(header).map(String::as_str).unwrap_or("")
            }))?;
        }
        writer.flush()?;
    }
    publish_temp(&tmp, out)?;
    Ok(records.len())
}

const PIPELINE_FLAGS: &[&str] = &[
    "--include-single-domain-fallback",
    "--reuse-exact-pwm",
    "--phase3-artifacts",
    "--sites-only",
];

pub(super) fn run(args: &[String]) -> Result<()> {
    let opts = parse_options_with_flags(args, PIPELINE_FLAGS)?;
    for key in opts.keys() {
        if !matches!(
            key.as_str(),
            "--reference-pack"
                | "--genome"
                | "--gff"
                | "--proteins"
                | "--cds"
                | "--annotation"
                | "--out-dir"
                | "--include-single-domain-fallback"
                | "--reuse-exact-pwm"
                | "--phase3-artifacts"
                | "--promoter-len"
                | "--downstream"
                | "--pvalue"
                | "--background"
                | "--background-mode"
                | "--scan-cache"
                | "--gene-biotype"
                | "--threads"
                | "--hmmsearch"
                | "--blastp"
                | "--mapping-mode"
                | "--sites-only"
        ) {
            bail!("unknown option: {key}");
        }
    }
    let required = |key: &str| {
        opts.get(key)
            .map(String::as_str)
            .with_context(|| format!("missing {key}"))
    };
    let pack_dir = PathBuf::from(required("--reference-pack")?);
    let genome = PathBuf::from(required("--genome")?);
    let gff = PathBuf::from(required("--gff")?);
    let out_dir = PathBuf::from(required("--out-dir")?);
    let (sequence, is_cds) = match (opts.get("--proteins"), opts.get("--cds")) {
        (Some(path), None) => (PathBuf::from(path), false),
        (None, Some(path)) => (PathBuf::from(path), true),
        _ => bail!("supply exactly one of --proteins or --cds"),
    };
    let include_fallback = bool_option(&opts, "--include-single-domain-fallback");
    let mapping_mode = opts
        .get("--mapping-mode")
        .cloned()
        .unwrap_or_else(|| "strict".to_string());
    if !matches!(mapping_mode.as_str(), "strict" | "balanced" | "exploratory") {
        bail!("--mapping-mode must be strict, balanced or exploratory");
    }
    let pack = load_pack(&pack_dir)?;
    let promoter_len = opts
        .get("--promoter-len")
        .map(|value| value.parse::<usize>().context("invalid --promoter-len"))
        .transpose()?
        .or(pack.manifest.defaults.promoter_len)
        .unwrap_or(1000);
    if promoter_len == 0 {
        bail!("--promoter-len must be > 0");
    }
    let pvalue = opts
        .get("--pvalue")
        .map(|value| value.parse::<f64>().context("invalid --pvalue"))
        .transpose()?
        .or(pack.manifest.defaults.pvalue)
        .unwrap_or(1e-4);
    if !(0.0..=1.0).contains(&pvalue) {
        bail!("--pvalue must be in [0, 1]");
    }
    let background = opts
        .get("--background")
        .map(|value| super::parse_background(value))
        .transpose()?
        .unwrap_or([0.25; 4]);
    // The promoter A/C/G/T composition is the correct null: a uniform i.i.d.
    // background mis-calibrates every p-value on species whose promoters are
    // A/T-rich (measured at r=0.24 against observed per-motif coverage on HF1,
    // versus r=0.94 for the composition-matched null). `uniform` stays available
    // for comparability with results produced before this default changed.
    let background_mode = match opts.get("--background-mode") {
        Some(mode) => mode.clone(),
        // An explicit --background supplies the whole null, so it keeps the
        // uniform model it has always described.
        None if opts.contains_key("--background") => "uniform".to_string(),
        None => "empirical".to_string(),
    };
    if !matches!(background_mode.as_str(), "uniform" | "empirical") {
        bail!("--background-mode must be uniform or empirical");
    }
    if background_mode == "empirical" && opts.contains_key("--background") {
        bail!("--background and --background-mode empirical are mutually exclusive");
    }
    let reuse_exact_pwm = bool_option(&opts, "--reuse-exact-pwm");
    let phase3_artifacts = bool_option(&opts, "--phase3-artifacts");
    let sites_only = bool_option(&opts, "--sites-only");
    if sites_only && !phase3_artifacts {
        bail!("--sites-only requires --phase3-artifacts");
    }
    let gene_biotype = opts
        .get("--gene-biotype")
        .cloned()
        .or_else(|| pack.manifest.defaults.gene_biotype.clone());
    let gene_biotype_from_pack =
        opts.get("--gene-biotype").is_none() && gene_biotype.is_some();
    let threads = opts
        .get("--threads")
        .map(|value| value.parse::<usize>().context("invalid --threads"))
        .transpose()?;
    if threads == Some(0) {
        bail!("--threads must be a positive integer");
    }
    let hmmsearch = opts
        .get("--hmmsearch")
        .cloned()
        .unwrap_or_else(|| pack.tool("hmmsearch"));
    let blastp = opts
        .get("--blastp")
        .cloned()
        .unwrap_or_else(|| pack.tool("blastp"));
    let hmmsearch = if Path::new(&hmmsearch).components().count() > 1 {
        resolve_tool(&pack.dir, &hmmsearch)
    } else {
        PathBuf::from(hmmsearch)
    };
    let blastp = if Path::new(&blastp).components().count() > 1 {
        resolve_tool(&pack.dir, &blastp)
    } else {
        PathBuf::from(blastp)
    };

    let outputs = output_paths(&out_dir);
    if outputs.values().any(|path| path.exists()) {
        bail!(
            "output directory already contains pipeline outputs; use a new --out-dir: {}",
            out_dir.display()
        );
    }
    fs::create_dir_all(&out_dir)
        .with_context(|| format!("create output directory {}", out_dir.display()))?;

    println!(
        "reference pack {} ({}) {}",
        pack.manifest.pack_version,
        pack.manifest.species,
        pack.manifest.assembly
    );

    report_pipeline_stage(&outputs["progress"], "identifying_tfs", "Identifying transcription factors")?;
    // 1. TF catalogue from protein/CDS + GFF; no user gene-map.
    let mut args = vec![
        "--gff".to_owned(),
        gff.to_string_lossy().into_owned(),
        "--hmm".to_owned(),
        pack.path("hmm").to_string_lossy().into_owned(),
        "--rules".to_owned(),
        pack.path("rules").to_string_lossy().into_owned(),
        "--hmmsearch".to_owned(),
        hmmsearch.to_string_lossy().into_owned(),
        "--out".to_owned(),
        outputs["tf_genes"].to_string_lossy().into_owned(),
    ];
    if is_cds {
        args.extend(["--cds".to_owned(), sequence.to_string_lossy().into_owned()]);
    } else {
        args.extend([
            "--proteins".to_owned(),
            sequence.to_string_lossy().into_owned(),
        ]);
    }
    if let Some(threads) = threads {
        args.extend(["--threads".to_owned(), threads.to_string()]);
    }
    plant_tfs::run(&args).context("TF identification failed")?;
    let genes = io::read_tf_genes_allow_empty(&outputs["tf_genes"])
        .with_context(|| format!("read {}", outputs["tf_genes"].display()))?;

    report_pipeline_stage(&outputs["progress"], "mapping_motifs", "Finding reference motif assignments")?;
    // 2. TF -> motif mapping from the same protein/CDS FASTA.
    let strict_map = out_dir.join("tf_map.strict.tsv");
    let mut args = vec![
        "--tf-genes".to_owned(),
        outputs["tf_genes"].to_string_lossy().into_owned(),
        "--reference-map".to_owned(),
        pack.path("reference_map").to_string_lossy().into_owned(),
        "--blastp".to_owned(),
        blastp.to_string_lossy().into_owned(),
        "--out".to_owned(),
        strict_map.to_string_lossy().into_owned(),
        "--diagnostics-out".to_owned(),
        outputs["mapping_diagnostics"].to_string_lossy().into_owned(),
        "--diagnostics-family-out".to_owned(),
        outputs["mapping_diagnostics_family"].to_string_lossy().into_owned(),
    ];
    if pack.paths.contains_key("reference_proteins") {
        args.extend(["--reference-proteins".to_owned(), pack.path("reference_proteins").display().to_string()]);
    } else {
        args.extend([
            "--reference-domains".to_owned(), pack.path("reference_domains").display().to_string(),
            "--dbd-models".to_owned(), pack.path("dbd_models").display().to_string(),
        ]);
    }
    if is_cds {
        args.extend(["--cds".to_owned(), sequence.to_string_lossy().into_owned()]);
    } else {
        args.extend([
            "--proteins".to_owned(),
            sequence.to_string_lossy().into_owned(),
        ]);
    }
    motif_transfer::run(&args).context("TF–motif mapping failed")?;
    if mapping_mode == "strict" {
        fs::rename(&strict_map, &outputs["tf_map"]).with_context(|| {
            format!(
                "move strict TF–motif mapping {} to {}",
                strict_map.display(),
                outputs["tf_map"].display()
            )
        })?;
    } else {
        report_pipeline_stage(&outputs["progress"], "mapping_dbd", "Evaluating domain-supported motif candidates")?;
        let calibrated_models = optional_pack_path(&pack, "calibrated_dbd_models")?;
        let calibrated_thresholds = optional_pack_path(&pack, "calibrated_thresholds")?;
        optional_pack_path(&pack, "specificity_groups")?;
        let calibrated_map = out_dir.join("tf_map.calibrated.tsv");
        run_dbd_mapping(
            &pack,
            &outputs["tf_genes"],
            &sequence,
            is_cds,
            &blastp,
            calibrated_models,
            calibrated_thresholds,
            &calibrated_map,
            false,
            true,
        )?;
        let mut extras = vec![calibrated_map];
        if mapping_mode == "exploratory" {
            let exploratory_thresholds = optional_pack_path(&pack, "exploratory_thresholds")?;
            let exploratory_map = out_dir.join("tf_map.exploratory.tsv");
            run_dbd_mapping(
                &pack,
                &outputs["tf_genes"],
                &sequence,
                is_cds,
                &blastp,
                pack.path("dbd_models"),
                exploratory_thresholds,
                &exploratory_map,
                true,
                false,
            )?;
            extras.push(exploratory_map);
        }
        merge_tf_maps(&strict_map, &extras, &outputs["tf_map"])?;
    }
    let tf_map = io::read_tf_map_allow_empty(&outputs["tf_map"]).with_context(|| {
        format!(
            "read {}; the pack produced no usable TF–motif mappings",
            outputs["tf_map"].display()
        )
    })?;
    let mapping_diagnostics_rows = read_tsv_rows(&outputs["mapping_diagnostics"])
        .with_context(|| format!("read {}", outputs["mapping_diagnostics"].display()))?;
    let mapping_family_rows = read_tsv_rows(&outputs["mapping_diagnostics_family"])
        .with_context(|| format!("read {}", outputs["mapping_diagnostics_family"].display()))?;
    let map_rows = tf_map.len();
    let map_tfs: BTreeSet<&str> = tf_map.iter().map(|tf| tf.tf_id.as_str()).collect();
    let full_rows = tf_map
        .iter()
        .filter(|tf| tf.domain_mode.eq_ignore_ascii_case("full_architecture"))
        .count();
    let fallback_rows = tf_map
        .iter()
        .filter(|tf| is_single_domain_fallback(tf))
        .count();
    let (default_map, default_stats) = filter_scan_mappings(tf_map.clone(), false);
    if !include_fallback && map_rows > 0 && default_stats.retained == 0 {
        bail!(
            "all {} TF–motif mapping row(s) are single_domain_fallback; rerun with --include-single-domain-fallback to scan them",
            map_rows
        );
    }
    let (actual_map, actual_stats) = if include_fallback {
        filter_scan_mappings(tf_map.clone(), true)
    } else {
        (default_map, default_stats)
    };

    // 3. Promoter scan with the pack motif library and the generated catalogue/map.
    let mut args = vec![
        "--genome".to_owned(),
        genome.to_string_lossy().into_owned(),
        "--gff".to_owned(),
        gff.to_string_lossy().into_owned(),
        "--motifs".to_owned(),
        pack.path("motifs").to_string_lossy().into_owned(),
        "--tf-genes".to_owned(),
        outputs["tf_genes"].to_string_lossy().into_owned(),
        "--tf-map".to_owned(),
        outputs["tf_map"].to_string_lossy().into_owned(),
        "--promoter-len".to_owned(),
        promoter_len.to_string(),
        "--pvalue".to_owned(),
        pvalue.to_string(),
        "--promoters-out".to_owned(),
        outputs["promoters"].to_string_lossy().into_owned(),
        "--progress-json".to_owned(),
        outputs["progress"].to_string_lossy().into_owned(),
        "--report-json".to_owned(),
        outputs["scan_report"].to_string_lossy().into_owned(),
        "--out".to_owned(),
        outputs["evidence"].to_string_lossy().into_owned(),
    ];
    if background_mode == "uniform" {
        args.extend([
            "--background".to_owned(),
            format!(
                "{},{},{},{}",
                background[0], background[1], background[2], background[3]
            ),
        ]);
    } else {
        args.extend([
            "--background-mode".to_owned(),
            background_mode.clone(),
        ]);
    }
    if reuse_exact_pwm {
        args.push("--reuse-exact-pwm".to_owned());
    }
    if let Some(cache) = opts.get("--scan-cache") {
        args.extend(["--scan-cache".to_owned(), cache.clone()]);
    }
    if phase3_artifacts {
        args.extend([
            "--sites-out".to_owned(),
            outputs["sites"].to_string_lossy().into_owned(),
            "--motif-testability-out".to_owned(),
            outputs["motif_testability"].to_string_lossy().into_owned(),
            "--promoter-summary-out".to_owned(),
            outputs["promoter_summary"].to_string_lossy().into_owned(),
            "--scan-manifest".to_owned(),
            outputs["scan_manifest"].to_string_lossy().into_owned(),
        ]);
    }
    if include_fallback {
        args.push("--include-single-domain-fallback".to_owned());
    }
    if sites_only {
        args.push("--sites-only".to_owned());
    }
    if let Some(downstream) = opts.get("--downstream") {
        args.extend(["--downstream".to_owned(), downstream.clone()]);
    }
    if let Some(biotype) = &gene_biotype {
        args.extend(["--gene-biotype".to_owned(), biotype.clone()]);
        if gene_biotype_from_pack {
            args.push("--gene-biotype-from-pack".to_owned());
        }
    }
    if let Some(threads) = threads {
        args.extend(["--threads".to_owned(), threads.to_string()]);
    }
    if let Some(annotation) = opts.get("--annotation") {
        args.extend(["--annotation".to_owned(), annotation.clone()]);
    }
    if tf_map.is_empty() {
        args.push("--allow-empty-tf-map".to_owned());
    }
    cmd_scan(&args).context("promoter scan failed")?;
    let scan_report: serde_json::Value = serde_json::from_reader(
        File::open(&outputs["scan_report"]).with_context(|| {
            format!("open scan report {}", outputs["scan_report"].display())
        })?,
    )
    .with_context(|| format!("parse scan report {}", outputs["scan_report"].display()))?;
    let scan_context = scan_summary_context(&scan_report);
    let effective_background = scan_context
        .get("background_values")
        .and_then(|value| parse_background_context(value))
        .unwrap_or(background);
    let used_phase3_scan = scan_context
        .get("phase3_scan")
        .is_some_and(|value| value == "true");
    let scan_manifest: serde_json::Value = if outputs["scan_manifest"].exists() {
        serde_json::from_reader(
            File::open(&outputs["scan_manifest"]).with_context(|| {
                format!("open scan manifest {}", outputs["scan_manifest"].display())
            })?,
        )
        .with_context(|| {
            format!(
                "parse scan manifest {}",
                outputs["scan_manifest"].display()
            )
        })?
    } else {
        serde_json::Value::Null
    };
    let gene_biotype_status =
        gene_biotype_status(&gene_biotype, gene_biotype_from_pack, &scan_report);

    // 4. SQLite evidence database in the hits/pairs schema used by prior.
    // Site-only jobs intentionally keep evidence empty and use the site table;
    // the schema-v3 background export is built from sites + policy mappings.
    if !sites_only {
        let args = vec![
            "--input".to_owned(),
            outputs["evidence"].to_string_lossy().into_owned(),
            "--out".to_owned(),
            outputs["database"].to_string_lossy().into_owned(),
            "--format".to_owned(),
            "sqlite".to_owned(),
        ];
        cmd_convert(&args).context("SQLite conversion failed")?;
    }

    let edges = io::read_edges(&outputs["evidence"])
        .with_context(|| format!("read {}", outputs["evidence"].display()))?;
    let motifs = pwm::read_jaspar(pack.path("motifs"))
        .with_context(|| format!("read {}", pack.path("motifs").display()))?;

    // 4b. Derived trace views: local TF -> motif provenance -> evidence.
    // Keep these small and derived; evidence.tsv remains the site-level table.
    let reference_map = io::read_tf_map(pack.path("reference_map"))
        .with_context(|| format!("read {}", pack.path("reference_map").display()))?;
    let scanned_pairs: BTreeSet<(String, String)> = actual_map
        .iter()
        .map(|mapping| (mapping.tf_id.clone(), mapping.motif_id.clone()))
        .collect();
    let motif_library = pack
        .path("motifs")
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let trace_summary = trace::write_trace_outputs(
        &out_dir,
        &trace::TraceInput {
            genes: &genes,
            mappings: &tf_map,
            scanned_pairs: &scanned_pairs,
            reference_map: &reference_map,
            motifs: &motifs,
            edges: &edges,
            include_fallback,
            pack_version: &pack.manifest.pack_version,
            motif_library: &motif_library,
        },
    )
    .context("write trace outputs")?;

    let pack_motif_ids: BTreeSet<String> = motifs.iter().map(|motif| motif.id.clone()).collect();
    let min_pvalue_by_motif: BTreeMap<String, f64> = motifs
        .iter()
        .map(|motif| {
            let min_p =
                scan::Scanner::new(motif, &effective_background, 10.0).min_pvalue();
            (motif.id.clone(), min_p)
        })
        .collect();
    let hit_motif_ids: BTreeSet<String> = edges
        .iter()
        .map(|edge| edge.motif_id.clone())
        .filter(|motif| !motif.is_empty())
        .collect();
    let scanned_motif_ids: BTreeSet<String> =
        actual_map.iter().map(|tf| tf.motif_id.clone()).collect();
    let pack_without_hits: Vec<_> = pack_motif_ids
        .difference(&hit_motif_ids)
        .cloned()
        .collect();
    let scanned_without_hits: Vec<_> = scanned_motif_ids
        .difference(&hit_motif_ids)
        .cloned()
        .collect();
    let excluded_map: Vec<&io::TfMotif> = if include_fallback {
        Vec::new()
    } else {
        tf_map
            .iter()
            .filter(|tf| is_single_domain_fallback(tf))
            .collect()
    };
    let excluded_motif_ids: BTreeSet<String> = excluded_map
        .iter()
        .map(|tf| tf.motif_id.clone())
        .filter(|motif| !scanned_motif_ids.contains(motif))
        .collect();
    let excluded_fallback_motifs: Vec<_> = excluded_motif_ids.into_iter().collect();
    let scanned_unreachable: Vec<_> = scanned_motif_ids
        .iter()
        .filter(|motif| {
            min_pvalue_by_motif
                .get(motif.as_str())
                .is_some_and(|min_p| *min_p > pvalue)
        })
        .cloned()
        .collect();
    let pack_unreachable: Vec<_> = pack_motif_ids
        .iter()
        .filter(|motif| {
            min_pvalue_by_motif
                .get(motif.as_str())
                .is_some_and(|min_p| *min_p > pvalue)
        })
        .cloned()
        .collect();
    let fallback_report: Vec<_> = tf_map
        .iter()
        .filter(|tf| is_single_domain_fallback(tf))
        .map(|tf| {
            serde_json::json!({
                "tf_id": tf.tf_id,
                "tf_symbol": tf.tf_symbol,
                "motif_id": tf.motif_id,
                "extraction_family": tf.extraction_family,
                "mapping_method": tf.mapping_method,
            })
        })
        .collect();
    let identified: Vec<_> = genes.values().filter(|gene| gene.is_identified()).collect();
    let tf_without_motif: Vec<_> = identified
        .iter()
        .filter(|gene| !map_tfs.contains(gene.gene_id.as_str()))
        .map(|gene| gene.gene_id.clone())
        .collect();
    let mut ambiguous_genes: Vec<String> = genes
        .values()
        .filter(|gene| gene.status == "ambiguous")
        .map(|gene| gene.gene_id.clone())
        .collect();
    ambiguous_genes.sort();
    let mut rejected_genes: Vec<String> = genes
        .values()
        .filter(|gene| gene.status == "rejected")
        .map(|gene| gene.gene_id.clone())
        .collect();
    rejected_genes.sort();
    let mut motifs_by_tf: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for mapping in &tf_map {
        motifs_by_tf
            .entry(mapping.tf_id.as_str())
            .or_default()
            .insert(mapping.motif_id.as_str());
    }
    let tf_with_multiple_motifs = motifs_by_tf.values().filter(|motifs| motifs.len() > 1).count();
    let pack_motif_id_set: BTreeSet<&str> = motifs.iter().map(|motif| motif.id.as_str()).collect();
    let mut scanned_tfs_by_motif: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for mapping in &actual_map {
        if !pack_motif_id_set.contains(mapping.motif_id.as_str()) {
            continue;
        }
        scanned_tfs_by_motif
            .entry(mapping.motif_id.as_str())
            .or_default()
            .insert(mapping.tf_id.as_str());
    }
    let scanned_motifs_shared = scanned_tfs_by_motif
        .values()
        .filter(|tfs| tfs.len() > 1)
        .count();
    let tied_examples: Vec<_> = tf_map
        .iter()
        .filter(|mapping| {
            mapping
                .mapping_method
                .to_ascii_lowercase()
                .contains("tied")
        })
        .take(20)
        .map(|mapping| {
            serde_json::json!({
                "tf_id": mapping.tf_id,
                "tf_symbol": mapping.tf_symbol,
                "motif_id": mapping.motif_id,
                "reference_tf_id": mapping.reference_tf_id,
                "mapping_method": mapping.mapping_method,
            })
        })
        .collect();
    let report = serde_json::json!({
        "reference_pack": {
            "version": pack.manifest.pack_version,
            "species": pack.manifest.species,
            "assembly": pack.manifest.assembly,
            "motif_scope": pack.manifest.motif_scope,
            "directory": pack.dir,
        },
        "inputs": {
            "genome": genome,
            "gff": gff,
            "sequences": sequence,
            "sequence_type": if is_cds { "cds" } else { "protein" },
        },
        "parameters": {
            "promoter_len": promoter_len,
            "pvalue": pvalue,
            "background": effective_background,
            "background_mode": background_mode,
            "exact_pwm_reuse": reuse_exact_pwm,
            "phase3_artifacts": phase3_artifacts,
            "include_single_domain_fallback": include_fallback,
            "gene_biotype": &gene_biotype,
            "gene_biotype_from_pack": gene_biotype_from_pack,
            "threads": threads,
        },
        "tools": {
            "hmmsearch": hmmsearch,
            "blastp": blastp,
        },
        "tf_catalogue": {
            "genes": genes.len(),
            "identified": identified.len(),
            "ambiguous": ambiguous_genes.len(),
            "rejected": rejected_genes.len(),
            "ambiguous_examples": ambiguous_genes.iter().take(100).collect::<Vec<_>>(),
            "rejected_examples": rejected_genes.iter().take(100).collect::<Vec<_>>(),
        },
        "tf_map": {
            "method": if pack.paths.contains_key("reference_proteins") { "same_family_reciprocal_protein" } else { "reciprocal_dbd" },
            "rows": map_rows,
            "unique_tfs": map_tfs.len(),
            "full_architecture_rows": full_rows,
            "single_domain_fallback_rows": fallback_rows,
            "tf_with_multiple_motifs": tf_with_multiple_motifs,
            "motifs_shared_by_multiple_tfs": scanned_motifs_shared,
            "motifs_shared_by_multiple_tfs_all_mappings": trace_summary.motifs_with_multiple_tfs,
            "tied_mapping_rows": trace_summary.tied_mapping_rows,
            "tied_examples": tied_examples,
            "tf_without_motif": {
                "count": tf_without_motif.len(),
                "examples": tf_without_motif.iter().take(100).collect::<Vec<_>>(),
            },
        },
        "mapping_diagnostics": {
            "tf_table": outputs["mapping_diagnostics"].display().to_string(),
            "family_table": outputs["mapping_diagnostics_family"].display().to_string(),
            "tf_rows": mapping_diagnostics_rows.len(),
            "family_rows": mapping_family_rows.len(),
            "families": mapping_family_rows,
        },
        "trace": {
            "tf_motif_rows": trace_summary.tf_motif_rows,
            "ambiguous_tf_motif_rows": trace_summary.ambiguous_tf_motif_rows,
            "motifs_with_multiple_tfs": trace_summary.motifs_with_multiple_tfs,
            "scanned_motifs_with_multiple_tfs": scanned_motifs_shared,
            "tied_mapping_rows": trace_summary.tied_mapping_rows,
            "site_tie_groups": trace_summary.site_tie_groups,
            "site_tie_rows": trace_summary.site_tie_rows,
            "outputs": {
                "tf_motif_trace": outputs["tf_motif_trace"].display().to_string(),
                "motif_ambiguity": outputs["motif_ambiguity"].display().to_string(),
                "target_candidate_ties": outputs["target_candidate_ties"].display().to_string(),
            },
        },
        "scan": {
            "mode": if pack.paths.contains_key("reference_proteins") { "representative_motif" } else if include_fallback { "include_single_domain_fallback" } else { "full_architecture_only" },
            "evidence_rows": edges.len(),
            "actual_retained_rows": actual_stats.retained,
            "actual_excluded_rows": actual_stats.excluded,
            "default_retained_rows": default_stats.retained,
            "default_excluded_rows": default_stats.excluded,
            "fallback_rows": fallback_report,
            "excluded_fallback_rows": if include_fallback { serde_json::json!([]) } else { serde_json::json!(fallback_report) },
            "gene_biotype": gene_biotype_status,
        },
        "phase3": {
            "enabled": used_phase3_scan,
            "background_mode": scan_context.get("background_mode"),
            "unique_site_hits": scan_context.get("unique_site_hits").and_then(|value| value.parse::<u64>().ok()),
            "motif_id_hits": scan_context.get("motif_id_hits").and_then(|value| value.parse::<u64>().ok()),
            "evidence_rows": scan_context.get("evidence_rows").and_then(|value| value.parse::<u64>().ok()),
            "pwm_groups": scan_context.get("pwm_groups").and_then(|value| value.parse::<u64>().ok()),
            "unique_pwms": scan_context.get("unique_pwms").and_then(|value| value.parse::<u64>().ok()),
            "duplicate_pwm_groups": scan_context.get("duplicate_pwm_groups").and_then(|value| value.parse::<u64>().ok()),
            "cache_hits": scan_context.get("cache_hits").and_then(|value| value.parse::<u64>().ok()),
            "cache_misses": scan_context.get("cache_misses").and_then(|value| value.parse::<u64>().ok()),
            "testable_motifs": scan_context.get("testable_motifs").and_then(|value| value.parse::<u64>().ok()),
            "unattainable_motifs": scan_context.get("unattainable_motifs").and_then(|value| value.parse::<u64>().ok()),
            "artifacts": {
                "sites": outputs["sites"].display().to_string(),
                "motif_testability": outputs["motif_testability"].display().to_string(),
                "promoter_summary": outputs["promoter_summary"].display().to_string(),
                "scan_manifest": outputs["scan_manifest"].display().to_string(),
            },
            "manifest": scan_manifest,
        },
        "motifs": {
            "pack_motifs": pack_motif_ids.len(),
            "with_hits": pack_motif_ids.len() - pack_without_hits.len(),
            "pack_without_hits": {
                "count": pack_without_hits.len(),
                "examples": &pack_without_hits,
            },
            "scanned_without_hits": {
                "count": scanned_without_hits.len(),
                "examples": &scanned_without_hits,
            },
            "excluded_fallback_motifs": {
                "count": excluded_fallback_motifs.len(),
                "examples": &excluded_fallback_motifs,
            },
            "threshold": {
                "requested_pvalue": pvalue,
                "scanned_unreachable": {
                    "count": scanned_unreachable.len(),
                    "examples": &scanned_unreachable,
                },
                "pack_unreachable": {
                    "count": pack_unreachable.len(),
                    "examples": &pack_unreachable,
                },
                "min_pvalue_by_motif": &min_pvalue_by_motif,
            },
        },
        "outputs": outputs,
    });
    fs::write(
        &outputs["pipeline_report"],
        serde_json::to_vec_pretty(&report)?,
    )
    .with_context(|| format!("write {}", outputs["pipeline_report"].display()))?;

    println!(
        "TF catalogue: {} identified / {} total; {} ambiguous, {} rejected",
        identified.len(),
        genes.len(),
        ambiguous_genes.len(),
        rejected_genes.len()
    );
    if !ambiguous_genes.is_empty() {
        println!(
            "ambiguous TF genes excluded from mapping: {}{}",
            ambiguous_genes.iter().take(12).cloned().collect::<Vec<_>>().join(", "),
            if ambiguous_genes.len() > 12 { ", ..." } else { "" }
        );
    }
    println!(
        "TF–motif: {} rows, {} TFs, {} full, {} fallback; scan mode {}",
        map_rows,
        map_tfs.len(),
        full_rows,
        fallback_rows,
        if include_fallback { "include fallback" } else { "full only" }
    );
    println!(
        "evidence: {} rows; default would retain {} and exclude {}",
        edges.len(),
        default_stats.retained,
        default_stats.excluded
    );
    println!(
        "trace: {} TF–motif rows; {} ambiguous TF–motif row(s); {} scanned motif(s) shared by >1 TF; {} tied mapping row(s); {} site tie group(s) / {} candidate row(s)",
        trace_summary.tf_motif_rows,
        trace_summary.ambiguous_tf_motif_rows,
        scanned_motifs_shared,
        trace_summary.tied_mapping_rows,
        trace_summary.site_tie_groups,
        trace_summary.site_tie_rows
    );
    println!(
        "motifs: {} in pack, {} with hits, {} pack-wide without hits; {} scanned motifs without hits; {} fallback motifs excluded from scan",
        pack_motif_ids.len(),
        pack_motif_ids.len() - pack_without_hits.len(),
        pack_without_hits.len(),
        scanned_without_hits.len(),
        excluded_fallback_motifs.len()
    );
    if !scanned_without_hits.is_empty() {
        println!(
            "scanned motifs without hits: {}",
            scanned_without_hits
                .iter()
                .take(12)
                .cloned()
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    if !excluded_fallback_motifs.is_empty() {
        println!(
            "fallback motifs excluded from scan: {}",
            excluded_fallback_motifs
                .iter()
                .take(12)
                .cloned()
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    if !scanned_unreachable.is_empty() {
        println!(
            "WARNING: requested p={} is below the minimum achievable p for {} scanned motif(s): {}; scan continues with the requested threshold and these motifs cannot hit",
            pvalue,
            scanned_unreachable.len(),
            scanned_unreachable
                .iter()
                .take(50)
                .cloned()
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    println!("outputs: {}", out_dir.display());
    Ok(())
}

fn resolve_tool(pack_dir: &Path, value: &str) -> PathBuf {
    let path = Path::new(value);
    if path.is_absolute() {
        path.to_owned()
    } else {
        pack_dir.join(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sha2::{Digest, Sha256};

    fn write_sha(path: &Path, value: &str) -> String {
        let mut file = File::create(path).unwrap();
        file.write_all(value.as_bytes()).unwrap();
        format!("{:x}", Sha256::digest(value.as_bytes()))
    }

    fn fixture_content(key: &'static str) -> &'static str {
        if key == "motifs" {
            ">MA0001.1 TEST\nA [ 10 0 ]\nC [ 0 10 ]\nG [ 0 0 ]\nT [ 0 0 ]\n"
        } else {
            key
        }
    }

    fn write_manifest(dir: &Path, files: BTreeMap<String, String>, checksums: BTreeMap<String, String>) {
        let manifest = serde_json::json!({
            "pack_version": "test",
            "species": "test",
            "assembly": "test",
            "files": files,
            "checksums": checksums,
            "tools": {},
        });
        std::fs::write(dir.join("manifest.json"), serde_json::to_vec(&manifest).unwrap()).unwrap();
    }

    fn temp_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-pack-{name}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn pack_requires_all_files_and_checksums() {
        let dir = temp_dir("missing");
        let mut files = BTreeMap::new();
        let mut checksums = BTreeMap::new();
        for key in REQUIRED_FILES {
            let name = format!("{key}.dat");
            let hash = write_sha(&dir.join(&name), fixture_content(key));
            files.insert(key.to_owned(), name.clone());
            checksums.insert(name, hash);
        }
        write_manifest(&dir, files.clone(), checksums.clone());
        assert!(load_pack(&dir).is_ok());

        let mut bad = files.clone();
        bad.remove("hmm");
        write_manifest(&dir, bad, checksums);
        assert!(load_pack(&dir).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn pack_rejects_checksum_mismatch_and_path_escape() {
        let dir = temp_dir("bad");
        let mut files = BTreeMap::new();
        let mut checksums = BTreeMap::new();
        for key in REQUIRED_FILES {
            let name = format!("{key}.dat");
            let hash = write_sha(&dir.join(&name), fixture_content(key));
            files.insert(key.to_owned(), name.clone());
            checksums.insert(name, hash);
        }
        let mut bad = checksums.clone();
        bad.insert("hmm.dat".to_owned(), "0".repeat(64));
        write_manifest(&dir, files.clone(), bad);
        assert!(load_pack(&dir).is_err());

        let mut escape = files;
        escape.insert("hmm".to_owned(), "../outside.dat".to_owned());
        write_manifest(&dir, escape, checksums);
        assert!(load_pack(&dir).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn pack_rejects_duplicate_motif_ids() {
        let dir = temp_dir("duplicate-motifs");
        let mut files = BTreeMap::new();
        let mut checksums = BTreeMap::new();
        for key in REQUIRED_FILES {
            let name = format!("{key}.dat");
            let hash = write_sha(&dir.join(&name), fixture_content(key));
            files.insert(key.to_owned(), name.clone());
            checksums.insert(name, hash);
        }
        let motifs = files["motifs"].clone();
        let duplicate = ">DUP1 FIRST\nA [ 10 ]\nC [ 0 ]\nG [ 0 ]\nT [ 0 ]\n\
                         >DUP1 SECOND\nA [ 0 ]\nC [ 10 ]\nG [ 0 ]\nT [ 0 ]\n";
        let hash = write_sha(&dir.join(&motifs), duplicate);
        checksums.insert(motifs, hash);
        write_manifest(&dir, files, checksums);
        let err = format!("{:#}", load_pack(&dir).unwrap_err());
        assert!(err.contains("duplicate motif ID DUP1"), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn report_records_requested_and_applied_gene_biotype() {
        let applied = serde_json::json!({
            "diagnostics": [{"code": "GENE_BIOTYPE_FILTER"}]
        });
        let status = gene_biotype_status(
            &Some("protein_coding".to_owned()),
            true,
            &applied,
        );
        assert_eq!(status["status"], "applied");
        assert_eq!(status["source"], "pack_default");
        assert_eq!(status["applied"], true);

        let skipped = serde_json::json!({
            "diagnostics": [{"code": "GENE_BIOTYPE_FILTER_SKIPPED"}]
        });
        let status = gene_biotype_status(
            &Some("protein_coding".to_owned()),
            true,
            &skipped,
        );
        assert_eq!(status["status"], "skipped");
        assert_eq!(status["applied"], false);

        let none = serde_json::json!({"diagnostics": []});
        let status = gene_biotype_status(&None, false, &none);
        assert_eq!(status["status"], "not_requested");
        assert_eq!(status["source"], "none");
    }

    #[test]
    fn pipeline_boolean_flags_do_not_consume_following_options() {
        let args = vec![
            "--phase3-artifacts".to_owned(),
            "--sites-only".to_owned(),
            "--reuse-exact-pwm".to_owned(),
            "--promoter-len".to_owned(),
            "1000".to_owned(),
        ];
        let opts = parse_options_with_flags(&args, PIPELINE_FLAGS).unwrap();
        assert_eq!(opts.get("--sites-only").map(String::as_str), Some("true"));
        assert_eq!(
            opts.get("--reuse-exact-pwm").map(String::as_str),
            Some("true")
        );
        assert_eq!(opts.get("--promoter-len").map(String::as_str), Some("1000"));
        assert!(bool_option(&opts, "--sites-only"));
        assert!(bool_option(&opts, "--reuse-exact-pwm"));
    }

    #[test]
    fn merged_strict_map_gets_explicit_mapping_status() {
        let dir = temp_dir("merge-status");
        let strict = dir.join("strict.tsv");
        let calibrated = dir.join("calibrated.tsv");
        let out = dir.join("merged.tsv");
        std::fs::write(
            &strict,
            "tf_id\ttf_symbol\tfamily\tmotif_id\treference_tf_id\tmapping_method\tdomain_mode\tdomain_count\n\
             TF1\tTF1\tNAC\tMA0001.1\tD1\treciprocal_protein_similarity\twhole_protein\t0\n",
        )
        .unwrap();
        std::fs::write(
            &calibrated,
            "tf_id\ttf_symbol\tfamily\tmotif_id\treference_tf_id\tmapping_method\tdomain_mode\tdomain_count\tmapping_status\n\
             TF2\tTF2\tNAC\tMA0002.1\tD2\tdbd_transfer_calibrated\tfull_architecture\t1\tdbd_calibrated\n",
        )
        .unwrap();
        let rows = merge_tf_maps(&strict, &[calibrated], &out).unwrap();
        assert_eq!(rows, 2);
        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.contains("\tstrict_rbh"), "{text}");
        assert!(text.contains("dbd_calibrated"), "{text}");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
