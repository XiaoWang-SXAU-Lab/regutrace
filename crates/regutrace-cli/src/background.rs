//! Build the frozen schema-v3 regulatory background SQLite from a completed
//! site-only job. Policy rows are supplied by the dashboard so the UI and the
//! exported database use exactly the same selection implementation.
use super::*;
use std::io::Read;
use rusqlite::{params, params_from_iter, Connection};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet, HashMap};

const TF_MOTIF_COLUMNS: [&str; 29] = [
    "tf_id", "tf_symbol", "family", "motif_id", "reference_tf_id", "mapping_method",
    "domain_mode", "extraction_family", "extraction_models", "domain_count",
    "mapping_status", "alternative_donor_count", "mapping_confidence", "specificity_groups",
    "dbd_identity", "query_dbd_coverage", "donor_dbd_coverage", "reason_code",
    "query_full_dbd_start", "query_full_dbd_end", "donor_full_dbd_start", "donor_full_dbd_end",
    "query_dbd_start", "query_dbd_end", "donor_dbd_start", "donor_dbd_end", "dbd_bitscore",
    "dbd_evalue", "alignment_length",
];

const HITS_SQL: &str = "
CREATE TABLE IF NOT EXISTS hits (
 tf_id TEXT NOT NULL, tf_symbol TEXT NOT NULL DEFAULT '',
 target_id TEXT NOT NULL, target_symbol TEXT NOT NULL DEFAULT '',
 evidence_type TEXT NOT NULL DEFAULT '', evidence_level TEXT NOT NULL DEFAULT '',
 method TEXT NOT NULL DEFAULT '', source TEXT NOT NULL DEFAULT '',
 pmid TEXT NOT NULL DEFAULT '', condition TEXT NOT NULL DEFAULT '',
 promoter_start INTEGER, promoter_end INTEGER,
 motif_id TEXT NOT NULL DEFAULT '', hit_start INTEGER, hit_end INTEGER,
 strand TEXT NOT NULL DEFAULT '', score REAL, p_value REAL,
 notes TEXT NOT NULL DEFAULT ''
);
CREATE TABLE IF NOT EXISTS pairs (
 tf_id TEXT NOT NULL, tf_symbol TEXT NOT NULL, target_id TEXT NOT NULL,
 target_symbol TEXT NOT NULL, n_sites INTEGER NOT NULL, n_motifs INTEGER NOT NULL,
 best_p REAL, best_score REAL, motifs TEXT, evidence_levels TEXT
);";

fn sha256(path: &Path) -> Result<String> {
    let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn read_tsv(path: &Path) -> Result<(Vec<String>, Vec<BTreeMap<String, String>>)> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("open {}", path.display()))?;
    let headers: Vec<String> = reader
        .headers()?
        .iter()
        .map(|value| value.trim().trim_start_matches('\u{feff}').to_string())
        .collect();
    let mut rows = Vec::new();
    for record in reader.records() {
        let record = record?;
        let mut row = BTreeMap::new();
        for (i, header) in headers.iter().enumerate() {
            row.insert(header.clone(), record.get(i).unwrap_or("").to_string());
        }
        rows.push(row);
    }
    Ok((headers, rows))
}

fn quote(name: &str) -> String {
    format!("\"{}\"", name.replace('"', "\"\""))
}

fn create_text_table(conn: &Connection, table: &str, columns: &[String]) -> Result<()> {
    conn.execute(&format!("DROP TABLE IF EXISTS {table}"), [])?;
    let cols = columns
        .iter()
        .map(|name| format!("{} TEXT", quote(name)))
        .collect::<Vec<_>>()
        .join(", ");
    conn.execute(&format!("CREATE TABLE {table} ({cols})"), [])?;
    Ok(())
}

fn insert_rows(
    conn: &Connection,
    table: &str,
    columns: &[String],
    rows: &[BTreeMap<String, String>],
) -> Result<()> {
    let placeholders = std::iter::repeat("?")
        .take(columns.len())
        .collect::<Vec<_>>()
        .join(",");
    let sql = format!(
        "INSERT INTO {table} ({}) VALUES ({placeholders})",
        columns.iter().map(|name| quote(name)).collect::<Vec<_>>().join(",")
    );
    let mut stmt = conn.prepare(&sql)?;
    for row in rows {
        let values = columns
            .iter()
            .map(|name| Value::Text(row.get(name).cloned().unwrap_or_default()))
            .collect::<Vec<_>>();
        stmt.execute(params_from_iter(values))?;
    }
    Ok(())
}

fn copy_table(conn: &Connection, path: &Path, table: &str, wanted: &[&str]) -> Result<usize> {
    let (headers, rows) = read_tsv(path)?;
    let missing: Vec<&str> = wanted
        .iter()
        .copied()
        .filter(|name| !headers.iter().any(|header| header == name))
        .collect();
    if !missing.is_empty() {
        bail!("{}: missing columns {:?}", path.display(), missing);
    }
    let columns: Vec<String> = wanted.iter().map(|name| (*name).to_string()).collect();
    create_text_table(conn, table, &columns)?;
    insert_rows(conn, table, &columns, &rows)?;
    Ok(rows.len())
}

fn policy_rows(conn: &Connection, path: &Path) -> Result<usize> {
    let (headers, rows) = read_tsv(path)?;
    let wanted = [
        "policy", "policy_version", "tf_id", "motif_id", "selection_role",
        "selection_source", "selection_group", "specificity_group", "reference_tf_id",
        "selection_order",
    ];
    let missing: Vec<&str> = wanted
        .iter()
        .copied()
        .filter(|name| !headers.iter().any(|header| header == name))
        .collect();
    if !missing.is_empty() {
        bail!("{}: missing policy columns {:?}", path.display(), missing);
    }
    conn.execute("DROP TABLE IF EXISTS tf_motif_policy", [])?;
    conn.execute(
        "CREATE TABLE tf_motif_policy (
            policy TEXT NOT NULL, policy_version TEXT NOT NULL, tf_id TEXT NOT NULL,
            motif_id TEXT NOT NULL, selection_role TEXT NOT NULL,
            selection_source TEXT NOT NULL, selection_group TEXT, specificity_group TEXT,
            reference_tf_id TEXT, selection_order INTEGER NOT NULL,
            PRIMARY KEY (policy, tf_id, motif_id)
        )",
        [],
    )?;
    let sql = "INSERT INTO tf_motif_policy VALUES (?,?,?,?,?,?,?,?,?,?)";
    let mut stmt = conn.prepare(sql)?;
    for row in &rows {
        stmt.execute(params![
            row.get("policy").cloned().unwrap_or_default(),
            row.get("policy_version").cloned().unwrap_or_default(),
            row.get("tf_id").cloned().unwrap_or_default(),
            row.get("motif_id").cloned().unwrap_or_default(),
            row.get("selection_role").cloned().unwrap_or_default(),
            row.get("selection_source").cloned().unwrap_or_default(),
            row.get("selection_group").cloned().unwrap_or_default(),
            row.get("specificity_group").cloned().unwrap_or_default(),
            row.get("reference_tf_id").cloned().unwrap_or_default(),
            row.get("selection_order").and_then(|v| v.parse::<i64>().ok()).unwrap_or(0),
        ])?;
    }
    conn.execute("CREATE INDEX IF NOT EXISTS tmpolicy_tf ON tf_motif_policy(policy, tf_id)", [])?;
    conn.execute("CREATE INDEX IF NOT EXISTS tmpolicy_motif ON tf_motif_policy(policy, motif_id)", [])?;
    Ok(rows.len())
}

fn policy_config(conn: &Connection, path: &Path) -> Result<usize> {
    let (headers, rows) = read_tsv(path)?;
    let wanted = [
        "policy", "policy_version", "parameters_json", "source_file", "source_sha256",
        "coordinate_source", "coordinate_sha256", "row_count",
    ];
    let missing: Vec<&str> = wanted
        .iter()
        .copied()
        .filter(|name| !headers.iter().any(|header| header == name))
        .collect();
    if !missing.is_empty() {
        bail!("{}: missing policy config columns {:?}", path.display(), missing);
    }
    conn.execute("DROP TABLE IF EXISTS tf_motif_policy_config", [])?;
    conn.execute(
        "CREATE TABLE tf_motif_policy_config (
            policy TEXT PRIMARY KEY, policy_version TEXT NOT NULL,
            parameters_json TEXT NOT NULL, source_file TEXT NOT NULL,
            source_sha256 TEXT NOT NULL, coordinate_source TEXT NOT NULL,
            coordinate_sha256 TEXT NOT NULL, row_count INTEGER NOT NULL
        )",
        [],
    )?;
    let mut stmt = conn.prepare(
        "INSERT INTO tf_motif_policy_config VALUES (?,?,?,?,?,?,?,?)",
    )?;
    for row in &rows {
        stmt.execute(params![
            row.get("policy").cloned().unwrap_or_default(),
            row.get("policy_version").cloned().unwrap_or_default(),
            row.get("parameters_json").cloned().unwrap_or_default(),
            row.get("source_file").cloned().unwrap_or_default(),
            row.get("source_sha256").cloned().unwrap_or_default(),
            row.get("coordinate_source").cloned().unwrap_or_default(),
            row.get("coordinate_sha256").cloned().unwrap_or_default(),
            row.get("row_count").and_then(|v| v.parse::<i64>().ok()).unwrap_or(0),
        ])?;
    }
    Ok(rows.len())
}

fn registry(conn: &Connection, motifs: &Path, metadata: Option<&Path>) -> Result<usize> {
    let pwms = pwm::read_jaspar(motifs)?;
    let mut meta = HashMap::new();
    if let Some(path) = metadata {
        if path.exists() {
            for row in read_tsv(path)?.1 {
                meta.insert(
                    row.get("motif_id").cloned().unwrap_or_default(),
                    (
                        row.get("pwm_hash").cloned().unwrap_or_default(),
                        row.get("family").cloned().unwrap_or_default(),
                        row.get("group_080").cloned().unwrap_or_default(),
                    ),
                );
            }
        }
    }
    conn.execute(
        "CREATE TABLE motif_registry (motif_id TEXT PRIMARY KEY, motif_name TEXT, family TEXT, pwm_hash TEXT, matrix_length INTEGER, specificity_group_080 TEXT, source TEXT)",
        [],
    )?;
    let mut stmt = conn.prepare("INSERT OR REPLACE INTO motif_registry VALUES (?,?,?,?,?,?,?)")?;
    for pwm in &pwms {
        let (hash, family, group) = meta.get(&pwm.id).cloned().unwrap_or_default();
        stmt.execute(params![
            pwm.id, pwm.name, family, hash, pwm.matrix.len() as i64, group, "JASPAR2024"
        ])?;
    }
    Ok(pwms.len())
}

fn scan_config(conn: &Connection, path: &Path) -> Result<usize> {
    let manifest: serde_json::Value =
        serde_json::from_reader(File::open(path).with_context(|| format!("open {}", path.display()))?)
            .with_context(|| format!("parse {}", path.display()))?;
    conn.execute("CREATE TABLE scan_config (key TEXT PRIMARY KEY, value TEXT NOT NULL)", [])?;
    let get = |pointer: &str| manifest.pointer(pointer).cloned().unwrap_or(serde_json::Value::Null);
    let flat = [
        ("scan_manifest_path", serde_json::Value::String(path.file_name().map(|v| v.to_string_lossy().into_owned()).unwrap_or_default())),
        ("background_mode", get("/background/mode")),
        ("background_values_acgt", get("/background/values_acgt")),
        // Record the library identity, not a machine-local path: this row is
        // copied into databases that get shared and imported elsewhere.
        ("motif_library", serde_json::Value::String(
            Path::new(manifest.pointer("/inputs/motifs").and_then(|v| v.as_str()).unwrap_or(""))
                .file_name()
                .map(|v| v.to_string_lossy().into_owned())
                .unwrap_or_default(),
        )),
        ("exact_pwm_reuse", get("/motifs/exact_pwm_reuse")),
        ("pwm_groups", get("/motifs/pwm_groups")),
        ("mapped_motif_ids", get("/motifs/mapped_motif_ids")),
        ("promoter_count", get("/promoters/count")),
        ("scan_seconds", get("/runtime/seconds")),
        ("cache_hits", get("/scan_jobs/cache_hits")),
        ("cache_misses", get("/scan_jobs/cache_misses")),
        ("strands", get("/scan_algorithm/strands")),
    ];
    let flat_len = flat.len();
    let mut stmt = conn.prepare("INSERT OR REPLACE INTO scan_config VALUES (?,?)")?;
    for (key, value) in flat {
        let text = match value {
            serde_json::Value::Null => String::new(),
            serde_json::Value::String(value) => value,
            value => serde_json::to_string(&value)?,
        };
        stmt.execute(params![key, text])?;
    }
    stmt.execute(params!["sites_only", "true"])?;
    Ok(flat_len + 1)
}

/// How one motif's sites are counted per target.
///
/// `independent` merges overlapping hit windows, so a degenerate match that slides one base at a
/// time is one site rather than thirty. `overlapping` is the pre-2026-11 `COUNT(*)` and exists
/// only so a frozen artifact can be rebuilt byte-for-byte.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SiteCountMode {
    Independent,
    Overlapping,
}

impl SiteCountMode {
    fn parse(raw: &str) -> Result<Self> {
        match raw {
            "independent" => Ok(Self::Independent),
            "overlapping" => Ok(Self::Overlapping),
            other => bail!("--site-count must be independent or overlapping, got {other}"),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Independent => "independent",
            Self::Overlapping => "overlapping",
        }
    }
}

/// Greedy interval scheduling: the largest set of mutually non-overlapping windows. Sorting is by
/// end coordinate (then start, for determinism), which is what `n_sites` means after this change —
/// "how many binding events", not "how many windows passed the threshold".
fn independent_site_count(intervals: &mut Vec<(i64, i64)>) -> usize {
    intervals.sort_unstable_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));
    let mut count = 0;
    let mut last_end = i64::MIN;
    for &(start, end) in intervals.iter() {
        if start > last_end {
            count += 1;
            last_end = end;
        }
    }
    count
}

fn aggregate_motif_target(
    conn: &Connection,
    sites: &Path,
    mode: SiteCountMode,
) -> Result<usize> {
    conn.execute(
        "CREATE TEMP TABLE raw_site (
            motif_id TEXT, target_id TEXT, start INTEGER, end INTEGER,
            p_value REAL, score REAL)",
        [],
    )?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(sites)
        .with_context(|| format!("open {}", sites.display()))?;
    let headers = reader.headers()?.clone();
    let index = |names: &[&str]| {
        headers
            .iter()
            .position(|header| names.iter().any(|name| header == *name))
    };
    let motif_col = index(&["motif_ids", "motif_id"]).context("sites table has no motif_ids")?;
    let target_col = index(&["gene_id", "target_id"]).context("sites table has no gene_id")?;
    let p_col = index(&["p_value"]).context("sites table has no p_value")?;
    let score_col = index(&["score"]);
    let start_col = index(&["hit_start", "start"]);
    let end_col = index(&["hit_end", "end"]);
    conn.execute("BEGIN", [])?;
    {
        let mut stmt = conn.prepare("INSERT INTO raw_site VALUES (?,?,?,?,?,?)")?;
        for record in reader.records() {
            let record = record?;
            let target = record.get(target_col).unwrap_or("").trim();
            if target.is_empty() {
                continue;
            }
            let p = record
                .get(p_col)
                .and_then(|value| value.trim().parse::<f64>().ok())
                .unwrap_or(1.0);
            let score = score_col
                .and_then(|col| record.get(col))
                .and_then(|value| value.trim().parse::<f64>().ok())
                .unwrap_or(0.0);
            let start = start_col
                .and_then(|col| record.get(col))
                .and_then(|value| value.trim().parse::<i64>().ok());
            let end = end_col
                .and_then(|col| record.get(col))
                .and_then(|value| value.trim().parse::<i64>().ok());
            for motif in record.get(motif_col).unwrap_or("").split(';') {
                if !motif.is_empty() {
                    stmt.execute(params![motif, target, start, end, p, score])?;
                }
            }
        }
    }
    conn.execute("COMMIT", [])?;

    match mode {
        SiteCountMode::Overlapping => {
            conn.execute(
                "CREATE TABLE motif_target AS
                 SELECT motif_id, target_id, COUNT(*) n_sites, MIN(p_value) best_p,
                        MAX(score) best_score
                 FROM raw_site GROUP BY motif_id, target_id",
                [],
            )?;
        }
        SiteCountMode::Independent => merge_independent_sites(conn)?,
    }

    conn.execute("DROP TABLE raw_site", [])?;
    conn.execute(
        "CREATE INDEX IF NOT EXISTS motif_target_motif ON motif_target(motif_id)",
        [],
    )?;
    conn.execute(
        "CREATE INDEX IF NOT EXISTS motif_target_target ON motif_target(target_id)",
        [],
    )?;
    let n: i64 = conn.query_row("SELECT COUNT(*) FROM motif_target", [], |row| row.get(0))?;
    Ok(n as usize)
}

/// Builds `motif_target` counting independent sites. `best_p`/`best_score` still come from the raw
/// windows: merging changes how many sites there are, not which one scored best.
fn merge_independent_sites(conn: &Connection) -> Result<()> {
    conn.execute(
        "CREATE TEMP TABLE merged_motif_target (
            motif_id TEXT, target_id TEXT, n_sites INTEGER, best_p REAL, best_score REAL)",
        [],
    )?;
    conn.execute(
        "CREATE INDEX raw_site_order ON raw_site(motif_id, target_id, start)",
        [],
    )?;
    let mut select = conn.prepare(
        "SELECT motif_id, target_id, start, end, p_value, score
         FROM raw_site ORDER BY motif_id, target_id, start",
    )?;
    let mut insert = conn.prepare("INSERT INTO merged_motif_target VALUES (?,?,?,?,?)")?;
    let mut rows = select.query([])?;

    let mut current: Option<(String, String)> = None;
    let mut intervals: Vec<(i64, i64)> = Vec::new();
    // Rows without usable coordinates cannot be placed on the promoter, so each stands alone.
    let mut unplaced = 0_usize;
    let mut best_p = f64::INFINITY;
    let mut best_score = f64::NEG_INFINITY;

    while let Some(row) = rows.next()? {
        let motif: String = row.get(0)?;
        let target: String = row.get(1)?;
        let start: Option<i64> = row.get(2)?;
        let end: Option<i64> = row.get(3)?;
        let p: f64 = row.get(4)?;
        let score: f64 = row.get(5)?;
        let changed = match &current {
            Some((m, t)) => m != &motif || t != &target,
            None => true,
        };
        if changed {
            if let Some((m, t)) = current.take() {
                insert.execute(params![
                    m,
                    t,
                    independent_site_count(&mut intervals) as i64 + unplaced as i64,
                    best_p,
                    best_score
                ])?;
            }
            intervals.clear();
            unplaced = 0;
            best_p = f64::INFINITY;
            best_score = f64::NEG_INFINITY;
            current = Some((motif, target));
        }
        best_p = best_p.min(p);
        best_score = best_score.max(score);
        match (start, end) {
            (Some(start), Some(end)) => intervals.push((start, end)),
            _ => unplaced += 1,
        }
    }
    drop(rows);
    if let Some((m, t)) = current.take() {
        insert.execute(params![
            m,
            t,
            independent_site_count(&mut intervals) as i64 + unplaced as i64,
            best_p,
            best_score
        ])?;
    }
    drop(insert);
    drop(select);

    conn.execute(
        "CREATE TABLE motif_target AS
         SELECT motif_id, target_id, n_sites, best_p, best_score FROM merged_motif_target",
        [],
    )?;
    conn.execute("DROP TABLE merged_motif_target", [])?;
    Ok(())
}

fn manifest_requested_pvalue(path: &Path) -> Result<Option<f64>> {
    let manifest: serde_json::Value =
        serde_json::from_reader(File::open(path).with_context(|| format!("open {}", path.display()))?)
            .with_context(|| format!("parse {}", path.display()))?;
    Ok(manifest
        .pointer("/thresholds/requested_pvalue")
        .and_then(|value| value.as_f64()))
}

fn approximately_equal(a: f64, b: f64) -> bool {
    let scale = a.abs().max(b.abs()).max(f64::MIN_POSITIVE);
    (a - b).abs() / scale <= 1e-9
}

fn parse_usize_option(opts: &HashMap<String, String>, key: &str, default: usize) -> Result<usize> {
    match opts.get(key) {
        Some(raw) => raw
            .parse()
            .with_context(|| format!("invalid {key} {raw}")),
        None => Ok(default),
    }
}

pub(super) fn run(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let out = PathBuf::from(
        opts.get("--out")
            .map(String::as_str)
            .context("missing --out")?,
    );
    let tmp = temp_output_path(&out);
    let result = run_inner(&opts, &out, &tmp);
    if result.is_err() {
        // A failed build must not leave a multi-gigabyte partial database (or
        // its WAL sidecars) in the job directory; every other command that
        // publishes through temp_output_path cleans up the same way.
        cleanup_output_temp(&tmp);
    }
    result
}

fn run_inner(opts: &HashMap<String, String>, out: &Path, tmp: &Path) -> Result<()> {
    let get = |key: &str| -> Result<&str> {
        opts.get(key).map(String::as_str).with_context(|| format!("missing {key}"))
    };
    let tf_genes = PathBuf::from(get("--tf-genes")?);
    let tf_map = PathBuf::from(get("--tf-map")?);
    let sites = PathBuf::from(get("--sites")?);
    let testability = PathBuf::from(get("--testability")?);
    let promoter_summary = PathBuf::from(get("--promoter-summary")?);
    let scan_manifest = PathBuf::from(get("--scan-manifest")?);
    let motifs = PathBuf::from(get("--motifs")?);
    let policy_rows_path = PathBuf::from(get("--policy-rows")?);
    let policy_config_path = PathBuf::from(get("--policy-config")?);
    let policy_version = get("--policy-version").unwrap_or("phase5-v1");
    let motif_metadata = opts.get("--motif-metadata").map(PathBuf::from);

    // The scan envelope is part of the schema-v3 contract: a shared database
    // that cannot state its own p-value and promoter window cannot be audited.
    let manifest_pvalue = manifest_requested_pvalue(&scan_manifest)?;
    let pvalue = match (opts.get("--pvalue"), manifest_pvalue) {
        (Some(raw), recorded) => {
            let requested: f64 = raw.parse().with_context(|| format!("invalid --pvalue {raw}"))?;
            if let Some(recorded) = recorded {
                if !approximately_equal(requested, recorded) {
                    bail!(
                        "--pvalue {requested} disagrees with the scan manifest ({recorded}); \
                         refusing to record an envelope the scan did not use"
                    );
                }
            }
            requested
        }
        (None, Some(recorded)) => recorded,
        (None, None) => 1e-4,
    };
    let promoter_upstream = parse_usize_option(opts, "--promoter-upstream", 1000)?;
    let promoter_downstream = parse_usize_option(opts, "--promoter-downstream", 0)?;
    let mapping_mode = match opts.get("--mapping-mode") {
        Some(mode) if matches!(mode.as_str(), "strict" | "balanced" | "exploratory") => {
            Some(mode.clone())
        }
        Some(other) => bail!("--mapping-mode must be strict, balanced or exploratory, got {other}"),
        None => None,
    };
    let site_count_mode = opts
        .get("--site-count")
        .map(|raw| SiteCountMode::parse(raw))
        .transpose()?
        .unwrap_or(SiteCountMode::Independent);
    let source_mode = match &mapping_mode {
        Some(mode) => format!("site_only_{mode}"),
        None => "site_only".to_string(),
    };

    if tmp.exists() {
        fs::remove_file(tmp)?;
    }
    let conn = Connection::open(tmp)?;
    conn.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL;")?;
    conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)", [])?;
    conn.execute_batch(HITS_SQL)?;

    let genes = read_tsv(&tf_genes)?.1;
    let (map_headers, map_rows) = read_tsv(&tf_map)?;
    if !map_headers.iter().any(|h| h == "tf_id") || !map_headers.iter().any(|h| h == "motif_id") {
        bail!("tf_map needs tf_id and motif_id");
    }
    let mut mapping_columns = map_headers.clone();
    for required in TF_MOTIF_COLUMNS {
        if !mapping_columns.iter().any(|h| h == required) {
            mapping_columns.push(required.to_string());
        }
    }
    create_text_table(&conn, "tf_motif", &mapping_columns)?;
    insert_rows(&conn, "tf_motif", &mapping_columns, &map_rows)?;

    let coordinate_headers: Vec<String> = mapping_columns.clone();
    create_text_table(&conn, "tf_motif_coordinate", &coordinate_headers)?;
    let coordinate_rows: Vec<BTreeMap<String, String>> = map_rows
        .iter()
        .filter(|row| {
            ["query_full_dbd_start", "query_full_dbd_end", "query_dbd_start", "query_dbd_end"]
                .iter()
                .any(|key| row.get(*key).is_some_and(|value| !value.trim().is_empty()))
        })
        .cloned()
        .collect();
    insert_rows(&conn, "tf_motif_coordinate", &coordinate_headers, &coordinate_rows)?;
    conn.execute("CREATE INDEX IF NOT EXISTS tfm_motif ON tf_motif(motif_id)", [])?;
    conn.execute("CREATE INDEX IF NOT EXISTS tfm_tf ON tf_motif(tf_id)", [])?;
    conn.execute("CREATE INDEX IF NOT EXISTS tmc_tf_motif ON tf_motif_coordinate(tf_id, motif_id)", [])?;

    let mut per_tf_motifs: HashMap<String, BTreeSet<String>> = HashMap::new();
    for row in &map_rows {
        if let (Some(tf), Some(motif)) = (row.get("tf_id"), row.get("motif_id")) {
            per_tf_motifs.entry(tf.clone()).or_default().insert(motif.clone());
        }
    }
    conn.execute(
        "CREATE TABLE tf_catalogue (tf_id TEXT PRIMARY KEY, tf_symbol TEXT, family TEXT, identification_status TEXT, source TEXT, motif_count INTEGER, target_count INTEGER)",
        [],
    )?;
    {
        let mut stmt = conn.prepare("INSERT OR REPLACE INTO tf_catalogue VALUES (?,?,?,?,?,?,0)")?;
        for gene in &genes {
            let id = gene.get("gene_id").cloned().unwrap_or_default();
            if id.is_empty() {
                continue;
            }
            stmt.execute(params![
                id,
                gene.get("tf_symbol").cloned().unwrap_or_default(),
                gene.get("family").cloned().unwrap_or_default(),
                gene.get("identification_status").cloned().unwrap_or_default(),
                gene.get("source").cloned().unwrap_or_default(),
                per_tf_motifs.get(&id).map_or(0, |set| set.len() as i64),
            ])?;
        }
    }
    let policy_count = policy_rows(&conn, &policy_rows_path)?;
    let config_count = policy_config(&conn, &policy_config_path)?;
    conn.execute("CREATE TABLE tf_motif_selected (tf_id TEXT, motif_id TEXT)", [])?;
    conn.execute(
        "INSERT INTO tf_motif_selected SELECT tf_id, motif_id FROM tf_motif_policy WHERE policy='best'",
        [],
    )?;
    conn.execute("CREATE INDEX IF NOT EXISTS tms_motif ON tf_motif_selected(motif_id)", [])?;
    conn.execute("CREATE INDEX IF NOT EXISTS tms_tf ON tf_motif_selected(tf_id)", [])?;

    let testability_count = copy_table(&conn, &testability, "motif_testability", &[
        "motif_id", "pwm_hash", "matrix_length", "information_content", "max_score",
        "min_pvalue", "requested_pvalue", "attainable", "status", "mapped",
        "motif_id_hit_count", "tf_count", "background_mode",
    ])?;
    let promoter_count = copy_table(&conn, &promoter_summary, "promoter", &[
        "gene_id", "chrom", "promoter_start", "promoter_end", "promoter_strand",
        "promoter_length", "scanned_pwms", "testable_pwms", "pwms_with_hits",
        "sites", "best_pvalue", "best_score", "best_pwm_hash", "best_motif_ids",
    ])?;
    conn.execute("CREATE INDEX IF NOT EXISTS motif_testability_id ON motif_testability(motif_id)", [])?;
    conn.execute("CREATE INDEX IF NOT EXISTS promoter_gene ON promoter(gene_id)", [])?;
    let registry_count = registry(&conn, &motifs, motif_metadata.as_deref())?;
    let scan_config_count = scan_config(&conn, &scan_manifest)?;
    let motif_target_count = aggregate_motif_target(&conn, &sites, site_count_mode)?;

    conn.execute(
        "INSERT INTO pairs
         SELECT s.tf_id, COALESCE(c.tf_symbol, s.tf_id), mt.target_id, '',
                SUM(mt.n_sites), COUNT(DISTINCT mt.motif_id), MIN(mt.best_p),
                MAX(mt.best_score), GROUP_CONCAT(DISTINCT mt.motif_id), 'D'
         FROM tf_motif_selected s JOIN motif_target mt ON mt.motif_id=s.motif_id
         LEFT JOIN tf_catalogue c ON c.tf_id=s.tf_id
         GROUP BY s.tf_id, mt.target_id",
        [],
    )?;
    conn.execute(
        "CREATE TEMP TABLE tf_counts AS
         SELECT s.tf_id, COUNT(DISTINCT mt.target_id) n
         FROM tf_motif_selected s JOIN motif_target mt ON mt.motif_id=s.motif_id
         GROUP BY s.tf_id",
        [],
    )?;
    conn.execute(
        "UPDATE tf_catalogue SET target_count=COALESCE((SELECT n FROM tf_counts WHERE tf_counts.tf_id=tf_catalogue.tf_id),0)",
        [],
    )?;
    conn.execute("DROP TABLE tf_counts", [])?;
    for (table, index, column) in [
        ("pairs", "pairs_tf", "tf_id"), ("pairs", "pairs_target", "target_id"),
        ("hits", "hits_tf", "tf_id"), ("hits", "hits_target", "target_id"),
        ("hits", "hits_motif", "motif_id"), ("hits", "hits_p", "p_value"),
    ] {
        conn.execute(&format!("CREATE INDEX IF NOT EXISTS {index} ON {table}({column})"), [])?;
    }

    let tables = [
        "tf_catalogue", "tf_motif", "tf_motif_coordinate", "tf_motif_selected",
        "tf_motif_policy", "tf_motif_policy_config", "motif_target", "motif_registry",
        "motif_testability", "promoter", "scan_config", "pairs", "hits",
    ];
    let mut counts = BTreeMap::new();
    for table in tables {
        counts.insert(table, conn.query_row(&format!("SELECT COUNT(*) FROM {table}"), [], |r| r.get::<_, i64>(0))?);
    }
    let inputs = serde_json::json!({
        "tf_genes": {"file": tf_genes.file_name().map(|v| v.to_string_lossy().into_owned()), "sha256": sha256(&tf_genes)?},
        "tf_map": {"file": tf_map.file_name().map(|v| v.to_string_lossy().into_owned()), "sha256": sha256(&tf_map)?},
        "sites": {"file": sites.file_name().map(|v| v.to_string_lossy().into_owned()), "sha256": sha256(&sites)?},
        "testability": {"file": testability.file_name().map(|v| v.to_string_lossy().into_owned()), "sha256": sha256(&testability)?},
        "promoter_summary": {"file": promoter_summary.file_name().map(|v| v.to_string_lossy().into_owned()), "sha256": sha256(&promoter_summary)?},
        "scan_manifest": {"file": scan_manifest.file_name().map(|v| v.to_string_lossy().into_owned()), "sha256": sha256(&scan_manifest)?},
        "motifs": {"file": motifs.file_name().map(|v| v.to_string_lossy().into_owned()), "sha256": sha256(&motifs)?},
    });
    let mut meta = vec![
        ("schema_version", serde_json::json!("regutrace-regulatory-background-db-v3")),
        ("policy_version", serde_json::json!(policy_version)),
        ("policy", serde_json::json!("best@0.8")),
        ("evidence_level", serde_json::json!("D")),
        ("source_mode", serde_json::json!(source_mode)),
        ("promoter_upstream_bp", serde_json::json!(promoter_upstream)),
        ("promoter_downstream_bp", serde_json::json!(promoter_downstream)),
        ("pvalue", serde_json::json!(pvalue)),
        ("coordinate_system", serde_json::json!(
            "protein coordinates for query DBD; donor coordinates per reference protein"
        )),
        ("normalized_layers", serde_json::json!({
            "scan_config": "scan parameters and background",
            "motif_testability": "attainable p/status per motif",
            "promoter": "per-promoter scanned/testable/site summary",
            "site_table": "external optional attachment (not materialized in this DB)"
        })),
        ("tables", serde_json::json!(tables)),
        ("counts", serde_json::json!(counts)),
        ("inputs", inputs),
        ("limitations", serde_json::json!([
            "All PWM motif-target associations are evidence level D; they are not experimental binding.",
            "pairs is the best-supported policy; all/best/reps associations are stored in tf_motif_policy.",
            "Unobserved genes are not negatives in any downstream benchmark."
        ])),
    ];
    if let Some(mode) = &mapping_mode {
        meta.push(("mapping_mode", serde_json::json!(mode)));
    }
    // Optional key: pre-2026-11 artifacts have no `site_count_mode` and must keep validating, so
    // this is recorded for new builds without entering the required-key contract.
    meta.push(("site_count_mode", serde_json::json!(site_count_mode.as_str())));
    for (key, value) in meta {
        conn.execute(
            "INSERT OR REPLACE INTO meta VALUES (?,?)",
            params![key, if let serde_json::Value::String(value) = &value { value.clone() } else { serde_json::to_string(&value)? }],
        )?;
    }
    conn.execute_batch("PRAGMA user_version=3; PRAGMA wal_checkpoint(TRUNCATE); PRAGMA journal_mode=DELETE;")?;
    conn.execute("VACUUM", [])?;
    drop(conn);
    publish_temp(tmp, out)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({
            "out": out,
            "policy_rows": policy_count,
            "policy_configs": config_count,
            "tf_genes": genes.len(),
            "tf_motif": map_rows.len(),
            "tf_motif_coordinate": coordinate_rows.len(),
            "motif_target": motif_target_count,
            "motif_registry": registry_count,
            "motif_testability": testability_count,
            "promoter": promoter_count,
            "scan_config": scan_config_count,
        }))?
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    type Row<'a> = (&'a str, &'a str, &'a str, &'a str, &'a str, &'a str);

    /// Runs `aggregate_motif_target` over a throwaway sites table.
    /// Row order: motif_ids, gene_id, hit_start, hit_end, score, p_value.
    fn aggregate(rows: &[Row], mode: SiteCountMode, tag: &str) -> Vec<(String, String, i64, f64, f64)> {
        let dir = std::env::temp_dir().join(format!("regutrace_bg_{}_{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("sites.tsv");
        {
            let mut file = std::fs::File::create(&path).unwrap();
            writeln!(file, "motif_ids\tgene_id\thit_start\thit_end\tscore\tp_value").unwrap();
            for (motif, gene, start, end, score, p) in rows {
                writeln!(file, "{motif}\t{gene}\t{start}\t{end}\t{score}\t{p}").unwrap();
            }
        }
        let conn = Connection::open_in_memory().unwrap();
        aggregate_motif_target(&conn, &path, mode).unwrap();
        let mut stmt = conn
            .prepare("SELECT motif_id, target_id, n_sites, best_p, best_score FROM motif_target ORDER BY motif_id, target_id")
            .unwrap();
        let out = stmt
            .query_map([], |row| {
                Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?, row.get(4)?))
            })
            .unwrap()
            .collect::<std::result::Result<Vec<_>, _>>()
            .unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        out
    }

    #[test]
    fn site_count_mode_rejects_unknown_values() {
        assert!(SiteCountMode::parse("nonsense").is_err());
        assert!(SiteCountMode::parse("independent").is_ok());
        assert!(SiteCountMode::parse("overlapping").is_ok());
    }

    #[test]
    fn overlapping_windows_become_one_independent_site() {
        // A degenerate motif sliding one base at a time is one binding event, not three.
        let rows = [
            ("MA1", "G1", "51", "56", "12.0", "1e-6"),
            ("MA1", "G1", "52", "57", "12.0", "1e-6"),
            ("MA1", "G1", "53", "58", "12.0", "1e-6"),
            // Two genuinely separate sites in another promoter stay separate.
            ("MA1", "G2", "51", "56", "12.0", "1e-6"),
            ("MA1", "G2", "60", "65", "12.0", "1e-6"),
        ];
        let independent = aggregate(&rows, SiteCountMode::Independent, "overlap");
        assert_eq!(independent[0], ("MA1".into(), "G1".into(), 1, 1e-6, 12.0));
        assert_eq!(independent[1], ("MA1".into(), "G2".into(), 2, 1e-6, 12.0));

        // The legacy mode is what the frozen Phase 5 artifact used; it must still be exact.
        let legacy = aggregate(&rows, SiteCountMode::Overlapping, "overlap-legacy");
        assert_eq!(legacy[0].2, 3);
        assert_eq!(legacy[1].2, 2);
    }

    #[test]
    fn opposite_strands_at_one_position_count_once() {
        let rows = [
            ("MA1", "G1", "51", "56", "+", "1e-6"),
            ("MA1", "G1", "51", "56", "-", "1e-6"),
        ];
        // The strand column is not part of the count: one position is one site.
        let independent = aggregate(
            &rows.map(|(m, g, s, e, _strand, p)| (m, g, s, e, "12.0", p)),
            SiteCountMode::Independent,
            "strand",
        );
        assert_eq!(independent[0].2, 1);
    }

    #[test]
    fn best_p_and_best_score_are_taken_over_raw_windows() {
        let rows = [
            ("MA1", "G1", "51", "56", "7.5", "1e-3"),
            ("MA1", "G1", "52", "57", "21.9", "1e-7"),
        ];
        let independent = aggregate(&rows, SiteCountMode::Independent, "best");
        let (_, _, sites, best_p, best_score) = independent[0].clone();
        assert_eq!(sites, 1);
        assert!((best_p - 1e-7).abs() < 1e-18);
        assert!((best_score - 21.9).abs() < 1e-9);
    }

    #[test]
    fn shared_pwm_rows_count_for_every_motif_and_blank_coordinates_stand_alone() {
        let rows = [
            // One physical site shared by two motif IDs (bit-identical PWMs).
            ("MA1;MA2", "G1", "51", "56", "12.0", "1e-6"),
            // A site row without usable coordinates cannot be merged with anything.
            ("MA1", "G1", "", "", "9.0", "1e-4"),
        ];
        let independent = aggregate(&rows, SiteCountMode::Independent, "shared");
        assert_eq!(independent[0], ("MA1".into(), "G1".into(), 2, 1e-6, 12.0));
        assert_eq!(independent[1], ("MA2".into(), "G1".into(), 1, 1e-6, 12.0));
    }
}
