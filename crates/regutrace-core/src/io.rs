//! TSV input/output for TF–Target evidence records.

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{bail, Context, Result};
use csv::{ReaderBuilder, WriterBuilder};

use crate::fasta::open_reader;

#[derive(Clone, Debug)]
pub struct TfMotif {
    pub tf_id: String,
    pub tf_symbol: String,
    pub family: String,
    pub motif_id: String,
}

#[derive(Clone, Debug, Default)]
pub struct EdgeRecord {
    pub tf_id: String,
    pub tf_symbol: String,
    pub target_id: String,
    pub target_symbol: String,
    pub evidence_type: String,
    pub evidence_level: String,
    pub method: String,
    pub source: String,
    pub pmid: String,
    pub condition: String,
    pub promoter_start: Option<u64>,
    pub promoter_end: Option<u64>,
    pub motif_id: String,
    pub hit_start: Option<u64>,
    pub hit_end: Option<u64>,
    pub strand: String,
    pub score: Option<f64>,
    pub p_value: Option<f64>,
    pub notes: String,
}

fn normalize_header(name: &str) -> String {
    name.trim()
        .trim_start_matches('\u{feff}')
        .to_ascii_lowercase()
        .replace([' ', '-'], "_")
}

fn find_col(headers: &csv::StringRecord, names: &[&str]) -> Option<usize> {
    headers.iter().position(|h| {
        let h = normalize_header(h);
        names.iter().any(|n| h == *n)
    })
}

fn cell(rec: &csv::StringRecord, idx: Option<usize>) -> String {
    idx.and_then(|i| rec.get(i))
        .map(|s| s.trim().to_string())
        .unwrap_or_default()
}

fn parse_level(raw: &str, context: &str) -> Result<String> {
    let level = raw.trim().to_ascii_uppercase();
    if level.is_empty() {
        bail!("{context}: missing evidence_level (A/B/C/D/E)");
    }
    if !matches!(level.as_str(), "A" | "B" | "C" | "D" | "E") {
        bail!("{context}: invalid evidence_level '{raw}' (expected A/B/C/D/E)");
    }
    Ok(level)
}

fn parse_u64_field(raw: String) -> Result<Option<u64>> {
    let t = raw.trim();
    if t.is_empty() {
        return Ok(None);
    }
    Ok(Some(
        t.parse::<u64>()
            .with_context(|| format!("invalid integer '{t}'"))?,
    ))
}

fn parse_f64_field(raw: String) -> Result<Option<f64>> {
    let t = raw.trim();
    if t.is_empty() {
        return Ok(None);
    }
    let v = t
        .parse::<f64>()
        .with_context(|| format!("invalid number '{t}'"))?;
    if !v.is_finite() {
        bail!("non-finite number '{t}'");
    }
    Ok(Some(v))
}

pub fn read_tf_map(path: &Path) -> Result<Vec<TfMotif>> {
    let reader = open_reader(path)?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(reader);
    let headers = rdr.headers().with_context(|| format!("read header {}", path.display()))?.clone();
    let tf = find_col(&headers, &["tf_id", "tf"])
        .with_context(|| format!("tf-map missing tf_id column: {}", path.display()))?;
    let motif = find_col(&headers, &["motif_id", "motif"])
        .with_context(|| format!("tf-map missing motif_id column: {}", path.display()))?;
    let symbol = find_col(&headers, &["tf_symbol", "symbol"]);
    let family = find_col(&headers, &["family", "tf_family"]);
    let mut out = Vec::new();
    let mut row = 1usize;
    for rec in rdr.records() {
        let rec = rec?;
        row += 1;
        let tf_id = cell(&rec, Some(tf));
        let motif_id = cell(&rec, Some(motif));
        if tf_id.is_empty() || motif_id.is_empty() {
            bail!("tf-map row {row}: empty tf_id or motif_id");
        }
        out.push(TfMotif {
            tf_symbol: if symbol.is_some() {
                cell(&rec, symbol)
            } else {
                tf_id.clone()
            },
            family: cell(&rec, family),
            tf_id,
            motif_id,
        });
    }
    if out.is_empty() {
        bail!("tf-map has no usable rows: {}", path.display());
    }
    Ok(out)
}

pub fn read_annotation(path: &Path) -> Result<HashMap<String, String>> {
    let reader = open_reader(path)?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(reader);
    let headers = rdr.headers()?.clone();
    let gene = find_col(&headers, &["gene_id", "gene", "id"])
        .with_context(|| format!("annotation missing gene_id column: {}", path.display()))?;
    let symbol = find_col(&headers, &["symbol", "gene_symbol", "name"]);
    let mut out = HashMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        let id = cell(&rec, Some(gene));
        if id.is_empty() {
            continue;
        }
        let sym = cell(&rec, symbol);
        if !sym.is_empty() {
            out.insert(id, sym);
        }
    }
    Ok(out)
}

pub fn read_user_list(path: &Path) -> Result<Vec<EdgeRecord>> {
    let reader = open_reader(path)?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(reader);
    let headers = rdr.headers()?.clone();
    let tf = find_col(&headers, &["tf", "tf_id"])
        .with_context(|| format!("user list missing tf/tf_id column: {}", path.display()))?;
    let target = find_col(&headers, &["target", "target_id"])
        .with_context(|| format!("user list missing target/target_id column: {}", path.display()))?;
    let tf_symbol = find_col(&headers, &["tf_symbol"]);
    let target_symbol = find_col(&headers, &["target_symbol"]);
    let evidence_type = find_col(&headers, &["evidence_type"]);
    let evidence_level = find_col(&headers, &["evidence_level", "evidence"]);
    let method = find_col(&headers, &["method"]);
    let source = find_col(&headers, &["source"]);
    let pmid = find_col(&headers, &["pmid", "pubmed"]);
    let condition = find_col(&headers, &["condition"]);
    let score = find_col(&headers, &["score"]);
    let notes = find_col(&headers, &["notes", "note"]);

    let mut out = Vec::new();
    let mut row = 1usize;
    for rec in rdr.records() {
        let rec = rec?;
        row += 1;
        let tf_id = cell(&rec, Some(tf));
        let target_id = cell(&rec, Some(target));
        if tf_id.is_empty() || target_id.is_empty() {
            bail!("user list row {row}: empty tf or target");
        }
        let ev_type = {
            let v = cell(&rec, evidence_type);
            if v.is_empty() { "user_curated".to_string() } else { v }
        };
        let ev_level = parse_level(
            &cell(&rec, evidence_level),
            &format!("user list row {row}"),
        )?;
        let score_v = parse_f64_field(cell(&rec, score))
            .with_context(|| format!("user list row {row}"))?;
        out.push(EdgeRecord {
            tf_symbol: cell(&rec, tf_symbol),
            target_symbol: cell(&rec, target_symbol),
            evidence_type: ev_type,
            evidence_level: ev_level,
            method: cell(&rec, method),
            source: cell(&rec, source),
            pmid: cell(&rec, pmid),
            condition: cell(&rec, condition),
            score: score_v,
            notes: cell(&rec, notes),
            tf_id,
            target_id,
            ..Default::default()
        });
    }
    Ok(out)
}

fn fmt_opt<T: std::fmt::Display>(v: Option<T>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}
fn fmt_f64(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.6}")).unwrap_or_default()
}
fn fmt_p(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.6e}")).unwrap_or_default()
}

fn write_header<W: std::io::Write>(wtr: &mut csv::Writer<W>) -> Result<()> {
    wtr.write_record([
        "tf_id", "tf_symbol", "target_id", "target_symbol", "evidence_type", "evidence_level",
        "method", "source", "pmid", "condition", "promoter_start", "promoter_end", "motif_id",
        "hit_start", "hit_end", "strand", "score", "p_value", "notes",
    ])?;
    Ok(())
}

fn write_edge_record<W: std::io::Write>(
    wtr: &mut csv::Writer<W>,
    e: &EdgeRecord,
) -> Result<()> {
    let promoter_start = fmt_opt(e.promoter_start);
    let promoter_end = fmt_opt(e.promoter_end);
    let hit_start = fmt_opt(e.hit_start);
    let hit_end = fmt_opt(e.hit_end);
    let score = fmt_f64(e.score);
    let p_value = fmt_p(e.p_value);
    wtr.write_record([
        e.tf_id.as_str(),
        e.tf_symbol.as_str(),
        e.target_id.as_str(),
        e.target_symbol.as_str(),
        e.evidence_type.as_str(),
        e.evidence_level.as_str(),
        e.method.as_str(),
        e.source.as_str(),
        e.pmid.as_str(),
        e.condition.as_str(),
        promoter_start.as_str(),
        promoter_end.as_str(),
        e.motif_id.as_str(),
        hit_start.as_str(),
        hit_end.as_str(),
        e.strand.as_str(),
        score.as_str(),
        p_value.as_str(),
        e.notes.as_str(),
    ])?;
    Ok(())
}

/// Incremental evidence writer used for whole-genome scans that can emit
/// millions of rows. Rows are written in scan order rather than sorted, so the
/// caller should only use this when streaming output matters.
pub struct EdgeStream {
    wtr: csv::Writer<BufWriter<File>>,
}

impl EdgeStream {
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
        let mut wtr = WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(BufWriter::new(file));
        write_header(&mut wtr)?;
        Ok(Self { wtr })
    }

    pub fn write(&mut self, e: &EdgeRecord) -> Result<()> {
        write_edge_record(&mut self.wtr, e)
    }

    pub fn flush(&mut self) -> Result<()> {
        self.wtr.flush()?;
        Ok(())
    }
}

pub fn write_edges(path: &Path, edges: &mut Vec<EdgeRecord>) -> Result<()> {
    edges.sort_by(|a, b| {
        a.tf_id
            .cmp(&b.tf_id)
            .then(a.target_id.cmp(&b.target_id))
            .then(a.evidence_level.cmp(&b.evidence_level))
            .then(a.motif_id.cmp(&b.motif_id))
            .then(a.promoter_start.cmp(&b.promoter_start))
            .then(a.hit_start.cmp(&b.hit_start))
            .then(a.strand.cmp(&b.strand))
    });
    let mut wtr = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("create {}", path.display()))?;
    write_header(&mut wtr)?;
    for e in edges {
        write_edge_record(&mut wtr, e)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Read the unified evidence TSV produced by `write_edges` (or by the scanner).
pub fn read_edges(path: &Path) -> Result<Vec<EdgeRecord>> {
    let reader = open_reader(path)?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(reader);
    let headers = rdr.headers()?.clone();
    let tf = find_col(&headers, &["tf_id", "tf"]).context("evidence TSV missing tf_id")?;
    let target =
        find_col(&headers, &["target_id", "target"]).context("evidence TSV missing target_id")?;
    let tf_symbol = find_col(&headers, &["tf_symbol"]);
    let target_symbol = find_col(&headers, &["target_symbol"]);
    let evidence_type = find_col(&headers, &["evidence_type"]);
    let evidence_level = find_col(&headers, &["evidence_level", "evidence"]);
    let method = find_col(&headers, &["method"]);
    let source = find_col(&headers, &["source"]);
    let pmid = find_col(&headers, &["pmid", "pubmed"]);
    let condition = find_col(&headers, &["condition"]);
    let promoter_start = find_col(&headers, &["promoter_start"]);
    let promoter_end = find_col(&headers, &["promoter_end"]);
    let motif_id = find_col(&headers, &["motif_id", "motif"]);
    let hit_start = find_col(&headers, &["hit_start"]);
    let hit_end = find_col(&headers, &["hit_end"]);
    let strand = find_col(&headers, &["strand"]);
    let score = find_col(&headers, &["score"]);
    let p_value = find_col(&headers, &["p_value", "pvalue"]);
    let notes = find_col(&headers, &["notes", "note"]);

    let mut out = Vec::new();
    let mut row = 1usize;
    for rec in rdr.records() {
        let rec = rec?;
        row += 1;
        let ctx = || format!("evidence TSV row {row}");
        let tf_id = cell(&rec, Some(tf));
        let target_id = cell(&rec, Some(target));
        if tf_id.is_empty() || target_id.is_empty() {
            bail!("{}: empty tf_id or target_id", ctx());
        }
        let ev_type = {
            let v = cell(&rec, evidence_type);
            if v.is_empty() {
                bail!("{}: missing evidence_type", ctx());
            } else {
                v
            }
        };
        let ev_level = parse_level(&cell(&rec, evidence_level), &ctx())?;
        out.push(EdgeRecord {
            tf_id,
            tf_symbol: cell(&rec, tf_symbol),
            target_id,
            target_symbol: cell(&rec, target_symbol),
            evidence_type: ev_type,
            evidence_level: ev_level,
            method: cell(&rec, method),
            source: cell(&rec, source),
            pmid: cell(&rec, pmid),
            condition: cell(&rec, condition),
            promoter_start: parse_u64_field(cell(&rec, promoter_start))
                .with_context(ctx)?,
            promoter_end: parse_u64_field(cell(&rec, promoter_end)).with_context(ctx)?,
            motif_id: cell(&rec, motif_id),
            hit_start: parse_u64_field(cell(&rec, hit_start)).with_context(ctx)?,
            hit_end: parse_u64_field(cell(&rec, hit_end)).with_context(ctx)?,
            strand: cell(&rec, strand),
            score: parse_f64_field(cell(&rec, score)).with_context(ctx)?,
            p_value: parse_f64_field(cell(&rec, p_value)).with_context(ctx)?,
            notes: cell(&rec, notes),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn reads_tf_map_and_user_list() {
        let dir = std::env::temp_dir();
        let map = dir.join(format!("regutrace_map_{}.tsv", std::process::id()));
        let user = dir.join(format!("regutrace_user_{}.tsv", std::process::id()));
        {
            let mut f = std::fs::File::create(&map).unwrap();
            writeln!(f, "tf_id\ttf_symbol\tfamily\tmotif_id").unwrap();
            writeln!(f, "AT1G00010\tTF1\tHSF\tMA0001.1").unwrap();
        }
        {
            let mut f = std::fs::File::create(&user).unwrap();
            writeln!(f, "tf\ttarget\tevidence_level\tevidence_type\tsource").unwrap();
            writeln!(f, "AT1G00010\tAT2G00020\tA\tpoint_validation_emsa\tPMID:1").unwrap();
        }
        let m = read_tf_map(&map).unwrap();
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].motif_id, "MA0001.1");
        let u = read_user_list(&user).unwrap();
        assert_eq!(u.len(), 1);
        assert_eq!(u[0].evidence_level, "A");
        let _ = std::fs::remove_file(map);
        let _ = std::fs::remove_file(user);
    }

    #[test]
    fn rejects_missing_tf_map_ids() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_tfmap_missing_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "tf_id\tmotif_id").unwrap();
            writeln!(f, "\tMA0001.1").unwrap();
        }
        let err = read_tf_map(&path).unwrap_err().to_string();
        assert!(err.contains("row 2"), "{err}");
        assert!(err.contains("empty tf_id or motif_id"), "{err}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_invalid_user_evidence_level() {
        let user = std::env::temp_dir().join(format!(
            "regutrace_badlevel_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&user).unwrap();
            writeln!(f, "tf\ttarget\tevidence_level").unwrap();
            writeln!(f, "AT1G00010\tAT2G00020\tuser_curated").unwrap();
        }
        let err = read_user_list(&user).unwrap_err().to_string();
        assert!(err.contains("invalid evidence_level"), "{err}");
        let _ = std::fs::remove_file(user);
    }

    #[test]
    fn rejects_missing_user_evidence_level() {
        let user = std::env::temp_dir().join(format!(
            "regutrace_misslevel_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&user).unwrap();
            writeln!(f, "tf\ttarget").unwrap();
            writeln!(f, "AT1G00010\tAT2G00020").unwrap();
        }
        let err = read_user_list(&user).unwrap_err().to_string();
        assert!(err.contains("missing evidence_level"), "{err}");
        let _ = std::fs::remove_file(user);
    }

    #[test]
    fn rejects_missing_level_in_evidence_tsv() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_evidence_nolevel_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "tf_id\ttarget_id\tevidence_type").unwrap();
            writeln!(f, "TF1\tG1\tprediction_motif").unwrap();
        }
        let err = read_edges(&path).unwrap_err().to_string();
        assert!(err.contains("missing evidence_level"), "{err}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_missing_evidence_target() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_evidence_notarget_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "tf_id\ttarget_id\tevidence_type\tevidence_level").unwrap();
            writeln!(f, "TF\t\tcurated\tA").unwrap();
        }
        let err = read_edges(&path).unwrap_err().to_string();
        assert!(err.contains("empty tf_id or target_id"), "{err}");
        assert!(err.contains("row 2"), "{err}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_missing_user_target() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_user_notarget_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "tf\ttarget\tevidence_level").unwrap();
            writeln!(f, "TF\t\tA").unwrap();
        }
        let err = read_user_list(&path).unwrap_err().to_string();
        assert!(err.contains("empty tf or target"), "{err}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_missing_evidence_type() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_evidence_notype_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "tf_id\ttarget_id\tevidence_level").unwrap();
            writeln!(f, "TF1\tG1\tD").unwrap();
        }
        let err = read_edges(&path).unwrap_err().to_string();
        assert!(err.contains("missing evidence_type"), "{err}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_invalid_evidence_coordinates() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_evidence_badcoord_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(
                f,
                "tf_id\ttarget_id\tevidence_type\tevidence_level\tpromoter_start"
            )
            .unwrap();
            writeln!(f, "TF1\tG1\tprediction_motif\tD\tabc").unwrap();
        }
        let err = format!("{:#}", read_edges(&path).unwrap_err());
        assert!(err.contains("invalid integer"), "{err}");
        assert!(err.contains("row 2"), "{err}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_invalid_user_score() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_user_badscore_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "tf\ttarget\tevidence_level\tscore").unwrap();
            writeln!(f, "TF1\tG1\tA\tabc").unwrap();
        }
        let err = format!("{:#}", read_user_list(&path).unwrap_err());
        assert!(err.contains("invalid number"), "{err}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn round_trips_tiny_p_value() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_pvalue_{}.tsv",
            std::process::id()
        ));
        let mut edges = vec![EdgeRecord {
            tf_id: "TF1".into(),
            target_id: "G1".into(),
            evidence_type: "prediction_motif".into(),
            evidence_level: "D".into(),
            p_value: Some(1.23e-9),
            ..Default::default()
        }];
        write_edges(&path, &mut edges).unwrap();
        let back = read_edges(&path).unwrap();
        let p = back[0].p_value.unwrap();
        assert!((p - 1.23e-9).abs() < 1e-15, "p={p}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn normalizes_lowercase_user_evidence_level() {
        let user = std::env::temp_dir().join(format!(
            "regutrace_lowerlevel_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&user).unwrap();
            writeln!(f, "tf\ttarget\tevidence_level").unwrap();
            writeln!(f, "AT1G00010\tAT2G00020\ta").unwrap();
        }
        let u = read_user_list(&user).unwrap();
        assert_eq!(u[0].evidence_level, "A");
        let _ = std::fs::remove_file(user);
    }
}
