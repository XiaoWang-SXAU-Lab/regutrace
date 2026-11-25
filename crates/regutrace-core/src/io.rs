//! TSV input/output for TF–Target evidence records.

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{bail, Context, Result};
use csv::{ReaderBuilder, WriterBuilder};

use crate::fasta::open_reader;

#[derive(Clone, Debug, Default)]
pub struct TfMotif {
    pub tf_id: String,
    pub tf_symbol: String,
    pub family: String,
    pub motif_id: String,
    pub reference_tf_id: String,
    pub mapping_method: String,
    pub domain_mode: String,
    pub extraction_family: String,
    pub extraction_models: String,
    pub domain_count: String,
}

/// One row of a species TF catalogue. `gene_id` is the stable join key used by
/// mapped scans; motif-header names are never treated as gene identities.
#[derive(Clone, Debug)]
pub struct TfGene {
    pub gene_id: String,
    pub tf_symbol: String,
    pub family: String,
    pub status: String,
    pub source: String,
}

impl TfGene {
    pub fn is_identified(&self) -> bool {
        self.status == "identified"
    }
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
    let mappings = read_tf_map_allow_empty(path)?;
    if mappings.is_empty() {
        bail!("tf-map has no usable rows: {}", path.display());
    }
    Ok(mappings)
}

pub fn read_tf_map_allow_empty(path: &Path) -> Result<Vec<TfMotif>> {
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
    let reference_tf = find_col(&headers, &["reference_tf_id"]);
    let mapping_method = find_col(&headers, &["mapping_method"]);
    let domain_mode = find_col(&headers, &["domain_mode"]);
    let extraction_family = find_col(&headers, &["extraction_family"]);
    let extraction_models = find_col(&headers, &["extraction_models"]);
    let domain_count = find_col(&headers, &["domain_count"]);
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
            reference_tf_id: cell(&rec, reference_tf),
            mapping_method: cell(&rec, mapping_method),
            domain_mode: cell(&rec, domain_mode),
            extraction_family: cell(&rec, extraction_family),
            extraction_models: cell(&rec, extraction_models),
            domain_count: cell(&rec, domain_count),
        });
    }
    Ok(out)
}

pub fn read_tf_genes(path: &Path) -> Result<HashMap<String, TfGene>> {
    let genes = read_tf_genes_allow_empty(path)?;
    if genes.is_empty() {
        bail!("tf-genes has no usable rows: {}", path.display());
    }
    Ok(genes)
}

pub fn read_tf_genes_allow_empty(path: &Path) -> Result<HashMap<String, TfGene>> {
    let reader = open_reader(path)?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(reader);
    let headers = rdr
        .headers()
        .with_context(|| format!("read header {}", path.display()))?
        .clone();
    let gene = find_col(&headers, &["gene_id", "gene", "id"])
        .with_context(|| format!("tf-genes missing gene_id column: {}", path.display()))?;
    let symbol = find_col(&headers, &["tf_symbol", "symbol", "name"]);
    let family = find_col(&headers, &["family", "tf_family"]);
    let status = find_col(&headers, &["identification_status", "status"]);
    let source = find_col(&headers, &["source", "source_release", "catalogue"]);
    let mut out = HashMap::new();
    let mut row = 1usize;
    for rec in rdr.records() {
        let rec = rec?;
        row += 1;
        let gene_id = cell(&rec, Some(gene));
        if gene_id.is_empty() {
            bail!("tf-genes row {row}: empty gene_id");
        }
        let status = if status.is_some() {
            let s = cell(&rec, status).to_ascii_lowercase();
            if s.is_empty() {
                "identified".to_string()
            } else {
                s
            }
        } else {
            "identified".to_string()
        };
        if !matches!(status.as_str(), "identified" | "ambiguous" | "rejected") {
            bail!(
                "tf-genes row {row}: invalid identification_status '{status}' (expected identified/ambiguous/rejected)"
            );
        }
        let entry = TfGene {
            tf_symbol: if symbol.is_some() {
                let s = cell(&rec, symbol);
                if s.is_empty() { gene_id.clone() } else { s }
            } else {
                gene_id.clone()
            },
            family: cell(&rec, family),
            status,
            source: cell(&rec, source),
            gene_id: gene_id.clone(),
        };
        if out.insert(gene_id.clone(), entry).is_some() {
            bail!("tf-genes duplicate gene_id: {gene_id}");
        }
    }
    Ok(out)
}

fn strip_numeric_version(id: &str) -> Option<&str> {
    let (base, version) = id.rsplit_once('.')?;
    if base.is_empty() || version.is_empty() || !version.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    Some(base)
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
    let mut exact = HashMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        let id = cell(&rec, Some(gene));
        if id.is_empty() || id.starts_with('#') {
            continue;
        }
        let sym = cell(&rec, symbol);
        if !sym.is_empty() && sym != "-" {
            exact.insert(id, sym);
        }
    }
    // EggNOG often annotates protein/transcript IDs (Fta_1G0000020.1) while
    // the GFF gene model uses gene IDs (Fta_1G0000020). Exact keys win; add
    // gene-level aliases only when the gene key is absent.
    let mut out = exact.clone();
    for (id, sym) in &exact {
        let mut base = id.as_str();
        while let Some(next) = strip_numeric_version(base) {
            if exact.contains_key(next) {
                break;
            }
            out.entry(next.to_owned()).or_insert_with(|| sym.clone());
            base = next;
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

    /// Open an existing edge table for appending without writing a second
    /// header. Used when a streamed scan is followed by user-curated rows.
    pub fn append(path: &Path) -> Result<Self> {
        let file = std::fs::OpenOptions::new()
            .append(true)
            .open(path)
            .with_context(|| format!("open {} for append", path.display()))?;
        let wtr = WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(BufWriter::new(file));
        Ok(Self { wtr })
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
            writeln!(f, "tf_id\ttf_symbol\tfamily\tmotif_id\tdomain_mode\textraction_family\textraction_models\tdomain_count").unwrap();
            writeln!(f, "AT1G00010\tTF1\tHSF\tMA0001.1\tfull_architecture\tHSF\tHSF_DNA-bind\t1").unwrap();
        }
        {
            let mut f = std::fs::File::create(&user).unwrap();
            writeln!(f, "tf\ttarget\tevidence_level\tevidence_type\tsource").unwrap();
            writeln!(f, "AT1G00010\tAT2G00020\tA\tpoint_validation_emsa\tPMID:1").unwrap();
        }
        let m = read_tf_map(&map).unwrap();
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].motif_id, "MA0001.1");
        assert_eq!(m[0].domain_mode, "full_architecture");
        assert_eq!(m[0].extraction_models, "HSF_DNA-bind");
        assert_eq!(m[0].domain_count, "1");
        let u = read_user_list(&user).unwrap();
        assert_eq!(u.len(), 1);
        assert_eq!(u[0].evidence_level, "A");
        let _ = std::fs::remove_file(map);
        let _ = std::fs::remove_file(user);
    }

    #[test]
    fn reads_annotation_and_adds_gene_level_aliases() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_annotation_alias_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "gene_id\tsymbol\tdescription").unwrap();
            writeln!(f, "Fta_1G0000020.1\tFta_1G0000020\tPhagocyte signaling-impaired protein").unwrap();
            writeln!(f, "Fta_1G0000030.1\tPHOT\tphototropin").unwrap();
            writeln!(f, "G1\tEXACT\texact gene").unwrap();
            writeln!(f, "G1.1\tISOFORM\tisoform").unwrap();
            writeln!(f, "## 3 queries scanned\t\t").unwrap();
        }
        let ann = read_annotation(&path).unwrap();
        assert_eq!(ann["Fta_1G0000020"], "Fta_1G0000020");
        assert_eq!(ann["Fta_1G0000030"], "PHOT");
        assert_eq!(ann["G1"], "EXACT");
        assert!(!ann.contains_key("## 3 queries scanned"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn reads_tf_genes_and_rejects_duplicates() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_tfgenes_{}.tsv",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "gene_id\ttf_symbol\tfamily\tidentification_status\tsource").unwrap();
            writeln!(f, "AT1G01010\tNAC001\tNAC\tidentified\ttest").unwrap();
            writeln!(f, "AT1G01020\tAT1G01020\t\tambiguous\ttest").unwrap();
        }
        let genes = read_tf_genes(&path).unwrap();
        assert_eq!(genes.len(), 2);
        assert!(genes["AT1G01010"].is_identified());
        assert!(!genes["AT1G01020"].is_identified());
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "gene_id").unwrap();
            writeln!(f, "AT1G01010").unwrap();
            writeln!(f, "AT1G01010").unwrap();
        }
        let err = read_tf_genes(&path).unwrap_err().to_string();
        assert!(err.contains("duplicate gene_id"), "{err}");
        let _ = std::fs::remove_file(path);
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
