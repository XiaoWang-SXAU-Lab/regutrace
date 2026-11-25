//! Minimal, explicitly exploratory motif enrichment over gene sets.
//!
//! The test is a one-sided Fisher exact test on a binary promoter-hit matrix:
//! for each motif, a gene is a "hit" if the unique site table contains at least
//! one site for that motif in its promoter. The universe is the set of scanned
//! promoters from `promoter-summary.tsv`. BH q-values are computed across motifs
//! within one gene set.
//!
//! This is not a calibrated regulatory enrichment model. It ignores promoter
//! length, GC, overlapping sites, motif similarity and gene-set overlap. Those
//! limitations are recorded in `docs/development/phase3-enrichment-design.md`.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;

use anyhow::{bail, Context, Result};
use csv::ReaderBuilder;

use super::parse_options;

#[derive(Clone, Debug)]
struct TestRow {
    gene_set: String,
    description: String,
    set_size: usize,
    universe_size: usize,
    motif_id: String,
    motif_hit_genes: usize,
    overlap: usize,
    expected_overlap: f64,
    fold_enrichment: f64,
    odds_ratio: f64,
    p_value: f64,
    q_value: f64,
    tested_motifs: usize,
}

fn find_col(headers: &csv::StringRecord, names: &[&str]) -> Option<usize> {
    names
        .iter()
        .find_map(|name| headers.iter().position(|header| header == *name))
}

fn read_site_genes(path: &Path) -> Result<BTreeMap<String, BTreeSet<String>>> {
    let file = File::open(path).with_context(|| format!("open sites {}", path.display()))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(BufReader::new(file));
    let headers = reader
        .headers()
        .with_context(|| format!("read sites header {}", path.display()))?
        .clone();
    let motif_col = find_col(&headers, &["motif_ids", "motif_id"])
        .with_context(|| format!("sites table lacks motif_ids: {}", path.display()))?;
    let gene_col = find_col(&headers, &["gene_id", "target_id"])
        .with_context(|| format!("sites table lacks gene_id: {}", path.display()))?;
    let mut out: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for (line, record) in reader.records().enumerate() {
        let record = record.with_context(|| format!("{}: row {}", path.display(), line + 2))?;
        let gene = record.get(gene_col).unwrap_or("").trim();
        if gene.is_empty() {
            continue;
        }
        let motifs = record.get(motif_col).unwrap_or("");
        for motif in motifs.split(';') {
            let motif = motif.trim();
            if !motif.is_empty() {
                out.entry(motif.to_string())
                    .or_default()
                    .insert(gene.to_string());
            }
        }
    }
    if out.is_empty() {
        bail!("sites table has no motif–gene rows: {}", path.display());
    }
    Ok(out)
}

fn read_universe(path: &Path) -> Result<BTreeSet<String>> {
    let file =
        File::open(path).with_context(|| format!("open promoter summary {}", path.display()))?;
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(BufReader::new(file));
    let headers = reader
        .headers()
        .with_context(|| format!("read promoter summary header {}", path.display()))?
        .clone();
    let gene_col = find_col(&headers, &["gene_id", "target_id"])
        .with_context(|| format!("promoter summary lacks gene_id: {}", path.display()))?;
    let mut out = BTreeSet::new();
    for record in reader.records() {
        let record = record?;
        let gene = record.get(gene_col).unwrap_or("").trim();
        if !gene.is_empty() {
            out.insert(gene.to_string());
        }
    }
    if out.is_empty() {
        bail!("promoter summary has no genes: {}", path.display());
    }
    Ok(out)
}

fn read_gmt(path: &Path) -> Result<Vec<(String, String, Vec<String>)>> {
    let file = File::open(path).with_context(|| format!("open GMT {}", path.display()))?;
    let mut out = Vec::new();
    for (line_no, line) in BufReader::new(file).lines().enumerate() {
        let line = line.with_context(|| format!("read GMT line {}", line_no + 1))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            bail!("GMT line {} has fewer than two fields", line_no + 1);
        }
        let name = fields[0].trim().to_string();
        if name.is_empty() {
            bail!("GMT line {} has an empty set name", line_no + 1);
        }
        let (description, genes) = if fields.len() >= 3 {
            (fields[1].to_string(), fields[2..].to_vec())
        } else {
            (String::new(), fields[1..].to_vec())
        };
        let genes: Vec<String> = genes
            .into_iter()
            .map(str::trim)
            .filter(|gene| !gene.is_empty())
            .map(str::to_string)
            .collect();
        out.push((name, description, genes));
    }
    if out.is_empty() {
        bail!("GMT has no gene sets: {}", path.display());
    }
    Ok(out)
}

fn log_choose(n: u64, k: u64) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let n = n as f64;
    let k = k as f64;
    // lgamma is available through the standard library? No: implement via
    // Stirling-independent product would underflow for large n. Use the
    // Lanczos approximation below.
    ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0)
}

/// Lanczos approximation (g=7, n=9), accurate to ~1e-14 for the integer
/// arguments used here.
fn ln_gamma(x: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = COEFFS[0];
    let t = x + 7.5;
    for (index, coefficient) in COEFFS.iter().enumerate().skip(1) {
        a += coefficient / (x + index as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// P(X >= observed) for X ~ Hypergeometric(N, K, n).
fn hypergeometric_sf(observed: u64, draws: u64, successes: u64, population: u64) -> f64 {
    if observed == 0 {
        return 1.0;
    }
    let minimum = draws.saturating_sub(population.saturating_sub(successes));
    let maximum = draws.min(successes);
    if observed > maximum {
        return 0.0;
    }
    let start = observed.max(minimum);
    let log_denominator = log_choose(population, draws);
    let mut probability = (log_choose(successes, start)
        + log_choose(population - successes, draws - start)
        - log_denominator)
        .exp();
    let mut total = probability;
    for x in start..maximum {
        // P(x+1)/P(x). All factors are positive on the support.
        let numerator = (successes - x) as f64 * (draws - x) as f64;
        let denominator =
            (x + 1) as f64 * (population - successes - draws + x + 1) as f64;
        if denominator <= 0.0 || !numerator.is_finite() {
            break;
        }
        probability *= numerator / denominator;
        total += probability;
        if probability == 0.0 {
            break;
        }
    }
    total.clamp(0.0, 1.0)
}

fn bh_qvalues(pvalues: &[f64]) -> Vec<f64> {
    let m = pvalues.len();
    if m == 0 {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|a, b| {
        pvalues[*a]
            .partial_cmp(&pvalues[*b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut q: Vec<f64> = vec![1.0; m];
    let mut running: f64 = 1.0;
    for rank in (0..m).rev() {
        let index = order[rank];
        let value = pvalues[index] * m as f64 / (rank + 1) as f64;
        running = running.min(value);
        q[index] = running.clamp(0.0, 1.0);
    }
    q
}

pub(super) fn run(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    for key in opts.keys() {
        if !matches!(
            key.as_str(),
            "--sites"
                | "--promoter-summary"
                | "--gene-sets"
                | "--out"
                | "--summary-out"
                | "--min-set-size"
                | "--max-set-size"
        ) {
            bail!("unknown option: {key}");
        }
    }
    let sites_path = opts.get("--sites").context("missing --sites")?;
    let summary_path = opts
        .get("--promoter-summary")
        .context("missing --promoter-summary")?;
    let gmt_path = opts.get("--gene-sets").context("missing --gene-sets")?;
    let out_path = opts.get("--out").context("missing --out")?;
    let min_set_size: usize = opts
        .get("--min-set-size")
        .map(|value| value.parse().context("invalid --min-set-size"))
        .transpose()?
        .unwrap_or(5);
    let max_set_size: usize = opts
        .get("--max-set-size")
        .map(|value| value.parse().context("invalid --max-set-size"))
        .transpose()?
        .unwrap_or(usize::MAX);
    if min_set_size == 0 || min_set_size > max_set_size {
        bail!("require 0 < --min-set-size <= --max-set-size");
    }

    let motif_genes = read_site_genes(Path::new(sites_path))?;
    let universe = read_universe(Path::new(summary_path))?;
    let gene_sets = read_gmt(Path::new(gmt_path))?;
    let universe_size = universe.len();
    let motif_hit_genes: BTreeMap<String, BTreeSet<String>> = motif_genes
        .into_iter()
        .map(|(motif, genes)| {
            (
                motif,
                genes
                    .into_iter()
                    .filter(|gene| universe.contains(gene))
                    .collect(),
            )
        })
        .collect();
    let motif_ids: Vec<&String> = motif_hit_genes.keys().collect();

    let mut rows = Vec::new();
    for (set_name, description, genes) in gene_sets {
        let set_genes: BTreeSet<String> = genes
            .into_iter()
            .filter(|gene| universe.contains(gene))
            .collect();
        let set_size = set_genes.len();
        if set_size < min_set_size || set_size > max_set_size {
            continue;
        }
        let mut set_rows = Vec::new();
        for motif in &motif_ids {
            let hits = motif_hit_genes.get(*motif).expect("motif key");
            let overlap = set_genes.intersection(hits).count();
            let motif_hit_count = hits.len();
            let expected = if universe_size > 0 {
                set_size as f64 * motif_hit_count as f64 / universe_size as f64
            } else {
                0.0
            };
            let fold = if expected > 0.0 {
                overlap as f64 / expected
            } else {
                f64::NAN
            };
            let a = overlap as f64;
            let b = (set_size - overlap) as f64;
            let c = (motif_hit_count - overlap) as f64;
            let d = (universe_size - set_size - motif_hit_count + overlap) as f64;
            let odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5));
            let p_value = hypergeometric_sf(
                overlap as u64,
                set_size as u64,
                motif_hit_count as u64,
                universe_size as u64,
            );
            set_rows.push(TestRow {
                gene_set: set_name.clone(),
                description: description.clone(),
                set_size,
                universe_size,
                motif_id: (*motif).clone(),
                motif_hit_genes: motif_hit_count,
                overlap,
                expected_overlap: expected,
                fold_enrichment: fold,
                odds_ratio,
                p_value,
                q_value: 1.0,
                tested_motifs: motif_ids.len(),
            });
        }
        let pvalues: Vec<f64> = set_rows.iter().map(|row| row.p_value).collect();
        let qvalues = bh_qvalues(&pvalues);
        for (row, q) in set_rows.iter_mut().zip(qvalues) {
            row.q_value = q;
        }
        rows.extend(set_rows);
    }
    if rows.is_empty() {
        bail!("no gene set passed the size filters");
    }
    rows.sort_by(|a, b| {
        a.gene_set
            .cmp(&b.gene_set)
            .then(
                a.p_value
                    .partial_cmp(&b.p_value)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
            .then(a.motif_id.cmp(&b.motif_id))
    });

    let file = File::create(out_path).with_context(|| format!("create {}", out_path))?;
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(BufWriter::new(file));
    writer.write_record([
        "gene_set",
        "description",
        "set_size_in_universe",
        "universe_size",
        "motif_id",
        "motif_hit_genes",
        "overlap",
        "expected_overlap",
        "fold_enrichment",
        "odds_ratio",
        "fisher_p",
        "bh_q",
        "tested_motifs_in_set",
    ])?;
    for row in &rows {
        writer.write_record([
            row.gene_set.as_str(),
            row.description.as_str(),
            &row.set_size.to_string(),
            &row.universe_size.to_string(),
            row.motif_id.as_str(),
            &row.motif_hit_genes.to_string(),
            &row.overlap.to_string(),
            &format!("{:.6}", row.expected_overlap),
            &if row.fold_enrichment.is_finite() {
                format!("{:.6}", row.fold_enrichment)
            } else {
                String::new()
            },
            &format!("{:.6}", row.odds_ratio),
            &format!("{:.17e}", row.p_value),
            &format!("{:.17e}", row.q_value),
            &row.tested_motifs.to_string(),
        ])?;
    }
    writer.flush()?;

    if let Some(path) = opts.get("--summary-out") {
        let mut by_set: BTreeMap<&str, (usize, usize, f64, usize)> = BTreeMap::new();
        for row in &rows {
            let entry = by_set
                .entry(row.gene_set.as_str())
                .or_insert((row.set_size, row.tested_motifs, 1.0, 0));
            entry.2 = entry.2.min(row.q_value);
            if row.q_value <= 0.05 {
                entry.3 += 1;
            }
        }
        let file = File::create(path).with_context(|| format!("create {}", path))?;
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(BufWriter::new(file));
        writer.write_record([
            "gene_set",
            "set_size_in_universe",
            "tested_motifs",
            "min_bh_q",
            "motifs_q_le_0.05",
        ])?;
        for (name, (set_size, tested, min_q, significant)) in by_set {
            writer.write_record([
                name,
                &set_size.to_string(),
                &tested.to_string(),
                &format!("{:.17e}", min_q),
                &significant.to_string(),
            ])?;
        }
        writer.flush()?;
    }
    let n_significant = rows.iter().filter(|row| row.q_value <= 0.05).count();
    eprintln!(
        "enrichment: {} motif–gene-set tests; {} with BH q <= 0.05; universe {} genes",
        rows.len(),
        n_significant,
        universe_size
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hypergeometric_matches_known_values() {
        // N=10, K=4, n=3: P(X>=2) = (C(4,2)C(6,1)+C(4,3)C(6,0))/C(10,3)
        let expected = (6.0 * 6.0 + 4.0 * 1.0) / 120.0;
        assert!((hypergeometric_sf(2, 3, 4, 10) - expected).abs() < 1e-12);
        assert!((hypergeometric_sf(3, 3, 4, 10) - 4.0 / 120.0).abs() < 1e-12);
        assert_eq!(hypergeometric_sf(0, 3, 4, 10), 1.0);
    }

    #[test]
    fn bh_qvalues_are_monotone() {
        let q = bh_qvalues(&[0.01, 0.04, 0.03, 0.20]);
        assert!((q[0] - 0.04).abs() < 1e-12);
        assert!((q[1] - 0.04 * 4.0 / 3.0).abs() < 1e-12);
        assert!((q[2] - 0.04 * 4.0 / 3.0).abs() < 1e-12);
        assert!((q[3] - 0.20).abs() < 1e-12);
    }
}
