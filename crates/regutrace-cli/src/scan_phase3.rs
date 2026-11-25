//! Phase 3 scan architecture: unique PWM groups, exact reuse, optional
//! persistent scan cache, empirical zero-order background, and a separate
//! motif-site layer.
//!
//! The legacy `run_scan` path is intentionally kept byte-for-byte compatible.
//! This module is only entered when a Phase 3 option is requested.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use csv::WriterBuilder;
use rayon::prelude::*;
use regutrace_core::{io, promoter, pwm, scan};
use sha2::{Digest, Sha256};

use super::{Args, ProgressWriter};

const SCAN_ALGORITHM: &str = "regutrace-scan-v1";
const SCAN_SCALE: f64 = 10.0;
const CACHE_HEADER: &str = "# regutrace-scan-cache-v1";
const CACHE_VERSION: &str = "regutrace-scan-cache-v1";

#[derive(Clone, Debug)]
pub(crate) struct BackgroundModel {
    pub values: [f64; 4],
    pub mode: String,
    pub raw_counts: [u64; 4],
    pub ambiguous_bases: u64,
    pub pseudocount: f64,
}

#[derive(Clone, Debug)]
struct PwmGroup {
    digest: String,
    motif_ids: Vec<String>,
    pwm: pwm::Pwm,
    scan_key: String,
    min_pvalue: f64,
    tf_count: usize,
}

#[derive(Clone, Debug, Default)]
struct PromoterStats {
    sites: u64,
    pwms_with_hits: u32,
    best_pvalue: f64,
    best_score: f64,
    best_pwm_hash: String,
    best_motif_ids: String,
}

impl PromoterStats {
    fn new() -> Self {
        Self {
            sites: 0,
            pwms_with_hits: 0,
            best_pvalue: 1.0,
            best_score: f64::NEG_INFINITY,
            best_pwm_hash: String::new(),
            best_motif_ids: String::new(),
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Phase3Stats {
    pub unique_site_hits: u64,
    pub motif_id_hits: u64,
    pub evidence_rows: u64,
    pub pwm_groups: usize,
    pub mapped_motif_ids: usize,
    pub unique_pwms: usize,
    pub duplicate_pwm_groups: usize,
    pub cache_hits: u64,
    pub cache_misses: u64,
    pub testable_groups: u64,
    pub executed_promoter_scans: u64,
    pub skipped_unattainable_scans: u64,
    pub testable_motifs: usize,
    pub unattainable_motifs: usize,
    pub motifs_with_hits: usize,
    pub motifs_without_hits: usize,
    pub background: BackgroundModel,
    pub sites_path: Option<PathBuf>,
    pub testability_path: Option<PathBuf>,
    pub promoter_summary_path: Option<PathBuf>,
    pub manifest_path: Option<PathBuf>,
}

struct ScanKeyBase {
    prefix: Sha256,
}

impl ScanKeyBase {
    fn new(
        background: &[f64; 4],
        pvalue: f64,
        min_score: Option<f64>,
        both_strands: bool,
        promoter_hash: &str,
    ) -> Self {
        let mut prefix = Sha256::new();
        prefix.update(b"regutrace-scan-key-v2\0");
        prefix.update(SCAN_ALGORITHM.as_bytes());
        prefix.update(SCAN_SCALE.to_bits().to_le_bytes());
        for value in background {
            prefix.update(value.to_bits().to_le_bytes());
        }
        prefix.update(pvalue.to_bits().to_le_bytes());
        prefix.update(match min_score {
            Some(value) => {
                let mut buf = [0_u8; 9];
                buf[0] = 1;
                buf[1..].copy_from_slice(&value.to_bits().to_le_bytes());
                buf
            }
            None => [0_u8; 9],
        });
        prefix.update([u8::from(both_strands)]);
        prefix.update(promoter_hash.as_bytes());
        Self { prefix }
    }

    fn key(&self, pwm_digest: &str) -> String {
        let mut hasher = self.prefix.clone();
        hasher.update(pwm_digest.as_bytes());
        format!("{:x}", hasher.finalize())
    }
}

struct ScanCache {
    dir: PathBuf,
    hits: AtomicU64,
    misses: AtomicU64,
}

impl ScanCache {
    fn open(dir: &Path) -> Result<Self> {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("create scan cache directory {}", dir.display()))?;
        Ok(Self {
            dir: dir.to_path_buf(),
            hits: AtomicU64::new(0),
            misses: AtomicU64::new(0),
        })
    }

    fn path(&self, key: &str) -> PathBuf {
        self.dir.join(format!("{key}.tsv"))
    }

    fn load(&self, key: &str, promoter_count: usize) -> Option<Vec<Vec<scan::RawHit>>> {
        let path = self.path(key);
        if !path.exists() {
            self.misses.fetch_add(1, Ordering::Relaxed);
            return None;
        }
        match read_cache_file(&path, promoter_count) {
            Ok(hits) => {
                self.hits.fetch_add(1, Ordering::Relaxed);
                Some(hits)
            }
            Err(error) => {
                eprintln!(
                    "warning: ignoring invalid scan cache {}: {error:#}",
                    path.display()
                );
                self.misses.fetch_add(1, Ordering::Relaxed);
                None
            }
        }
    }

    fn store(&self, key: &str, hits: &[Vec<scan::RawHit>]) -> Result<()> {
        let path = self.path(key);
        let tmp = path.with_extension("tsv.tmp");
        let _ = std::fs::remove_file(&tmp);
        let file = File::create(&tmp)
            .with_context(|| format!("create scan cache temp {}", tmp.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "{CACHE_HEADER}")?;
        writeln!(
            writer,
            "promoter_index\tstart\tend\tstrand\tscore\tp_value"
        )?;
        for (promoter_index, promoter_hits) in hits.iter().enumerate() {
            for hit in promoter_hits {
                writeln!(
                    writer,
                    "{promoter_index}\t{}\t{}\t{}\t{:.17e}\t{:.17e}",
                    hit.start, hit.end, hit.strand, hit.score, hit.p_value
                )?;
            }
        }
        writer.flush()?;
        if path.exists() {
            std::fs::remove_file(&path)
                .with_context(|| format!("replace stale scan cache {}", path.display()))?;
        }
        std::fs::rename(&tmp, &path).with_context(|| {
            format!(
                "commit scan cache {} -> {}",
                tmp.display(),
                path.display()
            )
        })?;
        Ok(())
    }

    fn hit_count(&self) -> u64 {
        self.hits.load(Ordering::Relaxed)
    }

    fn miss_count(&self) -> u64 {
        self.misses.load(Ordering::Relaxed)
    }
}

fn read_cache_file(path: &Path, promoter_count: usize) -> Result<Vec<Vec<scan::RawHit>>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let header = lines
        .next()
        .transpose()?
        .with_context(|| format!("empty scan cache {}", path.display()))?;
    if header.trim() != CACHE_HEADER {
        bail!(
            "unexpected cache header in {}: {:?} (expected {})",
            path.display(),
            header,
            CACHE_HEADER
        );
    }
    let column_header = lines
        .next()
        .transpose()?
        .with_context(|| format!("scan cache {} has no column header", path.display()))?;
    if column_header
        != "promoter_index\tstart\tend\tstrand\tscore\tp_value"
    {
        bail!("unexpected cache columns in {}", path.display());
    }
    let mut out = vec![Vec::new(); promoter_count];
    for (line_no, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 6 {
            bail!(
                "{}:{}: expected 6 columns, found {}",
                path.display(),
                line_no + 3,
                fields.len()
            );
        }
        let promoter_index: usize = fields[0]
            .parse()
            .with_context(|| format!("{}:{}: promoter_index", path.display(), line_no + 3))?;
        if promoter_index >= promoter_count {
            bail!(
                "{}:{}: promoter_index {} outside 0..{}",
                path.display(),
                line_no + 3,
                promoter_index,
                promoter_count
            );
        }
        let start: u64 = fields[1]
            .parse()
            .with_context(|| format!("{}:{}: start", path.display(), line_no + 3))?;
        let end: u64 = fields[2]
            .parse()
            .with_context(|| format!("{}:{}: end", path.display(), line_no + 3))?;
        let strand = match fields[3] {
            "+" => '+',
            "-" => '-',
            other => bail!(
                "{}:{}: invalid strand {:?}",
                path.display(),
                line_no + 3,
                other
            ),
        };
        let score: f64 = fields[4]
            .parse()
            .with_context(|| format!("{}:{}: score", path.display(), line_no + 3))?;
        let p_value: f64 = fields[5]
            .parse()
            .with_context(|| format!("{}:{}: p_value", path.display(), line_no + 3))?;
        if start == 0 || end < start || !score.is_finite() || !(0.0..=1.0).contains(&p_value) {
            bail!(
                "{}:{}: invalid hit geometry or p-value",
                path.display(),
                line_no + 3
            );
        }
        out[promoter_index].push(scan::RawHit {
            start,
            end,
            strand,
            score,
            p_value,
        });
    }
    for hits in &mut out {
        hits.sort_by(|a, b| {
            a.start
                .cmp(&b.start)
                .then(a.end.cmp(&b.end))
                .then(a.strand.cmp(&b.strand))
        });
    }
    Ok(out)
}

fn base_index(byte: u8) -> Option<usize> {
    match byte {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

pub(crate) fn background_model(
    args: &Args,
    promoters: &[promoter::Promoter],
) -> Result<BackgroundModel> {
    match args.background_mode.as_str() {
        "uniform" => Ok(BackgroundModel {
            values: args.background,
            mode: "uniform".to_string(),
            raw_counts: [0; 4],
            ambiguous_bases: 0,
            pseudocount: 0.0,
        }),
        "empirical" => {
            let mut counts = [0_u64; 4];
            let mut ambiguous = 0_u64;
            for promoter in promoters {
                for &byte in &promoter.seq {
                    match base_index(byte) {
                        Some(index) => counts[index] += 1,
                        None => ambiguous += 1,
                    }
                }
            }
            let total: u64 = counts.iter().sum();
            if total == 0 {
                bail!("empirical background requested but promoters contain no A/C/G/T bases");
            }
            // Strand-symmetric zero-order model. Symmetrizing A/T and C/G
            // prevents the scan from depending on which strand a promoter was
            // extracted from. One pseudocount per base keeps every log-odds
            // finite.
            let at = (counts[0] + counts[3]) as f64 / 2.0;
            let cg = (counts[1] + counts[2]) as f64 / 2.0;
            let denom = 2.0 * (at + cg) + 4.0;
            let a = (at + 1.0) / denom;
            let c = (cg + 1.0) / denom;
            let values = [a, c, c, a];
            if values.iter().any(|value| !value.is_finite() || *value <= 0.0) {
                bail!("empirical background produced a non-finite base frequency");
            }
            Ok(BackgroundModel {
                values,
                mode: "empirical_promoter_strand_symmetric".to_string(),
                raw_counts: counts,
                ambiguous_bases: ambiguous,
                pseudocount: 1.0,
            })
        }
        other => bail!("--background-mode must be uniform or empirical, got {other}"),
    }
}

fn promoter_set_hash(promoters: &[promoter::Promoter]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(b"regutrace-promoter-set-v1\0");
    hasher.update((promoters.len() as u64).to_le_bytes());
    for promoter in promoters {
        hasher.update((promoter.gene_id.len() as u64).to_le_bytes());
        hasher.update(promoter.gene_id.as_bytes());
        hasher.update((promoter.seq.len() as u64).to_le_bytes());
        hasher.update(&promoter.seq);
    }
    format!("{:x}", hasher.finalize())
}

fn information_content(matrix: &[[f64; 4]]) -> f64 {
    matrix
        .iter()
        .map(|row| {
            let entropy: f64 = row
                .iter()
                .filter(|value| **value > 0.0)
                .map(|value| -value * value.log2())
                .sum();
            2.0 - entropy
        })
        .sum()
}

fn max_score(matrix: &[[f64; 4]], background: &[f64; 4]) -> f64 {
    matrix
        .iter()
        .map(|row| {
            (0..4)
                .map(|base| {
                    let probability = row[base].max(1e-9);
                    let background = background[base].max(1e-9);
                    (probability / background).log2()
                })
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .sum()
}

fn build_groups(
    args: &Args,
    motif_to_tfs: &BTreeMap<String, Vec<&io::TfMotif>>,
    motif_by_id: &HashMap<String, pwm::Pwm>,
    background: &[f64; 4],
    key_base: &ScanKeyBase,
) -> Vec<PwmGroup> {
    let mut groups = Vec::new();
    if args.reuse_exact_pwm {
        let mut by_digest: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for motif_id in motif_to_tfs.keys() {
            let digest = motif_by_id
                .get(motif_id)
                .expect("mapped motif exists")
                .digest();
            by_digest.entry(digest).or_default().push(motif_id.clone());
        }
        for (digest, mut motif_ids) in by_digest {
            motif_ids.sort();
            let representative = motif_ids[0].clone();
            let pwm = motif_by_id
                .get(&representative)
                .expect("motif id exists")
                .clone();
            let tf_ids: BTreeSet<&str> = motif_ids
                .iter()
                .flat_map(|id| {
                    motif_to_tfs
                        .get(id)
                        .into_iter()
                        .flatten()
                        .map(|tf| tf.tf_id.as_str())
                })
                .collect();
            let min_pvalue = scan::Scanner::new(&pwm, background, SCAN_SCALE).min_pvalue();
            groups.push(PwmGroup {
                scan_key: key_base.key(&digest),
                digest,
                motif_ids,
                pwm,
                min_pvalue,
                tf_count: tf_ids.len(),
            });
        }
        groups.sort_by(|a, b| {
            a.motif_ids[0]
                .cmp(&b.motif_ids[0])
                .then(a.digest.cmp(&b.digest))
        });
    } else {
        for motif_id in motif_to_tfs.keys() {
            let pwm = motif_by_id
                .get(motif_id)
                .expect("mapped motif exists")
                .clone();
            let digest = pwm.digest();
            let tf_ids: BTreeSet<&str> = motif_to_tfs
                .get(motif_id)
                .into_iter()
                .flatten()
                .map(|tf| tf.tf_id.as_str())
                .collect();
            let min_pvalue = scan::Scanner::new(&pwm, background, SCAN_SCALE).min_pvalue();
            groups.push(PwmGroup {
                scan_key: key_base.key(&digest),
                digest,
                motif_ids: vec![motif_id.clone()],
                pwm,
                min_pvalue,
                tf_count: tf_ids.len(),
            });
        }
    }
    groups
}

fn genomic_coordinates(
    promoter: &promoter::Promoter,
    hit_start: u64,
    hit_end: u64,
    hit_strand: char,
) -> Option<(u64, u64, char)> {
    if promoter.chrom.is_empty() || promoter.start == 0 || promoter.end == 0 {
        return None;
    }
    if hit_start == 0 || hit_end < hit_start {
        return None;
    }
    let (start, end, strand) = if promoter.strand == '-' {
        let start = promoter.end.checked_sub(hit_end)?.checked_add(1)?;
        let end = promoter.end.checked_sub(hit_start)?.checked_add(1)?;
        let strand = if hit_strand == '+' { '-' } else { '+' };
        (start, end, strand)
    } else {
        let start = promoter.start.checked_add(hit_start)?.checked_sub(1)?;
        let end = promoter.start.checked_add(hit_end)?.checked_sub(1)?;
        (start, end, hit_strand)
    };
    if end > promoter.end {
        return None;
    }
    Some((start, end, strand))
}

fn open_tsv(path: &Path) -> Result<csv::Writer<BufWriter<File>>> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("create output directory {}", parent.display()))?;
        }
    }
    let file = File::create(path).with_context(|| format!("create {}", path.display()))?;
    Ok(WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(BufWriter::new(file)))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run(
    args: &Args,
    promoters: &[promoter::Promoter],
    annotation: &HashMap<String, String>,
    motif_to_tfs: &BTreeMap<String, Vec<&io::TfMotif>>,
    motif_by_id: &HashMap<String, pwm::Pwm>,
    threads: usize,
    pool: Option<&rayon::ThreadPool>,
    progress: &mut ProgressWriter,
) -> Result<Phase3Stats> {
    let out_path = args
        .out
        .as_ref()
        .context("--out is required for the Phase 3 scan path")?;
    let started = Instant::now();
    let background = background_model(args, promoters)?;
    if background.values.iter().any(|value| !value.is_finite() || *value <= 0.0) {
        bail!("background values must be finite and positive");
    }
    let promoter_hash = promoter_set_hash(promoters);
    let key_base = ScanKeyBase::new(
        &background.values,
        args.pvalue,
        args.min_score,
        args.both_strands,
        &promoter_hash,
    );
    let groups = build_groups(
        args,
        motif_to_tfs,
        motif_by_id,
        &background.values,
        &key_base,
    );
    if groups.is_empty() {
        bail!("no PWM groups to scan");
    }
    let unique_pwms: BTreeSet<&str> = groups.iter().map(|group| group.digest.as_str()).collect();
    let duplicate_pwm_groups = groups
        .iter()
        .filter(|group| group.motif_ids.len() > 1)
        .count();
    let cache = match &args.scan_cache {
        Some(dir) => Some(ScanCache::open(dir)?),
        None => None,
    };

    let total_work = groups.len() as u64 * promoters.len() as u64;
    let mut done_work = 0_u64;
    progress.write(
        "scanning",
        "running",
        0,
        total_work,
        "scanning unique PWM groups",
        true,
    );

    let mut edge_stream = io::EdgeStream::create(out_path)?;
    let mut site_writer = match &args.sites_out {
        Some(path) => {
            let mut writer = open_tsv(path)?;
            writer.write_record([
                "pwm_hash",
                "motif_ids",
                "gene_id",
                "chrom",
                "promoter_start",
                "promoter_end",
                "promoter_strand",
                "hit_start",
                "hit_end",
                "hit_strand",
                "genomic_hit_start",
                "genomic_hit_end",
                "genomic_hit_strand",
                "score",
                "p_value",
                "tf_count",
                "scan_key",
            ])?;
            Some(writer)
        }
        None => None,
    };

    let mut promoter_stats: Vec<PromoterStats> =
        (0..promoters.len()).map(|_| PromoterStats::new()).collect();
    let mut motif_id_hits: HashMap<String, u64> =
        motif_to_tfs.keys().map(|id| (id.clone(), 0)).collect();
    let mut unique_site_hits = 0_u64;
    let mut motif_id_hit_total = 0_u64;
    let mut evidence_rows = 0_u64;
    let mut executed_promoter_scans = 0_u64;
    let mut skipped_unattainable_scans = 0_u64;
    let mut testable_motifs = 0_usize;
    let mut scanned_groups = 0_usize;

    for (group_index, group) in groups.iter().enumerate() {
        let testable = group.min_pvalue <= args.pvalue;
        if testable {
            testable_motifs += group.motif_ids.len();
            scanned_groups += 1;
        } else {
            skipped_unattainable_scans += promoters.len() as u64;
            done_work += promoters.len() as u64;
            progress.write(
                "scanning",
                "running",
                done_work,
                total_work,
                "skipping PWMs whose p-value floor exceeds the requested threshold",
                false,
            );
            continue;
        }

        let hits_by_promoter = if let Some(cache) = &cache {
            match cache.load(&group.scan_key, promoters.len()) {
                Some(hits) => hits,
                None => {
                    let hits =
                        scan_group(promoters, group, args, pool, &background.values)?;
                    cache.store(&group.scan_key, &hits)?;
                    hits
                }
            }
        } else {
            scan_group(promoters, group, args, pool, &background.values)?
        };
        executed_promoter_scans += promoters.len() as u64;

        for (promoter_index, hits) in hits_by_promoter.iter().enumerate() {
            if hits.is_empty() {
                continue;
            }
            let promoter = &promoters[promoter_index];
            let stats = &mut promoter_stats[promoter_index];
            stats.pwms_with_hits += 1;
            for hit in hits {
                unique_site_hits += 1;
                stats.sites += 1;
                if hit.p_value < stats.best_pvalue
                    || (hit.p_value == stats.best_pvalue && hit.score > stats.best_score)
                {
                    stats.best_pvalue = hit.p_value;
                    stats.best_score = hit.score;
                    stats.best_pwm_hash = group.digest.clone();
                    stats.best_motif_ids = group.motif_ids.join(";");
                }
                for motif_id in &group.motif_ids {
                    *motif_id_hits.entry(motif_id.clone()).or_insert(0) += 1;
                    motif_id_hit_total += 1;
                    let tfs = motif_to_tfs
                        .get(motif_id)
                        .expect("group motif is mapped");
                    if args.sites_only {
                        continue;
                    }
                    for tf in tfs {
                        evidence_rows += 1;
                        let target_symbol = annotation
                            .get(&promoter.gene_id)
                            .cloned()
                            .unwrap_or_else(|| promoter.symbol.clone());
                        let mut notes = format!("family={}", tf.family);
                        if !tf.reference_tf_id.is_empty()
                            || !tf.mapping_method.is_empty()
                            || !tf.domain_mode.is_empty()
                        {
                            notes.push_str(&format!(
                                ";motif_mapping={}",
                                serde_json::json!({
                                    "reference_tf_id": tf.reference_tf_id,
                                    "method": tf.mapping_method,
                                    "domain_mode": tf.domain_mode,
                                    "extraction_family": tf.extraction_family,
                                    "extraction_models": tf.extraction_models,
                                    "domain_count": tf.domain_count,
                                })
                            ));
                        }
                        edge_stream.write(&io::EdgeRecord {
                            tf_id: tf.tf_id.clone(),
                            tf_symbol: if tf.tf_symbol.is_empty() {
                                tf.tf_id.clone()
                            } else {
                                tf.tf_symbol.clone()
                            },
                            target_id: promoter.gene_id.clone(),
                            target_symbol,
                            evidence_type: "prediction_motif".into(),
                            evidence_level: "D".into(),
                            method: "jaspar_pwm".into(),
                            source: "promoter_scan".into(),
                            pmid: String::new(),
                            condition: String::new(),
                            promoter_start: (promoter.start > 0).then_some(promoter.start),
                            promoter_end: (promoter.end > 0).then_some(promoter.end),
                            motif_id: motif_id.clone(),
                            hit_start: Some(hit.start),
                            hit_end: Some(hit.end),
                            strand: hit.strand.to_string(),
                            score: Some(hit.score),
                            p_value: Some(hit.p_value),
                            notes,
                        })?;
                    }
                }
                if let Some(writer) = site_writer.as_mut() {
                    let genomic = genomic_coordinates(promoter, hit.start, hit.end, hit.strand);
                    writer.write_record([
                        group.digest.as_str(),
                        group.motif_ids.join(";").as_str(),
                        promoter.gene_id.as_str(),
                        promoter.chrom.as_str(),
                        &promoter.start.to_string(),
                        &promoter.end.to_string(),
                        &promoter.strand.to_string(),
                        &hit.start.to_string(),
                        &hit.end.to_string(),
                        &hit.strand.to_string(),
                        &genomic
                            .map(|value| value.0.to_string())
                            .unwrap_or_default(),
                        &genomic
                            .map(|value| value.1.to_string())
                            .unwrap_or_default(),
                        &genomic
                            .map(|value| value.2.to_string())
                            .unwrap_or_default(),
                        &format!("{:.6}", hit.score),
                        &format!("{:.17e}", hit.p_value),
                        &group.tf_count.to_string(),
                        group.scan_key.as_str(),
                    ])?;
                }
            }
        }

        done_work += promoters.len() as u64;
        if group_index % 5 == 0 || group_index + 1 == groups.len() {
            progress.write(
                "scanning",
                "running",
                done_work,
                total_work,
                "scanning unique PWM groups",
                false,
            );
        }
    }
    edge_stream.flush()?;
    if let Some(mut writer) = site_writer {
        writer.flush()?;
    }

    let mut testability_rows = Vec::new();
    let mut unattainable_motifs = 0_usize;
    let mut motifs_with_hits = 0_usize;
    let mut motifs_without_hits = 0_usize;
    let mut motif_ids: Vec<&String> = motif_by_id.keys().collect();
    motif_ids.sort();
    for motif_id in motif_ids {
        let motif = motif_by_id.get(motif_id).expect("motif exists");
        let scanner = scan::Scanner::new(motif, &background.values, SCAN_SCALE);
        let min_pvalue = scanner.min_pvalue();
        let mapped = motif_to_tfs.contains_key(motif_id);
        let hits = *motif_id_hits.get(motif_id).unwrap_or(&0);
        let status = if !mapped {
            "not_scanned"
        } else if min_pvalue > args.pvalue {
            unattainable_motifs += 1;
            "threshold_unattainable"
        } else if hits > 0 {
            motifs_with_hits += 1;
            "significant_hits"
        } else {
            motifs_without_hits += 1;
            "no_significant_hits"
        };
        let tf_count = motif_to_tfs
            .get(motif_id)
            .map(|tfs| {
                tfs.iter()
                    .map(|tf| tf.tf_id.as_str())
                    .collect::<BTreeSet<_>>()
                    .len()
            })
            .unwrap_or(0);
        testability_rows.push(serde_json::json!({
            "motif_id": motif_id,
            "pwm_hash": motif.digest(),
            "matrix_length": motif.matrix.len(),
            "information_content": information_content(&motif.matrix),
            "max_score": max_score(&motif.matrix, &background.values),
            "min_pvalue": min_pvalue,
            "requested_pvalue": args.pvalue,
            "attainable": min_pvalue <= args.pvalue,
            "status": status,
            "mapped": mapped,
            "motif_id_hit_count": hits,
            "tf_count": tf_count,
            "background_mode": background.mode.clone(),
        }));
    }

    if let Some(path) = &args.motif_testability_out {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("create output directory {}", parent.display()))?;
            }
        }
        let file = File::create(path)
            .with_context(|| format!("create motif testability table {}", path.display()))?;
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(BufWriter::new(file));
        writer.write_record([
            "motif_id",
            "pwm_hash",
            "matrix_length",
            "information_content",
            "max_score",
            "min_pvalue",
            "requested_pvalue",
            "attainable",
            "status",
            "mapped",
            "motif_id_hit_count",
            "tf_count",
            "background_mode",
        ])?;
        for row in &testability_rows {
            writer.write_record([
                row["motif_id"].as_str().unwrap_or(""),
                row["pwm_hash"].as_str().unwrap_or(""),
                &row["matrix_length"].to_string(),
                &format!("{:.6}", row["information_content"].as_f64().unwrap_or(0.0)),
                &format!("{:.6}", row["max_score"].as_f64().unwrap_or(0.0)),
                &format!("{:.17e}", row["min_pvalue"].as_f64().unwrap_or(0.0)),
                &format!("{:.17e}", row["requested_pvalue"].as_f64().unwrap_or(0.0)),
                if row["attainable"].as_bool().unwrap_or(false) {
                    "true"
                } else {
                    "false"
                },
                row["status"].as_str().unwrap_or(""),
                if row["mapped"].as_bool().unwrap_or(false) {
                    "true"
                } else {
                    "false"
                },
                &row["motif_id_hit_count"].to_string(),
                &row["tf_count"].to_string(),
                row["background_mode"].as_str().unwrap_or(""),
            ])?;
        }
        writer.flush()?;
    }

    if let Some(path) = &args.promoter_summary_out {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("create output directory {}", parent.display()))?;
            }
        }
        let file = File::create(path)
            .with_context(|| format!("create promoter summary {}", path.display()))?;
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(BufWriter::new(file));
        writer.write_record([
            "gene_id",
            "chrom",
            "promoter_start",
            "promoter_end",
            "promoter_strand",
            "promoter_length",
            "scanned_pwms",
            "testable_pwms",
            "pwms_with_hits",
            "sites",
            "best_pvalue",
            "best_score",
            "best_pwm_hash",
            "best_motif_ids",
        ])?;
        let testable_groups = groups
            .iter()
            .filter(|group| group.min_pvalue <= args.pvalue)
            .count();
        for (promoter, stats) in promoters.iter().zip(&promoter_stats) {
            writer.write_record([
                promoter.gene_id.as_str(),
                promoter.chrom.as_str(),
                &promoter.start.to_string(),
                &promoter.end.to_string(),
                &promoter.strand.to_string(),
                &promoter.seq.len().to_string(),
                &groups.len().to_string(),
                &testable_groups.to_string(),
                &stats.pwms_with_hits.to_string(),
                &stats.sites.to_string(),
                &if stats.sites == 0 {
                    "1".to_string()
                } else {
                    format!("{:.17e}", stats.best_pvalue)
                },
                &if stats.sites == 0 {
                    String::new()
                } else {
                    format!("{:.6}", stats.best_score)
                },
                stats.best_pwm_hash.as_str(),
                stats.best_motif_ids.as_str(),
            ])?;
        }
        writer.flush()?;
    }

    let cache_hits = cache.as_ref().map(ScanCache::hit_count).unwrap_or(0);
    let cache_misses = cache.as_ref().map(ScanCache::miss_count).unwrap_or(0);
    let elapsed = started.elapsed().as_secs_f64();
    let mapped_motif_ids = motif_to_tfs.len();
    let unique_pwms = unique_pwms.len();
    let testability_summary = serde_json::json!({
        "mapped_motif_ids": mapped_motif_ids,
        "unique_pwms": unique_pwms,
        "testable_motifs": testable_motifs,
        "unattainable_motifs": unattainable_motifs,
        "motifs_with_hits": motifs_with_hits,
        "motifs_without_hits": motifs_without_hits,
        "not_scanned_motifs": motif_by_id.len().saturating_sub(mapped_motif_ids),
    });

    if let Some(path) = &args.scan_manifest_out {
        let (min_len, max_len, mean_len) = if promoters.is_empty() {
            (0, 0, 0.0)
        } else {
            let min_len = promoters
                .iter()
                .map(|promoter| promoter.seq.len())
                .min()
                .unwrap_or(0);
            let max_len = promoters
                .iter()
                .map(|promoter| promoter.seq.len())
                .max()
                .unwrap_or(0);
            let mean_len = promoters
                .iter()
                .map(|promoter| promoter.seq.len() as f64)
                .sum::<f64>()
                / promoters.len() as f64;
            (min_len, max_len, mean_len)
        };
        let manifest = serde_json::json!({
            "schema": "regutrace-scan-manifest-v1",
            "tool": "regutrace",
            "created_unix_seconds": std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|value| value.as_secs())
                .unwrap_or(0),
            "scan_algorithm": {
                "name": SCAN_ALGORITHM,
                "scale": SCAN_SCALE,
                "ambiguous_bases": "skip_window",
                "strands": if args.both_strands { "both" } else { "forward" },
            },
            "inputs": {
                "motifs": args.motifs.as_ref().map(|path| path.display().to_string()),
                "genome": args.genome.as_ref().map(|path| path.display().to_string()),
                "gff": args.gff.as_ref().map(|path| path.display().to_string()),
                "promoters": args.promoters.as_ref().map(|path| path.display().to_string()),
                "tf_map": args.tf_map.as_ref().map(|path| path.display().to_string()),
            },
            "promoters": {
                "count": promoters.len(),
                "hash": promoter_hash,
                "min_length": min_len,
                "max_length": max_len,
                "mean_length": mean_len,
            },
            "background": {
                "mode": background.mode.clone(),
                "values_acgt": background.values,
                "raw_counts_acgt": background.raw_counts,
                "ambiguous_bases": background.ambiguous_bases,
                "pseudocount_per_base": background.pseudocount,
                "symmetrized": background.mode.starts_with("empirical"),
            },
            "thresholds": {
                "requested_pvalue": args.pvalue,
                "min_score": args.min_score,
            },
            "sites_only": args.sites_only,
            "motifs": {
                "mapped_motif_ids": mapped_motif_ids,
                "unique_pwms": unique_pwms,
                "pwm_groups": groups.len(),
                "duplicate_pwm_groups": duplicate_pwm_groups,
                "exact_pwm_reuse": args.reuse_exact_pwm,
            },
            "scan_jobs": {
                "group_promoter_pairs": total_work,
                "executed_group_promoter_pairs": executed_promoter_scans,
                "skipped_unattainable_group_promoter_pairs": skipped_unattainable_scans,
                "cache_hits": cache_hits,
                "cache_misses": cache_misses,
                "cache_version": CACHE_VERSION,
                "threads": threads,
            },
            "sites": {
                "unique_pwm_sites": unique_site_hits,
                "motif_id_site_hits": motif_id_hit_total,
                "evidence_rows": evidence_rows,
            },
            "testability": testability_summary,
            "runtime": {
                "seconds": elapsed,
                "groups_per_second": groups.len() as f64 / elapsed.max(1e-9),
            },
            "outputs": {
                "evidence": out_path.display().to_string(),
                "sites": args.sites_out.as_ref().map(|path| path.display().to_string()),
                "motif_testability": args.motif_testability_out.as_ref().map(|path| path.display().to_string()),
                "promoter_summary": args.promoter_summary_out.as_ref().map(|path| path.display().to_string()),
                "scan_manifest": path.display().to_string(),
            },
        });
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("create output directory {}", parent.display()))?;
            }
        }
        std::fs::write(path, serde_json::to_vec_pretty(&manifest)?)
            .with_context(|| format!("write scan manifest {}", path.display()))?;
    }

    progress.write(
        "scanning",
        "running",
        total_work,
        total_work,
        "unique PWM scan complete",
        true,
    );
    Ok(Phase3Stats {
        unique_site_hits,
        motif_id_hits: motif_id_hit_total,
        evidence_rows,
        pwm_groups: groups.len(),
        mapped_motif_ids,
        unique_pwms,
        duplicate_pwm_groups,
        cache_hits,
        cache_misses,
        testable_groups: scanned_groups as u64,
        executed_promoter_scans,
        skipped_unattainable_scans,
        testable_motifs,
        unattainable_motifs,
        motifs_with_hits,
        motifs_without_hits,
        background,
        sites_path: args.sites_out.clone(),
        testability_path: args.motif_testability_out.clone(),
        promoter_summary_path: args.promoter_summary_out.clone(),
        manifest_path: args.scan_manifest_out.clone(),
    })
}

fn scan_group(
    promoters: &[promoter::Promoter],
    group: &PwmGroup,
    args: &Args,
    pool: Option<&rayon::ThreadPool>,
    background: &[f64; 4],
) -> Result<Vec<Vec<scan::RawHit>>> {
    let scanner = scan::Scanner::new(&group.pwm, background, SCAN_SCALE);
    Ok(match pool {
        Some(pool) => pool.install(|| {
            promoters
                .par_iter()
                .map(|promoter| {
                    scanner.scan(
                        &promoter.seq,
                        args.pvalue,
                        args.min_score,
                        args.both_strands,
                    )
                })
                .collect()
        }),
        None => promoters
            .iter()
            .map(|promoter| {
                scanner.scan(
                    &promoter.seq,
                    args.pvalue,
                    args.min_score,
                    args.both_strands,
                )
            })
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_promoter(gene_id: &str, strand: char, seq: &[u8]) -> promoter::Promoter {
        promoter::Promoter {
            gene_id: gene_id.to_string(),
            symbol: gene_id.to_string(),
            chrom: "chr1".to_string(),
            strand,
            start: 101,
            end: 500,
            seq: seq.to_vec(),
        }
    }

    #[test]
    fn pwm_digest_is_matrix_content_only() {
        let a = pwm::Pwm {
            id: "A".into(),
            name: "first".into(),
            matrix: vec![[0.7, 0.1, 0.1, 0.1], [0.1, 0.7, 0.1, 0.1]],
        };
        let mut b = a.clone();
        b.id = "B".into();
        b.name = "second".into();
        assert_eq!(a.digest(), b.digest());
        b.matrix[0][0] = 0.700_000_1;
        assert_ne!(a.digest(), b.digest());
    }

    #[test]
    fn empirical_background_is_strand_symmetric() {
        let promoters = vec![test_promoter("g1", '+', b"AAAACCCCGGGGTTTT")];
        let mut args = Args::default();
        args.background_mode = "empirical".into();
        let model = background_model(&args, &promoters).unwrap();
        assert_eq!(model.values[0], model.values[3]);
        assert_eq!(model.values[1], model.values[2]);
        let sum: f64 = model.values.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert_eq!(model.raw_counts, [4, 4, 4, 4]);
    }

    #[test]
    fn genomic_coordinates_follow_promoter_orientation() {
        let plus = test_promoter("g1", '+', &vec![b'A'; 400]);
        assert_eq!(
            genomic_coordinates(&plus, 51, 56, '+'),
            Some((151, 156, '+'))
        );
        let minus = test_promoter("g2", '-', &vec![b'A'; 400]);
        assert_eq!(
            genomic_coordinates(&minus, 1, 10, '+'),
            Some((491, 500, '-'))
        );
        assert_eq!(
            genomic_coordinates(&minus, 1, 10, '-'),
            Some((491, 500, '+'))
        );
        let no_coords = promoter::Promoter {
            chrom: String::new(),
            start: 0,
            end: 0,
            ..plus
        };
        assert_eq!(genomic_coordinates(&no_coords, 1, 10, '+'), None);
    }

    #[test]
    fn scan_cache_roundtrip_preserves_hits() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace_phase3_cache_{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        let cache = ScanCache::open(&dir).unwrap();
        let hits = vec![
            vec![scan::RawHit {
                start: 3,
                end: 8,
                strand: '+',
                score: 9.5,
                p_value: 1e-5,
            }],
            Vec::new(),
        ];
        cache.store("abc123", &hits).unwrap();
        let loaded = cache.load("abc123", 2).unwrap();
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded[0][0].start, 3);
        assert_eq!(loaded[0][0].end, 8);
        assert_eq!(loaded[0][0].strand, '+');
        assert_eq!(loaded[1].len(), 0);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
