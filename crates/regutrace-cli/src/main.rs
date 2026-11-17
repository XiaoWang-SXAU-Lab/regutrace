//! regutrace: local TF–Target candidate builder (Phase 1).
//!
//! Example:
//!   regutrace scan --genome genome.fa --gff annotation.gff3 \
//!     --motifs jaspar_plants.txt --tf-map tf_map.tsv \
//!     --promoter-len 1000 --pvalue 1e-4 --out evidence.tsv

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::{bail, Context, Result};
use csv::{ReaderBuilder, StringRecord, WriterBuilder};
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use regutrace_core::{diagnostics, fasta, gff, io, promoter, pwm, scan, store::Store};
use rustyhmmer::api::{Cutoff, HmmAnnotator};
use rustyhmmer::seqio::Seq;
use rusqlite::types::Value;
use rusqlite::{params_from_iter, params, Connection, OpenFlags};
use serde_json::Value as JsonValue;

mod motif_transfer;
mod pipeline;
mod plant_tfs;
mod trace;

struct Args {
    genome: Option<PathBuf>,
    gff: Option<PathBuf>,
    motifs: Option<PathBuf>,
    tf_map: Option<PathBuf>,
    tf_genes: Option<PathBuf>,
    include_single_domain_fallback: bool,
    out: Option<PathBuf>,
    annotation: Option<PathBuf>,
    allow_empty_tf_map: bool,
    user_list: Option<PathBuf>,
    promoters: Option<PathBuf>,
    promoters_out: Option<PathBuf>,
    promoter_len: usize,
    downstream_len: usize,
    pvalue: f64,
    min_score: Option<f64>,
    background: [f64; 4],
    both_strands: bool,
    max_genes: Option<usize>,
    gene_biotype: Option<String>,
    gene_biotype_from_pack: bool,
    report_json: Option<PathBuf>,
    progress_json: Option<PathBuf>,
    stream: bool,
    threads: Option<usize>,
    strict: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            genome: None,
            gff: None,
            motifs: None,
            tf_map: None,
            tf_genes: None,
            include_single_domain_fallback: false,
            out: None,
            annotation: None,
            allow_empty_tf_map: false,
            user_list: None,
            promoters: None,
            promoters_out: None,
            promoter_len: 1000,
            downstream_len: 0,
            pvalue: 1e-4,
            min_score: None,
            background: [0.25; 4],
            both_strands: true,
            max_genes: None,
            gene_biotype: None,
            gene_biotype_from_pack: false,
            report_json: None,
            progress_json: None,
            stream: false,
            threads: None,
            strict: false,
        }
    }
}

fn usage() {
    eprintln!(
        "regutrace (Phase 1) — local TF–Target candidate builder\n\
\n\
USAGE:\n\
  regutrace scan --genome <fa[.gz]> --gff <gff3[.gz]> --motifs <jaspar.txt> \\\n\
    [--tf-map <tf_map.tsv>] [--tf-genes <tf_genes.tsv>] --out <edges.tsv> [options]\n\
  regutrace pipeline --reference-pack <dir> --genome <fa> --gff <gff3> \\\n\
    (--proteins <faa> | --cds <fna>) --out-dir <dir> [options]\n\
  regutrace trace --tf-genes <tsv> --tf-map <tsv> --evidence <tsv> --motifs <jaspar> \\\n\
    [--reference-map <tsv>] [--pack-version <label>] [--include-single-domain-fallback] --out-dir <dir>\n\
  regutrace import --db <evidence.db> --evidence <edges.tsv>\n\
  regutrace export --db <evidence.db> --out <edges.tsv>\n\
  regutrace stats --db <evidence.db>\n\
  regutrace convert --input <edges.tsv|csv[.gz]> --out <file> \\\n\
    --format sqlite|csv|tsv|jsonl [--gzip]\n\
  regutrace prior --db <tf_target.db> --targets <genes.txt> --out <upstream.tsv>\n\
  regutrace identify-plant-tfs (--proteins <faa[.gz]> | --cds <fna[.gz]>) (--hmm <models.hmm> | --hmm-list <files.txt>) --rules <rules.json>\n\
    (--gene-map <map.tsv> | --gff <gff3> | --gene-ids) [--hmmsearch <exe>] [--threads <n>] --out <tf_genes.tsv>\n\
  regutrace map-tf-motifs --tf-genes <tsv> --reference-map <tsv> (--forward <blast.tsv> --reverse <blast.tsv> | --proteins <faa> --reference-domains <faa> --dbd-models <json>) --out <tf_map.tsv>\n\
  regutrace classify-tfs --domtblout <hits> --rules <rules.json> --out <proteins.tsv>\n\
  regutrace identify-tfs --proteins <representative.faa> --hmm <tf.hmm> \\\n\
    [--family-map <model_family.tsv>] [--evalue <e>] [--min-coverage <0..1>] --out <tf_genes.tsv>\n\
  regutrace validate-db --db <tf_target.db>\n\
  regutrace validate [--genome <fa> --gff <gff3>] [--motifs <jaspar>] \\\n\
    [--tf-map <tsv>] [--annotation <tsv>] [--user-list <tsv>] \\\n\
    [--report-json <report.json>] [--strict]\n\
\n\
PIPELINE OPTIONS:\n\
  --reference-pack <dir>   versioned HMM/rules/reference-motif pack\n\
  --out-dir <dir>          write tf_genes, tf_map, evidence, trace tables and SQLite here\n\
  --proteins <faa> | --cds <fna>  user protein or CDS FASTA\n\
  --include-single-domain-fallback  also scan mappings marked fallback\n\
  --background <A,C,G,T>  base probabilities (default 0.25,0.25,0.25,0.25)\n\
  --hmmsearch <exe>       HMMER executable (default: manifest or PATH)\n\
  --blastp <exe>          BLASTP executable (default: manifest or PATH)\n\
\n\
TRACE OPTIONS:\n\
  --tf-genes <tsv>        identified TF catalogue\n\
  --tf-map <tsv>          TF–motif mapping used for the scan\n\
  --evidence <tsv>        scan evidence table\n\
  --motifs <jaspar>       motif library used for the scan\n\
  --reference-map <tsv>   optional reference TF map (adds reference symbols/families)\n\
  --pack-version <label>  optional reference-pack label recorded in the trace\n\
  --out-dir <dir>         write tf_motif_trace, motif_ambiguity, target_candidate_ties\n\
\n\
SCAN OPTIONS:\n\
  --promoter-len <bp>     upstream length (default 1000)\n\
  --downstream <bp>       downstream window incl. TSS (default 0; PlantRegMap uses 100)\n\
  --tf-map <tsv>          TF–motif map (required with --tf-genes)\n\
  --include-single-domain-fallback  also scan mappings marked single_domain_fallback\n\
  --tf-genes <tsv>        identified species TF catalogue (gene_id, family, status);\n\
                          enables mapped-only validation of tf-map entries\n\
\n\
IDENTIFY-TFS OPTIONS:\n\
  --proteins <faa>        one representative protein per gene\n\
  --hmm <hmm>             HMMER3 TF/domain library\n\
  --family-map <tsv>      hmm_model/family mapping (optional; defaults to model name)\n\
  --evalue <e>            domain E-value cutoff (default 1e-5)\n\
  --min-coverage <0..1>   minimum HMM coverage (default 0.5)\n\
  --promoters <fa>        scan a promoter FASTA directly (skip genome/GFF)\n\
  --promoters-out <fa>    export extracted promoters as FASTA\n\
  --pvalue <p>            PWM hit p-value cutoff (default 1e-4)\n\
  --min-score <s>         minimum PWM score (log2 odds; for benchmarks)\n\
  --background <A,C,G,T>  base probabilities (default 0.25,0.25,0.25,0.25)\n\
  --annotation <tsv>      gene_id/symbol annotation for target symbols\n\
  --user-list <tsv>       user TF–target evidence list to append\n\
  --max-genes <n>         only process the first n genes (testing)\n\
  --gene-biotype <value>  only scan genes with this GFF gene_biotype\n\
  --one-strand            scan only the plus strand of each promoter\n\
  --both-strands          scan both strands (default)\n\
  --report-json <path>    write a machine-readable diagnostics report\n\
  --progress-json <path>  write live stage/progress JSON for frontends\n\
  --stream                write evidence rows incrementally (unsorted; large scans)\n\
  --threads <n>           scan worker threads (default CPU cores - 1; 1 = serial)\n\
  --strict                treat warnings as failure (validate)\n\
  -h, --help              show this help\n\
\n\
tf-genes columns: gene_id, tf_symbol, family, identification_status, source\n\
tf-map columns: tf_id, tf_symbol, family, motif_id (one row per TF–motif pair);\n\
if omitted, tf_id/tf_symbol are taken from the motif header name\n\
user-list columns: tf, target, evidence_type, evidence_level, method, source, pmid, condition, score, notes"
    );
}

fn parse_background(s: &str) -> Result<[f64; 4]> {
    let vals: Vec<f64> = s
        .split(',')
        .map(|x| x.trim().parse::<f64>())
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("invalid background: {s}"))?;
    if vals.len() != 4 || vals.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        bail!("background must be four finite positive numbers A,C,G,T");
    }
    // Scale by the maximum first so four large finite weights cannot overflow
    // the sum and normalize to all zeros.
    let max = vals.iter().copied().fold(0.0_f64, f64::max);
    let scaled: Vec<f64> = vals.iter().map(|v| v / max).collect();
    let sum: f64 = scaled.iter().sum();
    if !sum.is_finite() || sum <= 0.0 {
        bail!("background sum is not finite");
    }
    let out = [
        scaled[0] / sum,
        scaled[1] / sum,
        scaled[2] / sum,
        scaled[3] / sum,
    ];
    if out.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        bail!("background values are too disparate after normalization");
    }
    Ok(out)
}

fn parse_scan_args(args: &[String]) -> Result<Args> {
    let mut a = Args::default();
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--genome" => a.genome = Some(PathBuf::from(it.next().context("--genome needs a value")?)),
            "--gff" => a.gff = Some(PathBuf::from(it.next().context("--gff needs a value")?)),
            "--motifs" => a.motifs = Some(PathBuf::from(it.next().context("--motifs needs a value")?)),
            "--tf-map" => a.tf_map = Some(PathBuf::from(it.next().context("--tf-map needs a value")?)),
            "--tf-genes" => {
                a.tf_genes = Some(PathBuf::from(it.next().context("--tf-genes needs a value")?))
            }
            "--out" => a.out = Some(PathBuf::from(it.next().context("--out needs a value")?)),
            "--allow-empty-tf-map" => a.allow_empty_tf_map = true,
            "--annotation" => {
                a.annotation = Some(PathBuf::from(it.next().context("--annotation needs a value")?))
            }
            "--user-list" => {
                a.user_list = Some(PathBuf::from(it.next().context("--user-list needs a value")?))
            }
            "--promoters" => {
                a.promoters = Some(PathBuf::from(it.next().context("--promoters needs a value")?))
            }
            "--promoters-out" => {
                a.promoters_out =
                    Some(PathBuf::from(it.next().context("--promoters-out needs a value")?))
            }
            "--min-score" => {
                a.min_score = Some(
                    it.next()
                        .context("--min-score needs a value")?
                        .parse()
                        .context("invalid --min-score")?,
                )
            }
            "--promoter-len" => {
                a.promoter_len = it
                    .next()
                    .context("--promoter-len needs a value")?
                    .parse()
                    .context("invalid --promoter-len")?
            }
            "--downstream" => {
                a.downstream_len = it
                    .next()
                    .context("--downstream needs a value")?
                    .parse()
                    .context("invalid --downstream")?
            }
            "--pvalue" => {
                a.pvalue = it
                    .next()
                    .context("--pvalue needs a value")?
                    .parse()
                    .context("invalid --pvalue")?
            }
            "--background" => {
                a.background = parse_background(&it.next().context("--background needs a value")?)?
            }
            "--max-genes" => {
                a.max_genes = Some(
                    it.next()
                        .context("--max-genes needs a value")?
                        .parse()
                        .context("invalid --max-genes")?,
                )
            }
            "--gene-biotype" => {
                let value = it.next().context("--gene-biotype needs a value")?;
                if value.trim().is_empty() {
                    bail!("--gene-biotype must be nonempty");
                }
                a.gene_biotype = Some(value.to_owned());
            }
            "--gene-biotype-from-pack" => a.gene_biotype_from_pack = true,
            "--report-json" => {
                a.report_json = Some(PathBuf::from(
                    it.next().context("--report-json needs a value")?,
                ))
            }
            "--progress-json" => {
                a.progress_json = Some(PathBuf::from(
                    it.next().context("--progress-json needs a value")?,
                ))
            }
            "--stream" => a.stream = true,
            "--threads" => {
                let n: usize = it
                    .next()
                    .context("--threads needs a value")?
                    .parse()
                    .context("invalid --threads")?;
                if n == 0 {
                    a.threads = None;
                } else {
                    a.threads = Some(n);
                }
            }
            "--include-single-domain-fallback" => a.include_single_domain_fallback = true,
            "--strict" => a.strict = true,
            "--one-strand" => a.both_strands = false,
            "--both-strands" => a.both_strands = true,
            "-h" | "--help" => {
                usage();
                std::process::exit(0);
            }
            other => bail!("unknown argument: {other}"),
        }
    }
    if a.promoter_len == 0 {
        bail!("--promoter-len must be > 0");
    }
    if !(0.0..=1.0).contains(&a.pvalue) {
        bail!("--pvalue must be in [0, 1]");
    }
    if a.threads.is_some_and(|n| n > 256) {
        bail!("--threads must be between 1 and 256");
    }
    if a.gene_biotype_from_pack && a.gene_biotype.is_none() {
        bail!("--gene-biotype-from-pack requires --gene-biotype");
    }
    Ok(a)
}

fn required<T>(v: Option<T>, name: &str) -> Result<T> {
    v.with_context(|| format!("missing required argument {name}"))
}

fn effective_threads(requested: Option<usize>) -> usize {
    if let Some(n) = requested {
        if n > 0 {
            return n;
        }
    }
    std::thread::available_parallelism()
        .map(|n| n.get().saturating_sub(1).max(1))
        .unwrap_or(1)
}

/// Writes a small progress document for the dashboard. Updates are throttled
/// so scanning hundreds of thousands of promoter-motif pairs does not turn
/// progress reporting into the bottleneck. The file is replaced atomically so
/// a polling reader never sees a partial document.
struct ProgressWriter {
    path: Option<PathBuf>,
    last_write: Instant,
    last_progress: f64,
}

impl ProgressWriter {
    fn new(path: Option<PathBuf>) -> Self {
        Self {
            path,
            // Allow the first update immediately.
            last_write: Instant::now() - Duration::from_secs(1),
            last_progress: -1.0,
        }
    }

    fn write(
        &mut self,
        stage: &str,
        status: &str,
        done: u64,
        total: u64,
        message: &str,
        force: bool,
    ) {
        let Some(path) = &self.path else {
            return;
        };
        let progress = if total == 0 {
            0.0
        } else {
            (done as f64 / total as f64).clamp(0.0, 1.0)
        };
        if !force
            && self.last_write.elapsed() < Duration::from_millis(250)
            && (progress - self.last_progress).abs() < 0.005
        {
            return;
        }
        let doc = serde_json::json!({
            "stage": stage,
            "status": status,
            "done": done,
            "total": total,
            "progress": progress,
            "message": message,
        });
        let Ok(text) = serde_json::to_string_pretty(&doc) else {
            return;
        };
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                let _ = std::fs::create_dir_all(parent);
            }
        }
        let tmp = path.with_extension("json.tmp");
        if std::fs::write(&tmp, &text).is_ok() {
            // rename over an existing file is not portable on Windows.
            let _ = std::fs::remove_file(path);
            if std::fs::rename(&tmp, path).is_err() {
                let _ = std::fs::write(path, &text);
            }
        }
        self.last_write = Instant::now();
        self.last_progress = progress;
    }
}

struct ChunkScan {
    edges: Vec<io::EdgeRecord>,
    hits: usize,
}

enum EdgeSink {
    Stream(Mutex<io::EdgeStream>),
    Memory(Mutex<Vec<io::EdgeRecord>>),
}

fn scan_chunk(
    promoters: &[promoter::Promoter],
    motif_ids: &[String],
    motif_to_tfs: &BTreeMap<String, Vec<&io::TfMotif>>,
    scanners: &HashMap<String, scan::Scanner>,
    annotation: &HashMap<String, String>,
    pvalue: f64,
    min_score: Option<f64>,
    both_strands: bool,
    done: &AtomicU64,
) -> ChunkScan {
    let mut edges = Vec::new();
    let mut hits_count = 0usize;
    for prom in promoters {
        for motif_id in motif_ids {
            let tfs = motif_to_tfs
                .get(motif_id)
                .expect("motif id comes from the same map");
            let scanner = scanners
                .get(motif_id)
                .expect("scanner exists for every motif");
            let hits = scanner.scan(&prom.seq, pvalue, min_score, both_strands);
            for hit in hits {
                hits_count += 1;
                for tf in tfs {
                    edges.push(io::EdgeRecord {
                        tf_id: tf.tf_id.clone(),
                        tf_symbol: if tf.tf_symbol.is_empty() {
                            tf.tf_id.clone()
                        } else {
                            tf.tf_symbol.clone()
                        },
                        target_id: prom.gene_id.clone(),
                        target_symbol: annotation
                            .get(&prom.gene_id)
                            .cloned()
                            .unwrap_or_else(|| prom.symbol.clone()),
                        evidence_type: "prediction_motif".into(),
                        evidence_level: "D".into(),
                        method: "jaspar_pwm".into(),
                        source: "promoter_scan".into(),
                        pmid: String::new(),
                        condition: String::new(),
                        promoter_start: (prom.start > 0).then_some(prom.start),
                        promoter_end: (prom.end > 0).then_some(prom.end),
                        motif_id: motif_id.clone(),
                        hit_start: Some(hit.start),
                        hit_end: Some(hit.end),
                        strand: hit.strand.to_string(),
                        score: Some(hit.score),
                        p_value: Some(hit.p_value),
                        notes: {
                            let mut notes = format!("family={}", tf.family);
                            if !tf.reference_tf_id.is_empty()
                                || !tf.mapping_method.is_empty()
                                || !tf.domain_mode.is_empty()
                            {
                                notes.push_str(&format!(";motif_mapping={}", serde_json::json!({
                                    "reference_tf_id": tf.reference_tf_id,
                                    "method": tf.mapping_method,
                                    "domain_mode": tf.domain_mode,
                                    "extraction_family": tf.extraction_family,
                                    "extraction_models": tf.extraction_models,
                                    "domain_count": tf.domain_count,
                                })));
                            }
                            notes
                        },
                    });
                }
            }
        }
        done.fetch_add(motif_ids.len() as u64, Ordering::Relaxed);
    }
    ChunkScan {
        edges,
        hits: hits_count,
    }
}

fn write_promoters_fasta(path: &std::path::Path, promoters: &[promoter::Promoter]) -> Result<()> {
    use std::io::Write;
    let file = std::fs::File::create(path)
        .with_context(|| format!("create promoter FASTA {}", path.display()))?;
    let mut w = std::io::BufWriter::new(file);
    for p in promoters {
        if p.start > 0 && !p.chrom.is_empty() {
            writeln!(
                w,
                ">{}|{}:{}-{}({})",
                p.gene_id, p.chrom, p.start, p.end, p.strand
            )?;
        } else {
            writeln!(w, ">{}", p.gene_id)?;
        }
        for chunk in p.seq.chunks(60) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()?;
    Ok(())
}

fn report_path_from_args(raw: &[String]) -> Option<PathBuf> {
    let mut it = raw.iter();
    while let Some(arg) = it.next() {
        if arg == "--report-json" {
            return it.next().map(PathBuf::from);
        }
    }
    None
}

fn progress_path_from_args(raw: &[String]) -> Option<PathBuf> {
    let mut it = raw.iter();
    while let Some(arg) = it.next() {
        if arg == "--progress-json" {
            return it.next().map(PathBuf::from);
        }
    }
    None
}

fn write_argument_error_report(
    command: &str,
    path: &Option<PathBuf>,
    err: &anyhow::Error,
) -> Result<()> {
    let mut report = diagnostics::Report::new(command);
    report.push(diagnostics::Diagnostic::error(
        if command == "scan" {
            "SCAN_FAILED"
        } else {
            "VALIDATE_FAILED"
        },
        format!("{err:#}"),
    ));
    for line in report.human_lines() {
        eprintln!("{line}");
    }
    if let Some(path) = path {
        report.write_json(path)?;
    }
    Ok(())
}

fn cmd_validate(raw: &[String]) -> Result<()> {
    let report_path = report_path_from_args(raw);
    let progress_path = progress_path_from_args(raw);
    let args = match parse_scan_args(raw) {
        Ok(args) => args,
        Err(e) => {
            write_argument_error_report("validate", &report_path, &e)?;
            ProgressWriter::new(progress_path).write(
                "failed",
                "error",
                0,
                0,
                &format!("{e:#}"),
                true,
            );
            return Err(e);
        }
    };
    let mut progress = ProgressWriter::new(args.progress_json.clone());
    progress.write(
        "validating",
        "running",
        0,
        1,
        "validating inputs",
        true,
    );
    let mut report = diagnostics::Report::new("validate");
    let result = run_validate(&args, &mut report, &mut progress);
    if let Err(e) = &result {
        report.push(diagnostics::Diagnostic::error(
            "VALIDATE_FAILED",
            format!("{e:#}"),
        ));
    }
    if args.strict && report.summary.warnings > 0 && !report.has_errors() {
        report.error(
            "STRICT_WARNINGS",
            "validation failed because --strict is set",
        );
    }
    for line in report.human_lines() {
        eprintln!("{line}");
    }
    if let Some(path) = &args.report_json {
        report.write_json(path)?;
    }
    if report.has_errors() {
        progress.write("failed", "error", 0, 0, "input validation failed", true);
        return Err(anyhow::anyhow!("input validation failed"));
    }
    progress.write("done", "done", 1, 1, "validation complete", true);
    Ok(())
}

fn filter_gene_biotype(
    mut genes: Vec<gff::GeneModel>,
    wanted: Option<&str>,
) -> Vec<gff::GeneModel> {
    if let Some(wanted) = wanted {
        genes.retain(|gene| gene.biotype.as_deref() == Some(wanted));
    }
    genes
}

fn apply_gene_biotype_filter(
    genes: Vec<gff::GeneModel>,
    wanted: Option<&str>,
    from_pack_default: bool,
    report: &mut diagnostics::Report,
) -> Result<Vec<gff::GeneModel>> {
    let Some(wanted) = wanted else {
        return Ok(genes);
    };
    let total = genes.len();
    let with_field = genes.iter().filter(|gene| gene.biotype.is_some()).count();
    let matched = genes
        .iter()
        .filter(|gene| gene.biotype.as_deref() == Some(wanted))
        .count();
    if with_field == 0 {
        let message = format!(
            "GFF has no gene_biotype attribute, so the gene_biotype={wanted} filter cannot be applied"
        );
        if from_pack_default {
            report.push(
                diagnostics::Diagnostic::warning(
                    "GENE_BIOTYPE_FILTER_SKIPPED",
                    format!(
                        "{message}; pack default ignored and all {total} gene(s) will be scanned"
                    ),
                )
                .with_suggestion(
                    "add gene_biotype to the GFF, or pass --gene-biotype explicitly to enforce a target universe",
                ),
            );
            return Ok(genes);
        }
        bail!(
            "{message}; provide a GFF with gene_biotype or remove --gene-biotype"
        );
    }
    if matched == 0 {
        let available: Vec<_> = genes
            .iter()
            .filter_map(|gene| gene.biotype.as_deref())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if from_pack_default {
            report.push(
                diagnostics::Diagnostic::warning(
                    "GENE_BIOTYPE_FILTER_SKIPPED",
                    format!(
                        "no gene has gene_biotype={wanted}; pack default ignored and all {total} gene(s) will be scanned"
                    ),
                )
                .with_context("available", available.join(","))
                .with_suggestion(
                    "check the GFF gene_biotype values or pass --gene-biotype explicitly",
                ),
            );
            return Ok(genes);
        }
        bail!(
            "no gene has gene_biotype={wanted}; available gene_biotype values: {}",
            available.join(",")
        );
    }
    report.info(
        "GENE_BIOTYPE_FILTER",
        format!("kept {matched} of {total} gene(s) with gene_biotype={wanted}"),
    );
    Ok(filter_gene_biotype(genes, Some(wanted)))
}

#[derive(Clone, Copy, Debug, Default)]
struct MappingFilterStats {
    total: usize,
    fallback: usize,
    retained: usize,
    excluded: usize,
}

fn is_single_domain_fallback(tf: &io::TfMotif) -> bool {
    tf.domain_mode.trim().eq_ignore_ascii_case("single_domain_fallback")
}

fn filter_scan_mappings(
    mut tf_map: Vec<io::TfMotif>,
    include_fallback: bool,
) -> (Vec<io::TfMotif>, MappingFilterStats) {
    let total = tf_map.len();
    let fallback = tf_map.iter().filter(|tf| is_single_domain_fallback(tf)).count();
    if include_fallback {
        return (
            tf_map,
            MappingFilterStats {
                total,
                fallback,
                retained: total,
                excluded: 0,
            },
        );
    }
    tf_map.retain(|tf| !is_single_domain_fallback(tf));
    let retained = tf_map.len();
    (
        tf_map,
        MappingFilterStats {
            total,
            fallback,
            retained,
            excluded: total - retained,
        },
    )
}

fn report_mapping_filter(
    report: &mut diagnostics::Report,
    stats: MappingFilterStats,
    include_fallback: bool,
) {
    report.info(
        "TFMAP_DOMAIN_FILTER",
        format!(
            "tf-map rows total={}, single_domain_fallback={}, retained={}, excluded={}{}",
            stats.total,
            stats.fallback,
            stats.retained,
            stats.excluded,
            if include_fallback {
                " (fallback inclusion requested)"
            } else {
                " (use --include-single-domain-fallback to retain fallback rows)"
            }
        ),
    );
}

fn validate_tf_identity(
    tf_map: &[io::TfMotif],
    tf_genes: Option<&Path>,
    report: &mut diagnostics::Report,
) -> Result<()> {
    let Some(path) = tf_genes else {
        return Ok(());
    };
    let genes = if tf_map.is_empty() { io::read_tf_genes_allow_empty(path)? } else { io::read_tf_genes(path)? };
    let mut unknown = Vec::new();
    let mut not_identified = Vec::new();
    for tf in tf_map {
        match genes.get(&tf.tf_id) {
            None => unknown.push(tf.tf_id.clone()),
            Some(gene) if !gene.is_identified() => {
                not_identified.push(format!("{} ({})", tf.tf_id, gene.status))
            }
            Some(_) => {}
        }
    }
    unknown.sort();
    unknown.dedup();
    not_identified.sort();
    not_identified.dedup();
    if !unknown.is_empty() {
        bail!(
            "tf-map references {} TF gene(s) absent from --tf-genes: {}{}",
            unknown.len(),
            unknown.iter().take(10).cloned().collect::<Vec<_>>().join(", "),
            if unknown.len() > 10 { ", ..." } else { "" }
        );
    }
    if !not_identified.is_empty() {
        bail!(
            "tf-map references TF gene(s) not marked identified: {}{}",
            not_identified.iter().take(10).cloned().collect::<Vec<_>>().join(", "),
            if not_identified.len() > 10 { ", ..." } else { "" }
        );
    }
    report.add_input("tf_genes", path);
    report.info(
        "TF_IDENTITY_SUMMARY",
        format!(
            "mapped {} TF–motif row(s) to {} identified TF gene(s)",
            tf_map.len(),
            genes.values().filter(|gene| gene.is_identified()).count()
        ),
    );
    Ok(())
}

fn run_validate(
    args: &Args,
    report: &mut diagnostics::Report,
    progress: &mut ProgressWriter,
) -> Result<()> {
    let mut any = false;

    let motifs = if let Some(path) = &args.motifs {
        any = true;
        report.add_input("motifs", path);
        progress.write(
            "reading_motifs",
            "running",
            0,
            1,
            "reading motif library",
            true,
        );
        let motifs = pwm::read_jaspar(path)?;
        diagnostics::review_motifs(&motifs, path, report);
        Some(motifs)
    } else {
        None
    };

    if let Some(path) = &args.tf_map {
        any = true;
        report.add_input("tf_map", path);
        let tf_map = io::read_tf_map(path)?;
        let (tf_map, filter) =
            filter_scan_mappings(tf_map, args.include_single_domain_fallback);
        report_mapping_filter(report, filter, args.include_single_domain_fallback);
        if filter.retained == 0 {
            report.push(
                diagnostics::Diagnostic::error(
                    "TFMAP_ALL_FALLBACK_FILTERED",
                    "all tf-map rows are single_domain_fallback; validation fails because scan would have no eligible mappings (use --include-single-domain-fallback to retain them)",
                )
                .with_file(path),
            );
        } else {
            validate_tf_identity(&tf_map, args.tf_genes.as_deref(), report)?;
            if let Some(motifs) = &motifs {
                diagnostics::review_tf_map(&tf_map, motifs, path, report);
            } else {
                report.push(
                    diagnostics::Diagnostic::warning(
                        "TFMAP_WITHOUT_MOTIFS",
                        "tf-map provided without --motifs; motif matching was not checked",
                    )
                    .with_file(path),
                );
            }
        }
    } else if args.tf_genes.is_some() {
        bail!("--tf-genes requires --tf-map (motif-header auto mapping has no gene identities)");
    }

    let annotation = if let Some(path) = &args.annotation {
        any = true;
        report.add_input("annotation", path);
        Some(io::read_annotation(path)?)
    } else {
        None
    };

    if let Some(path) = &args.user_list {
        any = true;
        report.add_input("user_list", path);
        let rows = io::read_user_list(path)?;
        diagnostics::review_user_list(&rows, path, report);
    }

    if let Some(path) = &args.promoters {
        any = true;
        report.add_input("promoters", path);
        let fa = fasta::read_fasta(path)?;
        diagnostics::review_fasta(&fa, path, report);
        let mut promoters = Vec::new();
        for id in &fa.order {
            promoters.push(promoter::Promoter {
                gene_id: id.clone(),
                symbol: id.clone(),
                chrom: String::new(),
                strand: '+',
                start: 0,
                end: 0,
                seq: fa.seqs.get(id).unwrap().clone(),
            });
        }
        diagnostics::review_promoters(&promoters, None, report);
        if annotation.is_some() {
            report.info(
                "ANNOTATION_NOT_CHECKED",
                "annotation coverage needs gene models; provide --genome/--gff",
            );
        }
    } else if args.genome.is_some() || args.gff.is_some() {
        any = true;
        let genome_path = args
            .genome
            .as_ref()
            .context("--gff requires --genome")?;
        let gff_path = args.gff.as_ref().context("--genome requires --gff")?;
        report.add_input("genome", genome_path);
        report.add_input("gff", gff_path);
        progress.write("reading_genome", "running", 0, 1, "reading genome FASTA", true);
        let fa = fasta::read_fasta(genome_path)?;
        diagnostics::review_fasta(&fa, genome_path, report);
        progress.write("reading_gff", "running", 0, 1, "reading GFF3", true);
        let (genes, gff_stats) = gff::read_gff_with_stats(gff_path)?;
        diagnostics::review_gff(&genes, &gff_stats, gff_path, report);
        diagnostics::review_chromosomes(&genes, &fa, gff_path, genome_path, report);
        let genes = apply_gene_biotype_filter(
            genes,
            args.gene_biotype.as_deref(),
            args.gene_biotype_from_pack,
            report,
        )?;
        let mut promoters = Vec::new();
        let mut skipped = 0usize;
        for gene in &genes {
            match fa.seq(&gene.chrom) {
                Some(seq) => {
                    if let Some(p) =
                        promoter::extract_promoter(seq, gene, args.promoter_len, args.downstream_len)
                    {
                        promoters.push(p);
                    } else {
                        skipped += 1;
                    }
                }
                None => skipped += 1,
            }
        }
        if skipped > 0 {
            report.push(
                diagnostics::Diagnostic::warning(
                    "PROMOTER_SKIPPED",
                    format!("{skipped} gene(s) produced no promoter"),
                )
                .with_file(gff_path)
                .with_context("count", skipped)
                .with_suggestion("check chromosome names and TSS coordinates"),
            );
        }
        diagnostics::review_promoters(
            &promoters,
            Some(args.promoter_len + args.downstream_len),
            report,
        );
        if let (Some(annotation), Some(path)) = (&annotation, &args.annotation) {
            diagnostics::review_annotation(annotation, &genes, path, report);
        }
    }

    if !any {
        report.push(
            diagnostics::Diagnostic::error("NO_INPUTS", "no input files were provided")
                .with_suggestion(
                    "provide --genome/--gff, --promoters, --motifs, --tf-map, --annotation or --user-list",
                ),
        );
    }
    if !report.has_errors() {
        report.info("VALIDATE_OK", "input review completed without errors");
    }
    Ok(())
}

fn cmd_scan(raw: &[String]) -> Result<()> {
    let report_path = report_path_from_args(raw);
    let progress_path = progress_path_from_args(raw);
    let args = match parse_scan_args(raw) {
        Ok(args) => args,
        Err(e) => {
            write_argument_error_report("scan", &report_path, &e)?;
            ProgressWriter::new(progress_path).write(
                "failed",
                "error",
                0,
                0,
                &format!("{e:#}"),
                true,
            );
            return Err(e);
        }
    };
    let mut progress = ProgressWriter::new(args.progress_json.clone());
    progress.write("starting", "running", 0, 1, "starting scan", true);
    let mut report = diagnostics::Report::new("scan");
    let result = run_scan(&args, &mut report, &mut progress);
    match &result {
        Ok(()) => progress.write("done", "done", 1, 1, "scan complete", true),
        Err(e) => {
            report.push(diagnostics::Diagnostic::error(
                "SCAN_FAILED",
                format!("{e:#}"),
            ));
            progress.write("failed", "error", 0, 0, &format!("{e:#}"), true);
        }
    }
    for line in report.human_lines() {
        eprintln!("{line}");
    }
    if let Some(path) = &args.report_json {
        report.write_json(path)?;
    }
    result
}

fn run_scan(
    args: &Args,
    report: &mut diagnostics::Report,
    progress: &mut ProgressWriter,
) -> Result<()> {
    let motifs_path = required(args.motifs.clone(), "--motifs")?;
    let out_path = required(args.out.clone(), "--out")?;
    report.add_input("motifs", &motifs_path);

    progress.write("reading_motifs", "running", 0, 1, "reading motif library", true);
    eprintln!("reading motifs: {}", motifs_path.display());
    let motifs = pwm::read_jaspar(&motifs_path)?;
    diagnostics::review_motifs(&motifs, &motifs_path, report);
    if args.tf_genes.is_some() && args.tf_map.is_none() {
        bail!("--tf-genes requires --tf-map (motif-header auto mapping has no gene identities)");
    }
    let mut tf_map: Vec<io::TfMotif> = match &args.tf_map {
        Some(path) if args.allow_empty_tf_map => io::read_tf_map_allow_empty(path)?,
        Some(path) => io::read_tf_map(path)?,
        None => {
            eprintln!("no --tf-map given; deriving TF IDs from motif names");
            motifs
                .iter()
                .map(|m| {
                    let symbol = if m.name.trim().is_empty() {
                        m.id.clone()
                    } else {
                        m.name.trim().to_string()
                    };
                    io::TfMotif {
                        tf_id: symbol.clone(),
                        tf_symbol: symbol,
                        family: String::new(),
                        motif_id: m.id.clone(),
                        reference_tf_id: String::new(),
                        mapping_method: String::new(),
                        ..Default::default()
                    }
                })
                .collect()
        }
    };
    if let Some(path) = &args.tf_map {
        let (filtered, filter) =
            filter_scan_mappings(tf_map, args.include_single_domain_fallback);
        tf_map = filtered;
        report_mapping_filter(report, filter, args.include_single_domain_fallback);
        if filter.retained == 0 && !(args.allow_empty_tf_map && filter.total == 0) {
            report.push(
                diagnostics::Diagnostic::error(
                    "TFMAP_ALL_FALLBACK_FILTERED",
                    "all tf-map rows are single_domain_fallback; rerun with --include-single-domain-fallback to scan them",
                )
                .with_file(path),
            );
            bail!(
                "all {} tf-map row(s) are single_domain_fallback; rerun with --include-single-domain-fallback to scan them",
                filter.total
            );
        }
        report.add_input("tf_map", path);
        diagnostics::review_tf_map(&tf_map, &motifs, path, report);
        validate_tf_identity(&tf_map, args.tf_genes.as_deref(), report)?;
    } else {
        report.info(
            "TFMAP_AUTO",
            format!("derived {} TF–motif row(s) from motif names", tf_map.len()),
        );
    }
    let mut motif_by_id: HashMap<String, pwm::Pwm> = HashMap::new();
    for m in motifs {
        motif_by_id.entry(m.id.clone()).or_insert(m);
    }
    let annotation = match &args.annotation {
        Some(p) => {
            report.add_input("annotation", p);
            io::read_annotation(p)?
        }
        None => HashMap::new(),
    };

    let (promoters, gene_count, missing_chrom) = if let Some(promoters_path) = &args.promoters {
        progress.write(
            "reading_promoters",
            "running",
            0,
            1,
            "reading promoter FASTA",
            true,
        );
        eprintln!("reading promoter FASTA: {}", promoters_path.display());
        report.add_input("promoters", promoters_path);
        let fa = fasta::read_fasta(promoters_path)?;
        diagnostics::review_fasta(&fa, promoters_path, report);
        let mut v = Vec::new();
        for id in &fa.order {
            if args.max_genes.is_some_and(|n| v.len() >= n) {
                break;
            }
            v.push(promoter::Promoter {
                gene_id: id.clone(),
                symbol: id.clone(),
                chrom: String::new(),
                strand: '+',
                start: 0,
                end: 0,
                seq: fa.seqs.get(id).unwrap().clone(),
            });
        }
        diagnostics::review_promoters(&v, None, report);
        (v, fa.order.len(), 0usize)
    } else {
        let genome_path = required(args.genome.clone(), "--genome")?;
        let gff_path = required(args.gff.clone(), "--gff")?;
        report.add_input("genome", &genome_path);
        report.add_input("gff", &gff_path);
        progress.write(
            "reading_genome",
            "running",
            0,
            1,
            "reading genome FASTA",
            true,
        );
        eprintln!("reading genome: {}", genome_path.display());
        let genome = fasta::read_fasta(&genome_path)?;
        diagnostics::review_fasta(&genome, &genome_path, report);
        progress.write("reading_gff", "running", 0, 1, "reading GFF3", true);
        eprintln!("reading GFF: {}", gff_path.display());
        let (genes, gff_stats) = gff::read_gff_with_stats(&gff_path)?;
        diagnostics::review_gff(&genes, &gff_stats, &gff_path, report);
        diagnostics::review_chromosomes(&genes, &genome, &gff_path, &genome_path, report);
        if let Some(path) = &args.annotation {
            diagnostics::review_annotation(&annotation, &genes, path, report);
        }
        let genes = apply_gene_biotype_filter(
            genes,
            args.gene_biotype.as_deref(),
            args.gene_biotype_from_pack,
            report,
        )?;
        progress.write(
            "extracting_promoters",
            "running",
            0,
            genes.len() as u64,
            "extracting promoters",
            true,
        );
        let mut v = Vec::new();
        let mut missing = 0usize;
        for (i, gene) in genes.iter().enumerate() {
            if args.max_genes.is_some_and(|n| v.len() >= n) {
                break;
            }
            match genome.seq(&gene.chrom) {
                Some(seq) => {
                    if let Some(p) =
                        promoter::extract_promoter(seq, gene, args.promoter_len, args.downstream_len)
                    {
                        v.push(p);
                    }
                }
                None => missing += 1,
            }
            if i % 1000 == 0 {
                progress.write(
                    "extracting_promoters",
                    "running",
                    i as u64,
                    genes.len() as u64,
                    "extracting promoters",
                    false,
                );
            }
        }
        diagnostics::review_promoters(
            &v,
            Some(args.promoter_len + args.downstream_len),
            report,
        );
        (v, genes.len(), missing)
    };
    eprintln!(
        "genes: {}, promoters: {}, missing chromosomes: {}",
        gene_count,
        promoters.len(),
        missing_chrom
    );
    if promoters.is_empty() {
        bail!("no promoters available; check inputs and chromosome naming");
    }
    if let Some(path) = &args.promoters_out {
        progress.write(
            "writing_promoters",
            "running",
            0,
            1,
            "writing promoter FASTA",
            true,
        );
        write_promoters_fasta(path, &promoters)?;
        eprintln!("wrote promoter FASTA: {}", path.display());
    }

    // Motif -> TF mappings: scan each motif once per promoter, then emit all TFs
    // that use it.
    let mut motif_to_tfs: BTreeMap<String, Vec<&io::TfMotif>> = BTreeMap::new();
    let mut missing_motif = 0usize;
    for tf in &tf_map {
        if !motif_by_id.contains_key(&tf.motif_id) {
            missing_motif += 1;
            continue;
        }
        motif_to_tfs.entry(tf.motif_id.clone()).or_default().push(tf);
    }
    if motif_to_tfs.is_empty() && !(args.allow_empty_tf_map && tf_map.is_empty()) {
        bail!("no TF–motif rows matched the motif library");
    }
    if missing_motif > 0 {
        eprintln!("warning: {missing_motif} TF–motif rows reference motifs not in the library");
    }

    let mut scanners: HashMap<String, scan::Scanner> = HashMap::new();
    for motif_id in motif_to_tfs.keys() {
        let p = motif_by_id.get(motif_id).unwrap();
        scanners.insert(
            motif_id.clone(),
            scan::Scanner::new(p, &args.background, 10.0),
        );
    }

    let threads = effective_threads(args.threads);
    let pool = if threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .context("create scan thread pool")?,
        )
    } else {
        None
    };
    if args.stream {
        if let Some(parent) = out_path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("create output directory {}", parent.display()))?;
            }
        }
    }
    let sink = if args.stream {
        Arc::new(EdgeSink::Stream(Mutex::new(io::EdgeStream::create(
            &out_path,
        )?)))
    } else {
        Arc::new(EdgeSink::Memory(Mutex::new(Vec::new())))
    };
    let motif_ids: Vec<String> = motif_to_tfs.keys().cloned().collect();
    let total_work = promoters.len() as u64 * motif_ids.len() as u64;
    let chunk_size = if threads <= 1 {
        64
    } else {
        (promoters.len() / (threads * 8)).clamp(1, 64)
    };
    let chunks: Vec<&[promoter::Promoter]> = promoters.chunks(chunk_size).collect();
    let done = Arc::new(AtomicU64::new(0));
    let hits = Arc::new(AtomicU64::new(0));
    let edge_count = Arc::new(AtomicU64::new(0));
    progress.write(
        "scanning",
        "running",
        0,
        total_work,
        "scanning promoters",
        true,
    );

    let stop = Arc::new(AtomicBool::new(false));
    let reporter = if args.progress_json.is_some() {
        let done = done.clone();
        let stop = stop.clone();
        let path = args.progress_json.clone();
        Some(std::thread::spawn(move || {
            let mut p = ProgressWriter::new(path);
            while !stop.load(Ordering::Relaxed) {
                p.write(
                    "scanning",
                    "running",
                    done.load(Ordering::Relaxed),
                    total_work,
                    "scanning promoters",
                    false,
                );
                std::thread::sleep(Duration::from_millis(250));
            }
        }))
    } else {
        None
    };

    let handle_chunk = |chunk: &[promoter::Promoter]| -> Result<()> {
        let scanned = scan_chunk(
            chunk,
            &motif_ids,
            &motif_to_tfs,
            &scanners,
            &annotation,
            args.pvalue,
            args.min_score,
            args.both_strands,
            &done,
        );
        hits.fetch_add(scanned.hits as u64, Ordering::Relaxed);
        edge_count.fetch_add(scanned.edges.len() as u64, Ordering::Relaxed);
        match sink.as_ref() {
            EdgeSink::Stream(stream) => {
                let mut w = stream
                    .lock()
                    .map_err(|_| anyhow::anyhow!("evidence writer lock poisoned"))?;
                for e in &scanned.edges {
                    w.write(e)?;
                }
            }
            EdgeSink::Memory(memory) => {
                let mut v = memory
                    .lock()
                    .map_err(|_| anyhow::anyhow!("evidence buffer lock poisoned"))?;
                v.extend(scanned.edges);
            }
        }
        Ok(())
    };

    let scan_result = if let Some(pool) = &pool {
        pool.install(|| chunks.par_iter().try_for_each(|chunk| handle_chunk(chunk)))
    } else {
        chunks.iter().try_for_each(|chunk| handle_chunk(chunk))
    };
    stop.store(true, Ordering::Relaxed);
    if let Some(h) = reporter {
        let _ = h.join();
    }
    if scan_result.is_ok() {
        progress.write(
            "scanning",
            "running",
            total_work,
            total_work,
            "scan complete",
            true,
        );
    }
    scan_result?;
    let prediction_hits = hits.load(Ordering::Relaxed) as usize;
    let mut edge_count = edge_count.load(Ordering::Relaxed) as usize;

    let user_count = match &args.user_list {
        Some(path) => {
            report.add_input("user_list", path);
            let mut user = io::read_user_list(path)?;
            diagnostics::review_user_list(&user, path, report);
            let tf_symbol_by_id: HashMap<&str, &str> = tf_map
                .iter()
                .map(|t| (t.tf_id.as_str(), t.tf_symbol.as_str()))
                .collect();
            for e in &mut user {
                if e.tf_symbol.is_empty() {
                    if let Some(s) = tf_symbol_by_id.get(e.tf_id.as_str()) {
                        e.tf_symbol = (*s).to_string();
                    }
                }
                if e.target_symbol.is_empty() {
                    if let Some(s) = annotation.get(&e.target_id) {
                        e.target_symbol = s.clone();
                    }
                }
            }
            let n = user.len();
            edge_count += n;
            match sink.as_ref() {
                EdgeSink::Stream(stream) => {
                    let mut w = stream
                        .lock()
                        .map_err(|_| anyhow::anyhow!("evidence writer lock poisoned"))?;
                    for e in &user {
                        w.write(e)?;
                    }
                }
                EdgeSink::Memory(memory) => {
                    let mut v = memory
                        .lock()
                        .map_err(|_| anyhow::anyhow!("evidence buffer lock poisoned"))?;
                    v.extend(user);
                }
            }
            n
        }
        None => 0,
    };

    progress.write(
        "writing_evidence",
        "running",
        0,
        1,
        "writing evidence table",
        true,
    );
    match sink.as_ref() {
        EdgeSink::Stream(stream) => {
            let mut w = stream
                .lock()
                .map_err(|_| anyhow::anyhow!("evidence writer lock poisoned"))?;
            w.flush()?;
        }
        EdgeSink::Memory(memory) => {
            let mut edges = {
                let mut v = memory
                    .lock()
                    .map_err(|_| anyhow::anyhow!("evidence buffer lock poisoned"))?;
                std::mem::take(&mut *v)
            };
            if let Some(parent) = out_path.parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent)
                        .with_context(|| format!("create output directory {}", parent.display()))?;
                }
            }
            io::write_edges(&out_path, &mut edges)?;
        }
    }

    eprintln!(
        "motifs: {}, TF–motif pairs: {}, promoters: {}, prediction hits: {}, user edges: {}, total edges: {}, threads: {}",
        motif_to_tfs.len(),
        tf_map.len() - missing_motif,
        promoters.len(),
        prediction_hits,
        user_count,
        edge_count,
        threads
    );
    eprintln!("wrote {}", out_path.display());
    report.push(
        diagnostics::Diagnostic::info(
            "SCAN_SUMMARY",
            format!(
                "{} genes, {} promoters, {} motif hits, {} user rows, {} total evidence rows",
                gene_count,
                promoters.len(),
                prediction_hits,
                user_count,
                edge_count
            ),
        )
        .with_context("output", out_path.display().to_string())
        .with_context("genes", gene_count)
        .with_context("promoters", promoters.len())
        .with_context("motif_hits", prediction_hits)
        .with_context("user_rows", user_count)
        .with_context("total_edges", edge_count)
        .with_context("motifs_scanned", motif_to_tfs.len())
        .with_context("tf_motif_pairs", tf_map.len() - missing_motif)
        .with_context("threads", threads),
    );
    Ok(())
}

fn parse_options(args: &[String]) -> Result<HashMap<String, String>> {
    parse_options_with_flags(args, &[])
}

fn parse_options_with_flags(
    args: &[String],
    flags: &[&str],
) -> Result<HashMap<String, String>> {
    let mut out = HashMap::new();
    let mut i = 0;
    while i < args.len() {
        let key = &args[i];
        if flags.iter().any(|flag| *flag == key) {
            let value = args
                .get(i + 1)
                .filter(|value| {
                    matches!(
                        value.to_ascii_lowercase().as_str(),
                        "true" | "false" | "1" | "0"
                    )
                })
                .cloned()
                .unwrap_or_else(|| "true".to_string());
            if out.insert(key.clone(), value).is_some() {
                bail!("duplicate option: {key}");
            }
            i += if args
                .get(i + 1)
                .is_some_and(|value| {
                    matches!(
                        value.to_ascii_lowercase().as_str(),
                        "true" | "false" | "1" | "0"
                    )
                })
            {
                2
            } else {
                1
            };
            continue;
        }
        let value = args
            .get(i + 1)
            .with_context(|| format!("{key} needs a value"))?;
        out.insert(key.clone(), value.clone());
        i += 2;
    }
    Ok(out)
}

fn cmd_import(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let db = opts.get("--db").context("missing --db")?;
    let evidence = opts.get("--evidence").context("missing --evidence")?;
    let mut store = Store::open(std::path::Path::new(db))?;
    let edges = io::read_edges(std::path::Path::new(evidence))?;
    let total = edges.len();
    let inserted = store.insert_edges(&edges)?;
    println!("read {total} evidence rows; inserted {inserted} new rows into {db}");
    Ok(())
}

fn cmd_export(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let db = opts.get("--db").context("missing --db")?;
    let out = opts.get("--out").context("missing --out")?;
    let store = Store::open(std::path::Path::new(db))?;
    let mut edges = store.export_edges()?;
    let n = edges.len();
    io::write_edges(std::path::Path::new(out), &mut edges)?;
    println!("exported {n} evidence rows to {out}");
    Ok(())
}

fn cmd_stats(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let db = opts.get("--db").context("missing --db")?;
    let store = Store::open(std::path::Path::new(db))?;
    let stats = store.stats()?;
    println!("evidence_level\tevidence_type\tcount");
    for (level, kind, count) in stats {
        println!("{level}\t{kind}\t{count}");
    }
    Ok(())
}

fn cmd_validate_db(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let db = opts.get("--db").context("missing --db")?;
    let conn = Connection::open_with_flags(db, OpenFlags::SQLITE_OPEN_READ_ONLY)
        .with_context(|| format!("open {db}"))?;
    let pairs_required = [
        "tf_id",
        "tf_symbol",
        "target_id",
        "target_symbol",
        "n_sites",
        "n_motifs",
        "best_p",
        "best_score",
        "motifs",
        "evidence_levels",
    ];
    for (table, required) in [
        ("hits", EVIDENCE_COLUMNS.as_slice()),
        ("pairs", pairs_required.as_slice()),
    ] {
        let exists: i64 = conn.query_row(
            "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?1",
            params![table],
            |r| r.get(0),
        )?;
        if exists == 0 {
            bail!("database has no {table} table");
        }
        let mut stmt = conn.prepare(&format!("PRAGMA table_info({table})"))?;
        let names: HashSet<String> = stmt
            .query_map([], |r| r.get::<_, String>(1))?
            .collect::<std::result::Result<_, _>>()?;
        let missing: Vec<&str> = required
            .iter()
            .copied()
            .filter(|name| !names.contains(*name))
            .collect();
        if !missing.is_empty() {
            bail!(
                "database table {table} is missing columns: {}",
                missing.join(", ")
            );
        }
    }
    let mut integrity = conn.prepare("PRAGMA integrity_check")?;
    for message in integrity.query_map([], |r| r.get::<_, String>(0))? {
        let message = message?;
        if message != "ok" {
            bail!("database integrity check failed: {message}");
        }
    }
    for (table, index, column) in [
        ("pairs", "pairs_tf", "tf_id"),
        ("pairs", "pairs_target", "target_id"),
        ("hits", "hits_tf", "tf_id"),
        ("hits", "hits_target", "target_id"),
        ("hits", "hits_motif", "motif_id"),
        ("hits", "hits_p", "p_value"),
    ] {
        let mut indexes = conn.prepare(&format!("PRAGMA index_list({table})"))?;
        let found = indexes
            .query_map([], |r| {
                Ok((
                    r.get::<_, String>(1)?,
                    r.get::<_, i64>(4)?, // partial
                ))
            })?
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let Some((_, partial)) = found.iter().find(|(name, _)| name == index) else {
            bail!("database table {table} is missing index {index}");
        };
        if *partial != 0 {
            bail!("database index {index} must not be a partial index");
        }
        let mut info = conn.prepare(&format!("PRAGMA index_info({index})"))?;
        let columns: Vec<String> = info
            .query_map([], |r| r.get::<_, String>(2))?
            .collect::<std::result::Result<_, _>>()?;
        if columns.len() != 1 || columns[0] != column {
            bail!(
                "database index {index} must index {table}.{column}, found {}",
                columns.join(", ")
            );
        }
    }
    let bad_hits: i64 = conn.query_row(
        "SELECT COUNT(*) FROM hits
         WHERE typeof(tf_id) <> 'text'
            OR typeof(tf_symbol) <> 'text'
            OR typeof(target_id) <> 'text'
            OR typeof(target_symbol) <> 'text'
            OR typeof(evidence_type) <> 'text'
            OR typeof(evidence_level) <> 'text'
            OR typeof(method) <> 'text'
            OR typeof(source) <> 'text'
            OR typeof(pmid) <> 'text'
            OR typeof(condition) <> 'text'
            OR (promoter_start IS NOT NULL AND typeof(promoter_start) <> 'integer')
            OR (promoter_end IS NOT NULL AND typeof(promoter_end) <> 'integer')
            OR typeof(motif_id) <> 'text'
            OR (hit_start IS NOT NULL AND typeof(hit_start) <> 'integer')
            OR (hit_end IS NOT NULL AND typeof(hit_end) <> 'integer')
            OR typeof(strand) <> 'text'
            OR (score IS NOT NULL AND typeof(score) NOT IN ('integer', 'real'))
            OR (p_value IS NOT NULL AND typeof(p_value) NOT IN ('integer', 'real'))
            OR typeof(notes) <> 'text'",
        [],
        |r| r.get(0),
    )?;
    if bad_hits != 0 {
        bail!("database hits table has {bad_hits} row(s) with invalid column types");
    }
    let bad_pairs: i64 = conn.query_row(
        "SELECT COUNT(*) FROM pairs
         WHERE typeof(tf_id) <> 'text'
            OR typeof(tf_symbol) <> 'text'
            OR typeof(target_id) <> 'text'
            OR typeof(target_symbol) <> 'text'
            OR typeof(n_sites) <> 'integer'
            OR typeof(n_motifs) <> 'integer'
            OR (best_p IS NOT NULL AND typeof(best_p) NOT IN ('integer', 'real'))
            OR (best_score IS NOT NULL AND typeof(best_score) NOT IN ('integer', 'real'))
            OR (motifs IS NOT NULL AND typeof(motifs) <> 'text')
            OR (evidence_levels IS NOT NULL AND typeof(evidence_levels) <> 'text')",
        [],
        |r| r.get(0),
    )?;
    if bad_pairs != 0 {
        bail!("database pairs table has {bad_pairs} row(s) with invalid column types");
    }
    println!("valid TF-target prior database");
    Ok(())
}

fn cmd_classify_tfs(args: &[String]) -> Result<()> {
    use regutrace_core::{hmmer, tf_classification};

    let opts = parse_options(args)?;
    for option in opts.keys() {
        if !matches!(option.as_str(), "--domtblout" | "--rules" | "--out") {
            bail!("unknown option: {option}");
        }
    }
    let domtblout = Path::new(opts.get("--domtblout").context("missing --domtblout")?);
    let rules = Path::new(opts.get("--rules").context("missing --rules")?);
    let out = Path::new(opts.get("--out").context("missing --out")?);
    for input in ["--domtblout", "--rules"] {
        if same_path(&opts[input], &opts["--out"]) {
            bail!("{input} and --out must be different files");
        }
    }
    let pack: tf_classification::RulePack = serde_json::from_reader(fasta::open_reader(rules)?)
        .with_context(|| format!("read TF rules {}", rules.display()))?;
    let hits = hmmer::read_hmmsearch_domains(BufReader::new(fasta::open_reader(domtblout)?))?;
    let (classifications, unknown_models) = tf_classification::classify(&pack, &hits)?;

    // Keep protein identities until a separate gene mapping step can validate them.
    let tmp = temp_output_path(out);
    let result = (|| -> Result<()> {
        let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(&tmp)?;
        writer.write_record(["protein_id", "family", "subfamily", "rules_version", "domains_json"])?;
        for classification in &classifications {
            let domains = serde_json::to_string(&classification.domains)?;
            for rule in &classification.families {
                writer.write_record([
                    classification.protein_id,
                    &rule.family,
                    &rule.subfamily,
                    &pack.version,
                    &domains,
                ])?;
            }
        }
        writer.flush()?;
        drop(writer);
        publish_temp(&tmp, out)
    })();
    if result.is_err() {
        cleanup_output_temp(&tmp);
    }
    result?;
    if !unknown_models.is_empty() {
        eprintln!("models absent from the rule pack: {}", unknown_models.join(", "));
    }
    println!("classified {} protein(s)", classifications.len());
    Ok(())
}

fn cmd_identify_tfs(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let proteins = opts.get("--proteins").context("missing --proteins")?;
    let hmm = opts.get("--hmm").context("missing --hmm")?;
    let out = opts.get("--out").context("missing --out")?;
    let evalue = opts
        .get("--evalue")
        .map(|value| value.parse::<f64>())
        .transpose()
        .context("--evalue must be a number")?
        .unwrap_or(1e-5);
    if !evalue.is_finite() || evalue <= 0.0 {
        bail!("--evalue must be a positive finite number");
    }
    let min_coverage = opts
        .get("--min-coverage")
        .map(|value| value.parse::<f64>())
        .transpose()
        .context("--min-coverage must be a number")?
        .unwrap_or(0.5);
    if !(0.0..=1.0).contains(&min_coverage) {
        bail!("--min-coverage must be between 0 and 1");
    }
    let family_map = match opts.get("--family-map") {
        Some(path) => read_family_map(Path::new(path))?,
        None => HashMap::new(),
    };

    let fa = fasta::read_fasta(Path::new(proteins))?;
    if fa.order.is_empty() {
        bail!("protein FASTA is empty: {proteins}");
    }
    let seqs: Vec<Seq> = fa
        .order
        .iter()
        .filter_map(|name| {
            fa.seqs.get(name).map(|seq| {
                Seq::from_amino(name.clone(), String::new(), &String::from_utf8_lossy(seq))
            })
        })
        .collect();
    let annotator = HmmAnnotator::from_hmm_file(hmm)
        .map_err(|err| anyhow::anyhow!("read HMM database {hmm}: {err}"))?
        .with_cutoff(Cutoff::Evalue(evalue));
    let hits = annotator.search(&seqs);

    let mut best: HashMap<String, (String, f64, f64)> = HashMap::new();
    for hit in hits {
        let coverage = if hit.model_len == 0 {
            0.0
        } else {
            (hit.hmm_to - hit.hmm_from + 1).max(0) as f64 / hit.model_len as f64
        };
        if hit.best_dom_evalue > evalue || coverage < min_coverage {
            continue;
        }
        let candidate = (hit.query_name.clone(), hit.best_dom_evalue, coverage);
        match best.get(&hit.target_name) {
            Some((_, old_evalue, old_coverage))
                if (*old_evalue < candidate.1)
                    || (*old_evalue == candidate.1 && *old_coverage >= candidate.2) => {}
            _ => {
                best.insert(hit.target_name, candidate);
            }
        }
    }
    if best.is_empty() {
        bail!("no TF-domain hit passed --evalue {evalue} and --min-coverage {min_coverage}");
    }

    let identified = best.len();
    let tmp = temp_output_path(Path::new(out));
    let result = (|| -> Result<()> {
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&tmp)
            .with_context(|| format!("create {}", tmp.display()))?;
        writer.write_record([
            "gene_id",
            "tf_symbol",
            "family",
            "identification_status",
            "hmm_model",
            "evalue",
            "coverage",
            "source",
        ])?;
        let mut genes: Vec<_> = best.into_iter().collect();
        genes.sort_by(|a, b| a.0.cmp(&b.0));
        for (gene, (model, hit_evalue, coverage)) in genes {
            let family = family_map.get(&model).cloned().unwrap_or_else(|| model.clone());
            writer.write_record([
                gene.as_str(),
                gene.as_str(),
                family.as_str(),
                "identified",
                model.as_str(),
                &format!("{hit_evalue:.6e}"),
                &format!("{coverage:.3}"),
                "regutrace-identify-tfs",
            ])?;
        }
        writer.flush()?;
        drop(writer);
        Ok(())
    })();
    if let Err(err) = result {
        cleanup_output_temp(&tmp);
        return Err(err);
    }
    if let Err(err) = publish_temp(&tmp, Path::new(out)) {
        cleanup_output_temp(&tmp);
        return Err(err);
    }
    println!("identified {identified} TF gene(s) -> {out}");
    Ok(())
}

fn read_family_map(path: &Path) -> Result<HashMap<String, String>> {
    let reader = fasta::open_reader(path)?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(reader);
    let headers = rdr
        .headers()
        .with_context(|| format!("read header {}", path.display()))?
        .clone();
    let model = headers
        .iter()
        .position(|header| {
            matches!(
                header.trim().to_ascii_lowercase().as_str(),
                "hmm_model" | "model" | "hmm" | "name"
            )
        })
        .with_context(|| format!("family map missing model column: {}", path.display()))?;
    let family = headers
        .iter()
        .position(|header| {
            matches!(
                header.trim().to_ascii_lowercase().as_str(),
                "family" | "tf_family"
            )
        })
        .with_context(|| format!("family map missing family column: {}", path.display()))?;
    let mut out = HashMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        let model = rec.get(model).unwrap_or("").trim();
        let family = rec.get(family).unwrap_or("").trim();
        if !model.is_empty() && !family.is_empty() {
            out.insert(model.to_string(), family.to_string());
        }
    }
    Ok(out)
}

fn cmd_prior(args: &[String]) -> Result<()> {
    let opts = parse_options(args)?;
    let db = opts.get("--db").context("missing --db")?;
    let targets = opts.get("--targets").context("missing --targets")?;
    let out = opts.get("--out").context("missing --out")?;
    if same_path(targets, out) {
        bail!("--targets and --out must be different files");
    }
    if same_path(db, out) {
        bail!("--db and --out must be different files");
    }

    let target_file = File::open(targets).with_context(|| format!("open {targets}"))?;
    let mut ids: Vec<String> = Vec::new();
    let mut first_nonempty = true;
    for line in BufReader::new(target_file).lines() {
        let id = line?.trim().to_string();
        if id.is_empty() {
            continue;
        }
        if first_nonempty {
            first_nonempty = false;
            if matches!(
                id.to_ascii_lowercase().as_str(),
                "gene_id" | "target_id" | "gene" | "target"
            ) {
                continue;
            }
        }
        ids.push(id);
    }
    if ids.is_empty() {
        bail!("targets file is empty");
    }

    let mut conn = Connection::open_with_flags(db, OpenFlags::SQLITE_OPEN_READ_ONLY)
        .with_context(|| format!("open {db}"))?;
    let tx = conn.transaction()?;
    let has_pairs: i64 = tx.query_row(
        "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='pairs'",
        [],
        |r| r.get(0),
    )?;
    if has_pairs == 0 {
        bail!("database has no pairs table; run `regutrace convert --format sqlite` first");
    }
    tx.execute_batch("CREATE TEMP TABLE IF NOT EXISTS target_ids(id TEXT PRIMARY KEY);")?;
    {
        let mut insert = tx.prepare("INSERT OR IGNORE INTO target_ids(id) VALUES (?1)")?;
        for id in &ids {
            insert.execute(params![id])?;
        }
    }

    let tmp = temp_output_path(Path::new(out));
    let write_result = (|| -> Result<()> {
        let mut stmt = tx.prepare(
            "SELECT p.tf_id,
                    COALESCE(MAX(NULLIF(p.tf_symbol, '')), '') AS tf_symbol,
                    COUNT(DISTINCT p.target_id) AS module_targets,
                    SUM(p.n_sites) AS n_sites,
                    MIN(p.best_p) AS best_p,
                    GROUP_CONCAT(DISTINCT p.evidence_levels) AS evidence_levels
             FROM pairs p
             JOIN target_ids t ON t.id = p.target_id
             GROUP BY p.tf_id
             ORDER BY module_targets DESC, best_p ASC, p.tf_id ASC",
        )?;
        let rows = stmt.query_map([], |r| {
            Ok((
                r.get::<_, String>(0)?,
                r.get::<_, String>(1)?,
                r.get::<_, i64>(2)?,
                r.get::<_, Option<i64>>(3)?,
                r.get::<_, Option<f64>>(4)?,
                r.get::<_, Option<String>>(5)?,
            ))
        })?;
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&tmp)
            .with_context(|| format!("create {}", tmp.display()))?;
        writer.write_record([
            "tf_id",
            "tf_symbol",
            "module_targets",
            "n_sites",
            "best_p",
            "evidence_levels",
        ])?;
        for row in rows {
            let (tf_id, tf_symbol, module_targets, n_sites, best_p, levels) = row?;
            writer.write_record([
                tf_id,
                tf_symbol,
                module_targets.to_string(),
                n_sites.unwrap_or(0).to_string(),
                best_p.map(|v| format!("{v:.6e}")).unwrap_or_default(),
                normalize_levels(levels),
            ])?;
        }
        writer.flush()?;
        drop(writer);
        drop(stmt);
        tx.commit()?;
        Ok(())
    })();
    if let Err(err) = write_result {
        cleanup_output_temp(&tmp);
        return Err(err);
    }
    if let Err(err) = publish_temp(&tmp, Path::new(out)) {
        cleanup_output_temp(&tmp);
        return Err(err);
    }
    Ok(())
}

fn normalize_levels(raw: Option<String>) -> String {
    let mut levels = BTreeSet::new();
    if let Some(raw) = raw {
        for level in raw.split(',') {
            let level = level.trim();
            if !level.is_empty() {
                levels.insert(level.to_string());
            }
        }
    }
    levels.into_iter().collect::<Vec<_>>().join(",")
}

const INTEGER_COLUMNS: [&str; 4] = [
    "promoter_start",
    "promoter_end",
    "hit_start",
    "hit_end",
];

const FLOAT_COLUMNS: [&str; 2] = ["score", "p_value"];

const EVIDENCE_COLUMNS: [&str; 19] = [
    "tf_id",
    "tf_symbol",
    "target_id",
    "target_symbol",
    "evidence_type",
    "evidence_level",
    "method",
    "source",
    "pmid",
    "condition",
    "promoter_start",
    "promoter_end",
    "motif_id",
    "hit_start",
    "hit_end",
    "strand",
    "score",
    "p_value",
    "notes",
];

const SQLITE_SCHEMA: &str = r#"
PRAGMA user_version = 1;
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS hits (
    tf_id TEXT NOT NULL,
    tf_symbol TEXT NOT NULL DEFAULT '',
    target_id TEXT NOT NULL,
    target_symbol TEXT NOT NULL DEFAULT '',
    evidence_type TEXT NOT NULL DEFAULT '',
    evidence_level TEXT NOT NULL DEFAULT '',
    method TEXT NOT NULL DEFAULT '',
    source TEXT NOT NULL DEFAULT '',
    pmid TEXT NOT NULL DEFAULT '',
    condition TEXT NOT NULL DEFAULT '',
    promoter_start INTEGER,
    promoter_end INTEGER,
    motif_id TEXT NOT NULL DEFAULT '',
    hit_start INTEGER,
    hit_end INTEGER,
    strand TEXT NOT NULL DEFAULT '',
    score REAL,
    p_value REAL,
    notes TEXT NOT NULL DEFAULT ''
);
"#;

/// Convert the canonical evidence TSV (or a CSV/TSV.GZ/CSV.GZ file) into a
/// full SQLite database, or stream it to CSV/TSV/JSONL. This is intentionally
/// row-oriented: the full hit table is retained, while `pairs` is materialized
/// for module-level analysis.
fn cmd_convert(args: &[String]) -> Result<()> {
    let opts = parse_options_with_flags(args, &["--gzip"])?;
    let input = opts.get("--input").context("missing --input")?;
    let out = opts.get("--out").context("missing --out")?;
    let format = opts.get("--format").context("missing --format")?.as_str();
    if !matches!(format, "sqlite" | "csv" | "tsv" | "jsonl") {
        bail!("--format must be one of sqlite, csv, tsv, jsonl");
    }
    let gzip = opts
        .get("--gzip")
        .is_some_and(|value| value.eq_ignore_ascii_case("true") || value == "1");
    if gzip && format == "sqlite" {
        bail!("--gzip is not valid with --format sqlite");
    }
    if same_path(input, out) {
        bail!("--input and --out must be different files");
    }
    if input.to_ascii_lowercase().ends_with(".jsonl")
        || input.to_ascii_lowercase().ends_with(".jsonl.gz")
    {
        bail!("JSONL input is not supported by convert; import a TSV/CSV or SQLite file");
    }

    let delimiter = if input.to_ascii_lowercase().ends_with(".csv")
        || input.to_ascii_lowercase().ends_with(".csv.gz")
    {
        b','
    } else {
        b'\t'
    };
    let reader = fasta::open_reader(Path::new(input))?;
    let mut rdr = ReaderBuilder::new()
        .delimiter(delimiter)
        .flexible(false)
        .from_reader(reader);
    let headers = rdr.headers()?.clone();
    validate_headers(&headers)?;

    let tmp = temp_output_path(Path::new(out));
    let result = match format {
        "sqlite" => convert_to_sqlite(&mut rdr, &headers, &tmp),
        _ => convert_to_text(&mut rdr, &headers, &tmp, format, gzip),
    };
    let rows = match result {
        Ok(rows) => rows,
        Err(err) => {
            cleanup_output_temp(&tmp);
            return Err(err);
        }
    };
    if let Err(err) = publish_temp(&tmp, Path::new(out)) {
        cleanup_output_temp(&tmp);
        return Err(err);
    }
    println!(
        "converted {rows} evidence rows to {out} ({format}{})",
        if gzip { ", gzip" } else { "" }
    );
    Ok(())
}

fn validate_headers(headers: &StringRecord) -> Result<()> {
    let mut seen = HashSet::new();
    for header in headers {
        let name = header.trim();
        if name.is_empty() {
            bail!("input contains an empty header column");
        }
        let key = name.to_ascii_lowercase();
        if !seen.insert(key) {
            bail!("input contains duplicate header column: {name}");
        }
    }
    if !seen.contains("tf_id") || !seen.contains("target_id") {
        bail!("input needs tf_id and target_id columns");
    }
    Ok(())
}

fn cell<'a>(headers: &StringRecord, rec: &'a StringRecord, name: &str) -> &'a str {
    headers
        .iter()
        .position(|h| h.trim().eq_ignore_ascii_case(name))
        .and_then(|i| rec.get(i))
        .unwrap_or("")
}

fn is_integer_column(name: &str) -> bool {
    INTEGER_COLUMNS
        .iter()
        .any(|column| column.eq_ignore_ascii_case(name.trim()))
}

fn is_float_column(name: &str) -> bool {
    FLOAT_COLUMNS
        .iter()
        .any(|column| column.eq_ignore_ascii_case(name.trim()))
}

fn validate_record(headers: &StringRecord, rec: &StringRecord, row: u64) -> Result<()> {
    if cell(headers, rec, "tf_id").trim().is_empty()
        || cell(headers, rec, "target_id").trim().is_empty()
    {
        bail!("row {row}: empty tf_id or target_id");
    }
    for name in INTEGER_COLUMNS {
        optional_i64(headers, rec, name, row)?;
    }
    for name in FLOAT_COLUMNS {
        optional_f64(headers, rec, name, row)?;
    }
    Ok(())
}

fn optional_i64(headers: &StringRecord, rec: &StringRecord, name: &str, row: u64) -> Result<Option<i64>> {
    let v = cell(headers, rec, name).trim();
    if v.is_empty() {
        return Ok(None);
    }
    v.parse::<i64>()
        .map(Some)
        .with_context(|| format!("row {row}: invalid integer in {name}"))
}

fn optional_f64(headers: &StringRecord, rec: &StringRecord, name: &str, row: u64) -> Result<Option<f64>> {
    let v = cell(headers, rec, name).trim();
    if v.is_empty() {
        return Ok(None);
    }
    let n = v
        .parse::<f64>()
        .with_context(|| format!("row {row}: invalid number in {name}"))?;
    if !n.is_finite() {
        bail!("row {row}: non-finite number in {name}");
    }
    Ok(Some(n))
}

fn record_values(headers: &StringRecord, rec: &StringRecord, row: u64) -> Result<Vec<Value>> {
    if cell(headers, rec, "tf_id").trim().is_empty()
        || cell(headers, rec, "target_id").trim().is_empty()
    {
        bail!("row {row}: empty tf_id or target_id");
    }
    let mut vals = Vec::with_capacity(EVIDENCE_COLUMNS.len());
    for name in EVIDENCE_COLUMNS {
        let v = match name {
            "promoter_start" | "promoter_end" | "hit_start" | "hit_end" => {
                optional_i64(headers, rec, name, row)?.map(Value::Integer).unwrap_or(Value::Null)
            }
            "score" | "p_value" => {
                optional_f64(headers, rec, name, row)?.map(Value::Real).unwrap_or(Value::Null)
            }
            _ => Value::Text(cell(headers, rec, name).to_string()),
        };
        vals.push(v);
    }
    Ok(vals)
}

fn convert_to_sqlite<R: std::io::Read>(
    rdr: &mut csv::Reader<R>,
    headers: &StringRecord,
    out: &Path,
) -> Result<u64> {
    let mut conn = Connection::open(out)?;
    conn.execute_batch("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;")?;
    conn.execute_batch(SQLITE_SCHEMA)?;
    let tx = conn.transaction()?;
    let mut rows = 0u64;
    {
        let mut stmt = tx.prepare(
            "INSERT INTO hits (
                tf_id, tf_symbol, target_id, target_symbol, evidence_type, evidence_level,
                method, source, pmid, condition, promoter_start, promoter_end, motif_id,
                hit_start, hit_end, strand, score, p_value, notes
            ) VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18,?19)",
        )?;
        for rec in rdr.records() {
            let rec = rec?;
            rows += 1;
            stmt.execute(params_from_iter(record_values(headers, &rec, rows)?))?;
        }
    }
    tx.commit()?;
    conn.execute_batch(
        "CREATE TABLE pairs AS
         SELECT tf_id, tf_symbol, target_id, target_symbol,
                COUNT(*) AS n_sites,
                COUNT(DISTINCT CASE WHEN motif_id <> '' THEN motif_id END) AS n_motifs,
                MIN(p_value) AS best_p,
                MAX(score) AS best_score,
                GROUP_CONCAT(DISTINCT CASE WHEN motif_id <> '' THEN motif_id END) AS motifs,
                GROUP_CONCAT(DISTINCT CASE WHEN evidence_level <> '' THEN evidence_level END) AS evidence_levels
         FROM hits
         GROUP BY tf_id, target_id;
         CREATE INDEX hits_tf ON hits(tf_id);
         CREATE INDEX hits_target ON hits(target_id);
         CREATE INDEX hits_motif ON hits(motif_id);
         CREATE INDEX hits_p ON hits(p_value);
         CREATE INDEX pairs_tf ON pairs(tf_id);
         CREATE INDEX pairs_target ON pairs(target_id);",
    )?;
    let pairs: i64 = conn.query_row("SELECT COUNT(*) FROM pairs", [], |r| r.get(0))?;
    let tfs: i64 = conn.query_row("SELECT COUNT(DISTINCT tf_id) FROM hits", [], |r| r.get(0))?;
    let targets: i64 = conn.query_row("SELECT COUNT(DISTINCT target_id) FROM hits", [], |r| r.get(0))?;
    for (key, value) in [
        ("schema_version", "1".to_string()),
        ("tool_version", env!("CARGO_PKG_VERSION").to_string()),
        ("format", "regutrace-evidence".to_string()),
        ("n_hits", rows.to_string()),
        ("n_pairs", pairs.to_string()),
        ("n_tfs", tfs.to_string()),
        ("n_targets", targets.to_string()),
        ("created_at", time::now_utc()),
    ] {
        conn.execute(
            "INSERT OR REPLACE INTO meta(key, value) VALUES (?1, ?2)",
            params![key, value],
        )?;
    }
    drop(conn);
    Ok(rows)
}

enum OutputSink {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
}

impl OutputSink {
    fn new(file: File, gzip: bool) -> Self {
        let buffered = BufWriter::new(file);
        if gzip {
            Self::Gzip(GzEncoder::new(buffered, Compression::default()))
        } else {
            Self::Plain(buffered)
        }
    }

    fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut writer) => {
                writer.flush()?;
                Ok(())
            }
            Self::Gzip(writer) => {
                let mut inner = writer.finish()?;
                inner.flush()?;
                Ok(())
            }
        }
    }
}

impl Write for OutputSink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buf),
            Self::Gzip(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
        }
    }
}

fn json_cell(
    headers: &StringRecord,
    rec: &StringRecord,
    name: &str,
    row: u64,
) -> Result<JsonValue> {
    if is_integer_column(name) {
        return Ok(optional_i64(headers, rec, name, row)?
            .map(|v| JsonValue::Number(v.into()))
            .unwrap_or(JsonValue::Null));
    }
    if is_float_column(name) {
        return Ok(optional_f64(headers, rec, name, row)?
            .and_then(serde_json::Number::from_f64)
            .map(JsonValue::Number)
            .unwrap_or(JsonValue::Null));
    }
    Ok(JsonValue::String(cell(headers, rec, name).to_string()))
}

fn convert_to_text<R: std::io::Read>(
    rdr: &mut csv::Reader<R>,
    headers: &StringRecord,
    out: &Path,
    format: &str,
    gzip: bool,
) -> Result<u64> {
    let file = File::create(out)?;
    let sink = OutputSink::new(file, gzip);
    match format {
        "csv" | "tsv" => {
            let delimiter = if format == "csv" { b',' } else { b'\t' };
            let mut wtr = WriterBuilder::new().delimiter(delimiter).from_writer(sink);
            wtr.write_record(headers)?;
            let mut rows = 0u64;
            for rec in rdr.records() {
                let rec = rec?;
                rows += 1;
                validate_record(headers, &rec, rows)?;
                wtr.write_record(&rec)?;
            }
            wtr.flush()?;
            let sink = wtr.into_inner().map_err(|err| err.into_error())?;
            sink.finish()?;
            Ok(rows)
        }
        "jsonl" => {
            let mut sink = sink;
            let mut rows = 0u64;
            for rec in rdr.records() {
                let rec = rec?;
                rows += 1;
                if cell(headers, &rec, "tf_id").trim().is_empty()
                    || cell(headers, &rec, "target_id").trim().is_empty()
                {
                    bail!("row {rows}: empty tf_id or target_id");
                }
                let mut obj = serde_json::Map::new();
                for name in headers.iter() {
                    obj.insert(
                        name.trim().to_string(),
                        json_cell(headers, &rec, name, rows)?,
                    );
                }
                serde_json::to_writer(&mut sink, &obj)?;
                sink.write_all(b"\n")?;
            }
            sink.finish()?;
            Ok(rows)
        }
        _ => unreachable!(),
    }
}

fn temp_output_path(out: &Path) -> PathBuf {
    let parent = out
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = out.file_name().and_then(|name| name.to_str()).unwrap_or("output");
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    parent.join(format!(".{name}.regutrace-{}-{nanos}.tmp", std::process::id()))
}

fn cleanup_output_temp(path: &Path) {
    let _ = fs::remove_file(path);
    for suffix in ["-journal", "-wal", "-shm"] {
        let mut sidecar = path.as_os_str().to_os_string();
        sidecar.push(suffix);
        let _ = fs::remove_file(PathBuf::from(sidecar));
    }
}

fn publish_temp(tmp: &Path, out: &Path) -> Result<()> {
    if out.is_dir() {
        bail!("output path is a directory: {}", out.display());
    }
    #[cfg(windows)]
    if out.exists() {
        return replace_file_windows(tmp, out)
            .with_context(|| format!("replace {}", out.display()));
    }
    fs::rename(tmp, out).with_context(|| format!("publish {}", out.display()))
}

#[cfg(windows)]
fn replace_file_windows(replacement: &Path, replaced: &Path) -> std::io::Result<()> {
    use std::os::windows::ffi::OsStrExt;

    unsafe extern "system" {
        fn ReplaceFileW(
            replaced_file_name: *const u16,
            replacement_file_name: *const u16,
            backup_file_name: *const u16,
            replace_flags: u32,
            exclude: *mut std::ffi::c_void,
            reserved: *mut std::ffi::c_void,
        ) -> i32;
    }

    let replaced: Vec<u16> = replaced
        .as_os_str()
        .encode_wide()
        .chain(std::iter::once(0))
        .collect();
    let replacement: Vec<u16> = replacement
        .as_os_str()
        .encode_wide()
        .chain(std::iter::once(0))
        .collect();
    let ok = unsafe {
        ReplaceFileW(
            replaced.as_ptr(),
            replacement.as_ptr(),
            std::ptr::null(),
            0,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        )
    };
    if ok == 0 {
        Err(std::io::Error::last_os_error())
    } else {
        Ok(())
    }
}

fn same_path(left: &str, right: &str) -> bool {
    let left = Path::new(left);
    let right = Path::new(right);
    if left == right {
        return true;
    }
    if same_file_identity(left, right) {
        return true;
    }
    match (fs::canonicalize(left), fs::canonicalize(right)) {
        (Ok(left), Ok(right)) => left == right,
        _ => false,
    }
}

#[cfg(unix)]
fn same_file_identity(left: &Path, right: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;
    let (Ok(left), Ok(right)) = (fs::metadata(left), fs::metadata(right)) else {
        return false;
    };
    left.dev() == right.dev() && left.ino() == right.ino()
}

#[cfg(windows)]
fn same_file_identity(left: &Path, right: &Path) -> bool {
    use std::os::windows::io::AsRawHandle;

    #[repr(C)]
    struct FileTime {
        low: u32,
        high: u32,
    }
    #[repr(C)]
    struct FileInformation {
        attributes: u32,
        creation_time: FileTime,
        last_access_time: FileTime,
        last_write_time: FileTime,
        volume_serial_number: u32,
        size_high: u32,
        size_low: u32,
        links: u32,
        index_high: u32,
        index_low: u32,
    }
    unsafe extern "system" {
        fn GetFileInformationByHandle(
            file: *mut std::ffi::c_void,
            information: *mut FileInformation,
        ) -> i32;
    }
    fn identity(path: &Path) -> Option<(u32, u32, u32)> {
        let file = File::open(path).ok()?;
        let mut information = std::mem::MaybeUninit::<FileInformation>::zeroed();
        let ok = unsafe {
            GetFileInformationByHandle(
                file.as_raw_handle() as *mut std::ffi::c_void,
                information.as_mut_ptr(),
            )
        };
        if ok == 0 {
            return None;
        }
        let information = unsafe { information.assume_init() };
        Some((
            information.volume_serial_number,
            information.index_high,
            information.index_low,
        ))
    }
    match (identity(left), identity(right)) {
        (Some(left), Some(right)) => left == right,
        _ => false,
    }
}

#[cfg(not(any(unix, windows)))]
fn same_file_identity(_left: &Path, _right: &Path) -> bool {
    false
}

mod time {
    pub fn now_utc() -> String {
        // The database is a local metadata artifact; a fixed UTC timestamp is
        // clearer than pulling a date-time dependency into the CLI.
        use std::time::{SystemTime, UNIX_EPOCH};
        let secs = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        secs.to_string()
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let res = match args.first().map(String::as_str) {
        None | Some("-h") | Some("--help") => {
            usage();
            return;
        }
        Some("scan") => cmd_scan(&args[1..]),
        Some("validate") => cmd_validate(&args[1..]),
        Some("pipeline") => pipeline::run(&args[1..]),
        Some("trace") => trace::run(&args[1..]),
        Some("import") => cmd_import(&args[1..]),
        Some("export") => cmd_export(&args[1..]),
        Some("stats") => cmd_stats(&args[1..]),
        Some("convert") => cmd_convert(&args[1..]),
        Some("prior") => cmd_prior(&args[1..]),
        Some("identify-plant-tfs") => plant_tfs::run(&args[1..]),
        Some("map-tf-motifs") => motif_transfer::run(&args[1..]),
        Some("classify-tfs") => cmd_classify_tfs(&args[1..]),
        Some("identify-tfs") => cmd_identify_tfs(&args[1..]),
        Some("validate-db") => cmd_validate_db(&args[1..]),
        Some(other) => Err(anyhow::anyhow!("unknown command: {other}")),
    };
    if let Err(e) = res {
        eprintln!("error: {e:#}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convert_and_prior() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-convert-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let tsv = dir.join("evidence.tsv");
        let db = dir.join("evidence.db");
        let targets = dir.join("targets.txt");
        let upstream = dir.join("upstream.tsv");
        std::fs::write(
            &tsv,
            "tf_id\ttf_symbol\ttarget_id\ttarget_symbol\tevidence_type\tevidence_level\tmethod\tsource\tpmid\tcondition\tpromoter_start\tpromoter_end\tmotif_id\thit_start\thit_end\tstrand\tscore\tp_value\tnotes\n\
TF1\tTF1\tG1\tGENE1\tprediction_motif\tD\tJASPAR\t\t\theat\t101\t1100\tMA1\t10\t16\t+\t8.5\t1e-6\ta\n\
TF1\tTF1\tG1\tGENE1\tprediction_motif\tD\tJASPAR\t\t\theat\t101\t1100\tMA2\t20\t26\t-\t7.0\t2e-5\tb\n\
TF2\tTF2\tG2\tGENE2\tprediction_motif\tD\tJASPAR\t\t\theat\t201\t1200\tMA1\t15\t21\t+\t6.0\t3e-4\tc\n\
TF2\tTF2\tG2\tGENE2\tprediction_motif\tA\tJASPAR\t\t\theat\t201\t1200\tMA2\t25\t31\t-\t5.0\t4e-4\td\n\
TF2\tTF2\tG3\tGENE3\tprediction_motif\tD\tJASPAR\t\t\theat\t301\t1300\tMA3\t30\t36\t+\t4.0\t5e-5\te\n",
        )
        .unwrap();
        cmd_convert(&[
            "--input".into(),
            tsv.to_string_lossy().into_owned(),
            "--out".into(),
            db.to_string_lossy().into_owned(),
            "--format".into(),
            "sqlite".into(),
        ])
        .unwrap();
        cmd_convert(&[
            "--input".into(),
            tsv.to_string_lossy().into_owned(),
            "--out".into(),
            db.to_string_lossy().into_owned(),
            "--format".into(),
            "sqlite".into(),
        ])
        .unwrap();
        let conn = Connection::open(&db).unwrap();
        let hits: i64 = conn.query_row("SELECT COUNT(*) FROM hits", [], |r| r.get(0)).unwrap();
        let pairs: i64 = conn.query_row("SELECT COUNT(*) FROM pairs", [], |r| r.get(0)).unwrap();
        assert_eq!(hits, 5);
        assert_eq!(pairs, 3);
        drop(conn);
        cmd_validate_db(&["--db".into(), db.to_string_lossy().into_owned()]).unwrap();
        let conn = Connection::open(&db).unwrap();
        let (n_sites, n_motifs, best_p, motifs, levels): (i64, i64, f64, String, String) = conn
            .query_row(
                "SELECT n_sites, n_motifs, best_p, motifs, evidence_levels
                 FROM pairs WHERE tf_id='TF2' AND target_id='G2'",
                [],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?, r.get(4)?)),
            )
            .unwrap();
        assert_eq!(n_sites, 2);
        assert_eq!(n_motifs, 2);
        assert!((best_p - 3e-4).abs() < 1e-12);
        assert_eq!(motifs.split(',').collect::<BTreeSet<_>>(), ["MA1", "MA2"].into_iter().collect());
        assert_eq!(levels.split(',').collect::<BTreeSet<_>>(), ["A", "D"].into_iter().collect());
        drop(conn);

        std::fs::write(&targets, "gene_id\nG1\nG2\nG2\nG3\n").unwrap();
        let alias_error = cmd_prior(&[
            "--db".into(),
            db.to_string_lossy().into_owned(),
            "--targets".into(),
            targets.to_string_lossy().into_owned(),
            "--out".into(),
            db.to_string_lossy().into_owned(),
        ])
        .unwrap_err();
        assert!(alias_error.to_string().contains("different files"));
        let hardlink = dir.join("evidence-hardlink.db");
        std::fs::hard_link(&db, &hardlink).unwrap();
        let hardlink_error = cmd_prior(&[
            "--db".into(),
            db.to_string_lossy().into_owned(),
            "--targets".into(),
            targets.to_string_lossy().into_owned(),
            "--out".into(),
            hardlink.to_string_lossy().into_owned(),
        ])
        .unwrap_err();
        assert!(hardlink_error.to_string().contains("different files"));
        cmd_prior(&[
            "--db".into(),
            db.to_string_lossy().into_owned(),
            "--targets".into(),
            targets.to_string_lossy().into_owned(),
            "--out".into(),
            upstream.to_string_lossy().into_owned(),
        ])
        .unwrap();
        let text = std::fs::read_to_string(&upstream).unwrap();
        let rows: Vec<&str> = text.lines().collect();
        assert_eq!(rows.len(), 3);
        assert!(rows[1].starts_with("TF2\tTF2\t2\t3\t5.000000e-5\tA,D"));
        assert!(rows[2].starts_with("TF1\tTF1\t1\t2\t1.000000e-6\tD"));

        std::fs::write(&targets, "G1\n").unwrap();
        cmd_prior(&[
            "--db".into(),
            db.to_string_lossy().into_owned(),
            "--targets".into(),
            targets.to_string_lossy().into_owned(),
            "--out".into(),
            upstream.to_string_lossy().into_owned(),
        ])
        .unwrap();
        let text = std::fs::read_to_string(&upstream).unwrap();
        assert_eq!(text.lines().count(), 2);
        assert!(text.contains("TF1"));
        assert!(!text.contains("TF2"));

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn convert_jsonl_keeps_numeric_types_and_gzip_flag() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-jsonl-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let tsv = dir.join("evidence.tsv");
        let jsonl = dir.join("evidence.jsonl.gz");
        std::fs::write(
            &tsv,
            "tf_id\ttarget_id\tpromoter_start\tscore\tp_value\tmotif_id\n\
TF1\tG1\t101\t8.5\t1e-6\tMA1\n\
TF1\tG2\t\t\t\t \n",
        )
        .unwrap();
        cmd_convert(&[
            "--input".into(),
            tsv.to_string_lossy().into_owned(),
            "--out".into(),
            jsonl.to_string_lossy().into_owned(),
            "--format".into(),
            "jsonl".into(),
            "--gzip".into(),
        ])
        .unwrap();
        let decoder = flate2::read::GzDecoder::new(File::open(&jsonl).unwrap());
        let mut text = String::new();
        let mut reader = BufReader::new(decoder);
        std::io::Read::read_to_string(&mut reader, &mut text).unwrap();
        let lines: Vec<JsonValue> = text
            .lines()
            .map(|line| serde_json::from_str(line).unwrap())
            .collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0]["promoter_start"], JsonValue::Number(101.into()));
        assert_eq!(lines[0]["score"], JsonValue::Number(serde_json::Number::from_f64(8.5).unwrap()));
        assert_eq!(lines[1]["promoter_start"], JsonValue::Null);
        assert_eq!(lines[1]["p_value"], JsonValue::Null);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn failed_conversion_preserves_existing_output() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-failure-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let tsv = dir.join("evidence.tsv");
        let out = dir.join("evidence.db");
        std::fs::write(
            &tsv,
            "tf_id\ttarget_id\tscore\nTF1\tG1\tnot-a-number\n",
        )
        .unwrap();
        std::fs::write(&out, "keep-me").unwrap();
        let result = cmd_convert(&[
            "--input".into(),
            tsv.to_string_lossy().into_owned(),
            "--out".into(),
            out.to_string_lossy().into_owned(),
            "--format".into(),
            "sqlite".into(),
        ]);
        assert!(result.is_err());
        assert_eq!(std::fs::read_to_string(&out).unwrap(), "keep-me");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn tf_genes_gate_rejects_unknown_and_ambiguous() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-tf-gate-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let genes = dir.join("tf_genes.tsv");
        std::fs::write(
            &genes,
            "gene_id\ttf_symbol\tfamily\tidentification_status\n\
             AT1G01010\tNAC001\tNAC\tidentified\n\
             AT1G01020\tTEST\t\tambiguous\n",
        )
        .unwrap();
        let mut report = diagnostics::Report::new("scan");
        let map = vec![io::TfMotif {
            tf_id: "AT1G01010".into(),
            tf_symbol: "NAC001".into(),
            family: "NAC".into(),
            motif_id: "MA0001.1".into(),
            reference_tf_id: String::new(),
            mapping_method: String::new(),
            ..Default::default()
        }];
        validate_tf_identity(&map, Some(&genes), &mut report).unwrap();

        let unknown = vec![io::TfMotif {
            tf_id: "AT9G99999".into(),
            tf_symbol: "UNKNOWN".into(),
            family: String::new(),
            motif_id: "MA0001.1".into(),
            reference_tf_id: String::new(),
            mapping_method: String::new(),
            ..Default::default()
        }];
        let err = validate_tf_identity(&unknown, Some(&genes), &mut report).unwrap_err();
        assert!(err.to_string().contains("absent from --tf-genes"), "{err}");

        let ambiguous = vec![io::TfMotif {
            tf_id: "AT1G01020".into(),
            tf_symbol: "TEST".into(),
            family: String::new(),
            motif_id: "MA0001.1".into(),
            reference_tf_id: String::new(),
            mapping_method: String::new(),
            ..Default::default()
        }];
        let err = validate_tf_identity(&ambiguous, Some(&genes), &mut report).unwrap_err();
        assert!(err.to_string().contains("not marked identified"), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    fn mapping_row(tf_id: &str, motif_id: &str, domain_mode: &str) -> io::TfMotif {
        io::TfMotif {
            tf_id: tf_id.to_owned(),
            tf_symbol: tf_id.to_owned(),
            family: "bZIP".to_owned(),
            motif_id: motif_id.to_owned(),
            reference_tf_id: tf_id.to_owned(),
            mapping_method: "reciprocal_dbd_similarity".to_owned(),
            domain_mode: domain_mode.to_owned(),
            ..Default::default()
        }
    }

    #[test]
    fn gene_biotype_filter_keeps_only_requested_models() {
        let gene = |id: &str, biotype: Option<&str>| gff::GeneModel {
            gene_id: id.to_owned(),
            symbol: id.to_owned(),
            chrom: "Chr1".to_owned(),
            start: 1,
            end: 100,
            strand: '+',
            transcript_id: None,
            tss: 1,
            biotype: biotype.map(str::to_owned),
        };
        let genes = vec![
            gene("coding", Some("protein_coding")),
            gene("ncrna", Some("ncRNA")),
            gene("inferred", None),
        ];
        let filtered = filter_gene_biotype(genes.clone(), Some("protein_coding"));
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].gene_id, "coding");
        assert_eq!(filter_gene_biotype(genes, None).len(), 3);
    }

    #[test]
    fn mapping_filter_keeps_old_maps_without_labels() {
        let (map, stats) = filter_scan_mappings(
            vec![mapping_row("YEL009C", "MA0303.3", "")],
            false,
        );
        assert_eq!(map.len(), 1);
        assert_eq!(stats.total, 1);
        assert_eq!(stats.fallback, 0);
        assert_eq!(stats.retained, 1);
        assert_eq!(stats.excluded, 0);
    }

    #[test]
    fn mapping_filter_defaults_to_full_architecture_only() {
        let rows = vec![
            mapping_row("TF1", "MA1", "full_architecture"),
            mapping_row("TF2", "MA2", "single_domain_fallback"),
            mapping_row("TF3", "MA3", "external_alignment"),
        ];
        let (map, stats) = filter_scan_mappings(rows.clone(), false);
        assert_eq!(map.len(), 2);
        assert_eq!(stats.total, 3);
        assert_eq!(stats.fallback, 1);
        assert_eq!(stats.retained, 2);
        assert_eq!(stats.excluded, 1);

        let (all, stats) = filter_scan_mappings(rows, true);
        assert_eq!(all.len(), 3);
        assert_eq!(stats.retained, 3);
        assert_eq!(stats.excluded, 0);
        assert_eq!(stats.fallback, 1);
    }

    #[test]
    fn mapping_filter_reports_all_fallback_as_zero_retained() {
        let (map, stats) = filter_scan_mappings(
            vec![
                mapping_row("TF1", "MA1", "single_domain_fallback"),
                mapping_row("TF2", "MA2", " Single_Domain_Fallback "),
            ],
            false,
        );
        assert!(map.is_empty());
        assert_eq!(stats.total, 2);
        assert_eq!(stats.fallback, 2);
        assert_eq!(stats.retained, 0);
        assert_eq!(stats.excluded, 2);
    }

    #[test]
    fn mapping_filter_keeps_ties_and_multiple_motifs_per_gene() {
        let rows = vec![
            mapping_row("TF1", "MA1", "full_architecture"),
            mapping_row("TF1", "MA2", "full_architecture"),
            mapping_row("TF2", "MA1", "full_architecture"),
        ];
        let (map, stats) = filter_scan_mappings(rows, false);
        assert_eq!(map.len(), 3);
        assert_eq!(stats.excluded, 0);
    }

    #[test]
    fn all_fallback_map_controls_validate_and_scan() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-fallback-gate-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let motifs = dir.join("motifs.jaspar");
        std::fs::write(
            &motifs,
            ">MA0001.1 TEST
A  [ 10 0 0 0 ]
C  [ 0 10 0 0 ]
G  [ 0 0 10 0 ]
T  [ 0 0 0 10 ]
",
        )
        .unwrap();
        let promoters = dir.join("promoters.fa");
        std::fs::write(&promoters, ">G1
ACGTACGTACGT
").unwrap();
        let map = dir.join("map.tsv");
        std::fs::write(
            &map,
            "tf_id	tf_symbol	family	motif_id	reference_tf_id	mapping_method	domain_mode
             TF1	TF1	bZIP	MA0001.1	TF1	reciprocal_dbd_similarity	single_domain_fallback
",
        )
        .unwrap();
        let out = dir.join("out.tsv");
        let validate = |include: bool| {
            let mut args = vec![
                "--promoters".to_owned(),
                promoters.to_string_lossy().into_owned(),
                "--motifs".to_owned(),
                motifs.to_string_lossy().into_owned(),
                "--tf-map".to_owned(),
                map.to_string_lossy().into_owned(),
            ];
            if include {
                args.push("--include-single-domain-fallback".to_owned());
            }
            cmd_validate(&args)
        };
        let scan = |include: bool| {
            let mut args = vec![
                "--promoters".to_owned(),
                promoters.to_string_lossy().into_owned(),
                "--motifs".to_owned(),
                motifs.to_string_lossy().into_owned(),
                "--tf-map".to_owned(),
                map.to_string_lossy().into_owned(),
                "--out".to_owned(),
                out.to_string_lossy().into_owned(),
            ];
            if include {
                args.push("--include-single-domain-fallback".to_owned());
            }
            cmd_scan(&args)
        };
        assert!(validate(false).is_err(), "all-fallback validate must fail");
        assert!(scan(false).is_err(), "all-fallback scan must fail");
        assert!(validate(true).is_ok(), "included fallback validate must pass");
        assert!(scan(true).is_ok(), "included fallback scan must pass");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn converted_sqlite_passes_validate_db_and_prior() {
        let dir = std::env::temp_dir().join(format!(
            "regutrace-convert-db-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let tsv = dir.join("evidence.tsv");
        std::fs::write(
            &tsv,
            "tf_id	tf_symbol	target_id	target_symbol	evidence_type	evidence_level	method	source	pmid	condition	promoter_start	promoter_end	motif_id	hit_start	hit_end	strand	score	p_value	notes
             TF1	TF1	G1	G1	prediction_motif	D	jaspar_pwm	promoter_scan					MA1			+	1	1e-5	family=bZIP;motif_mapping={\"domain_mode\":\"full_architecture\"}
",
        )
        .unwrap();
        let db = dir.join("evidence.db");
        assert!(cmd_convert(&[
            "--input".to_owned(),
            tsv.to_string_lossy().into_owned(),
            "--out".to_owned(),
            db.to_string_lossy().into_owned(),
            "--format".to_owned(),
            "sqlite".to_owned(),
        ])
        .is_ok());
        assert!(cmd_validate_db(&[
            "--db".to_owned(),
            db.to_string_lossy().into_owned(),
        ])
        .is_ok());
        let targets = dir.join("targets.txt");
        std::fs::write(&targets, "G1
").unwrap();
        let prior = dir.join("prior.tsv");
        assert!(cmd_prior(&[
            "--db".to_owned(),
            db.to_string_lossy().into_owned(),
            "--targets".to_owned(),
            targets.to_string_lossy().into_owned(),
            "--out".to_owned(),
            prior.to_string_lossy().into_owned(),
        ])
        .is_ok());
        let text = std::fs::read_to_string(&prior).unwrap();
        assert!(text.contains("TF1"), "{text}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn background_must_be_positive_and_finite() {
        assert!(parse_background("0.25,0.25,0.25,0.25").is_ok());
        assert!(parse_background("2,1,1,0").is_err());
        assert!(parse_background("1,0,0,0").is_err());
        assert!(parse_background("NaN,0.25,0.25,0.25").is_err());
        assert!(parse_background("inf,0.25,0.25,0.25").is_err());
        let huge = parse_background("1e308,1e308,1e308,1e308").unwrap();
        for v in huge {
            assert!((v - 0.25).abs() < 1e-12, "v={v}");
        }
        let tiny = parse_background("1e-308,1e-308,1e-308,1e-308").unwrap();
        for v in tiny {
            assert!((v - 0.25).abs() < 1e-12, "v={v}");
        }
        assert!(parse_background("1e308,1e-308,1e-308,1e-308").is_err());
    }
}
