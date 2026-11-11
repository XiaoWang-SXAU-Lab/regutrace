//! regutrace: local TF–Target candidate builder (Phase 1).
//!
//! Example:
//!   regutrace scan --genome genome.fa --gff annotation.gff3 \
//!     --motifs jaspar_plants.txt --tf-map tf_map.tsv \
//!     --promoter-len 1000 --pvalue 1e-4 --out evidence.tsv

use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use regutrace_core::{diagnostics, fasta, gff, io, promoter, pwm, scan, store::Store};

struct Args {
    genome: Option<PathBuf>,
    gff: Option<PathBuf>,
    motifs: Option<PathBuf>,
    tf_map: Option<PathBuf>,
    out: Option<PathBuf>,
    annotation: Option<PathBuf>,
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
            out: None,
            annotation: None,
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
    [--tf-map <tf_map.tsv>] --out <edges.tsv> [options]\n\
  regutrace import --db <evidence.db> --evidence <edges.tsv>\n\
  regutrace export --db <evidence.db> --out <edges.tsv>\n\
  regutrace stats --db <evidence.db>\n\
  regutrace validate [--genome <fa> --gff <gff3>] [--motifs <jaspar>] \\\n\
    [--tf-map <tsv>] [--annotation <tsv>] [--user-list <tsv>] \\\n\
    [--report-json <report.json>] [--strict]\n\
\n\
SCAN OPTIONS:\n\
  --promoter-len <bp>     upstream length (default 1000)\n\
  --downstream <bp>       downstream window incl. TSS (default 0; PlantRegMap uses 100)\n\
  --tf-map <tsv>          TF–motif map (optional; derived from motif names if omitted)\n\
  --promoters <fa>        scan a promoter FASTA directly (skip genome/GFF)\n\
  --promoters-out <fa>    export extracted promoters as FASTA\n\
  --pvalue <p>            PWM hit p-value cutoff (default 1e-4)\n\
  --min-score <s>         minimum PWM score (log2 odds; for benchmarks)\n\
  --background <A,C,G,T>  base probabilities (default 0.25,0.25,0.25,0.25)\n\
  --annotation <tsv>      gene_id/symbol annotation for target symbols\n\
  --user-list <tsv>       user TF–target evidence list to append\n\
  --max-genes <n>         only process the first n genes (testing)\n\
  --one-strand            scan only the plus strand of each promoter\n\
  --both-strands          scan both strands (default)\n\
  --report-json <path>    write a machine-readable diagnostics report\n\
  --progress-json <path>  write live stage/progress JSON for frontends\n\
  --stream                write evidence rows incrementally (unsorted; large scans)\n\
  --threads <n>           scan worker threads (default CPU cores - 1; 1 = serial)\n\
  --strict                treat warnings as failure (validate)\n\
  -h, --help              show this help\n\
\n\
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
            "--out" => a.out = Some(PathBuf::from(it.next().context("--out needs a value")?)),
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
                        notes: format!("family={}", tf.family),
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
    let tf_map: Vec<io::TfMotif> = match &args.tf_map {
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
                    }
                })
                .collect()
        }
    };
    if let Some(path) = &args.tf_map {
        report.add_input("tf_map", path);
        diagnostics::review_tf_map(&tf_map, &motifs, path, report);
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
    if motif_to_tfs.is_empty() {
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
    let mut out = HashMap::new();
    let mut it = args.iter();
    while let Some(key) = it.next() {
        let value = it
            .next()
            .with_context(|| format!("{key} needs a value"))?;
        out.insert(key.clone(), value.clone());
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

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let res = match args.first().map(String::as_str) {
        None | Some("-h") | Some("--help") => {
            usage();
            return;
        }
        Some("scan") => cmd_scan(&args[1..]),
        Some("validate") => cmd_validate(&args[1..]),
        Some("import") => cmd_import(&args[1..]),
        Some("export") => cmd_export(&args[1..]),
        Some("stats") => cmd_stats(&args[1..]),
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
