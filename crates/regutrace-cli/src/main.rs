//! regutrace: local TF–Target candidate builder (Phase 1).
//!
//! Example:
//!   regutrace scan --genome genome.fa --gff annotation.gff3 \
//!     --motifs jaspar_plants.txt --tf-map tf_map.tsv \
//!     --promoter-len 1000 --pvalue 1e-4 --out evidence.tsv

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use regutrace_core::{fasta, gff, io, promoter, pwm, scan, store::Store};

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
        }
    }
}

fn usage() {
    eprintln!(
        "regutrace (Phase 1) — local TF–Target candidate builder\n\
\n\
USAGE:\n\
  regutrace scan --genome <fa[.gz]> --gff <gff3[.gz]> --motifs <jaspar.txt> \\\n\
    --tf-map <tf_map.tsv> --out <edges.tsv> [options]\n\
  regutrace import --db <evidence.db> --evidence <edges.tsv>\n\
  regutrace export --db <evidence.db> --out <edges.tsv>\n\
  regutrace stats --db <evidence.db>\n\
\n\
SCAN OPTIONS:\n\
  --promoter-len <bp>     upstream length (default 1000)\n\
  --downstream <bp>       downstream window incl. TSS (default 0; PlantRegMap uses 100)\n\
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
  -h, --help              show this help\n\
\n\
tf-map columns: tf_id, tf_symbol, family, motif_id (one row per TF–motif pair)\n\
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
    Ok(a)
}

fn required<T>(v: Option<T>, name: &str) -> Result<T> {
    v.with_context(|| format!("missing required argument {name}"))
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

fn cmd_scan(args: &[String]) -> Result<()> {
    let args = parse_scan_args(args)?;
    let motifs_path = required(args.motifs.clone(), "--motifs")?;
    let tf_map_path = required(args.tf_map.clone(), "--tf-map")?;
    let out_path = required(args.out.clone(), "--out")?;

    eprintln!("reading motifs: {}", motifs_path.display());
    let motifs = pwm::read_jaspar(&motifs_path)?;
    let mut motif_by_id: HashMap<String, pwm::Pwm> = HashMap::new();
    for m in motifs {
        motif_by_id.entry(m.id.clone()).or_insert(m);
    }
    let tf_map = io::read_tf_map(&tf_map_path)?;
    let annotation = match &args.annotation {
        Some(p) => io::read_annotation(p)?,
        None => HashMap::new(),
    };

    let (promoters, gene_count, missing_chrom) = if let Some(promoters_path) = &args.promoters {
        eprintln!("reading promoter FASTA: {}", promoters_path.display());
        let fa = fasta::read_fasta(promoters_path)?;
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
                seq: fa.seqs[id].clone(),
            });
        }
        (v, fa.order.len(), 0usize)
    } else {
        let genome_path = required(args.genome.clone(), "--genome")?;
        let gff_path = required(args.gff.clone(), "--gff")?;
        eprintln!("reading genome: {}", genome_path.display());
        let genome = fasta::read_fasta(&genome_path)?;
        eprintln!("reading GFF: {}", gff_path.display());
        let genes = gff::read_gff(&gff_path)?;
        let mut v = Vec::new();
        let mut missing = 0usize;
        for gene in &genes {
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
        }
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
        write_promoters_fasta(path, &promoters)?;
        eprintln!("wrote promoter FASTA: {}", path.display());
    }

    // Motif -> TF mappings: scan each motif once per promoter, then emit all TFs
    // that use it.
    let mut motif_to_tfs: HashMap<String, Vec<&io::TfMotif>> = HashMap::new();
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

    let mut edges: Vec<io::EdgeRecord> = Vec::new();
    let mut prediction_hits = 0usize;
    for prom in &promoters {
        for (motif_id, tfs) in &motif_to_tfs {
            let scanner = scanners.get(motif_id).unwrap();
            let hits = scanner.scan(&prom.seq, args.pvalue, args.min_score, args.both_strands);
            for hit in hits {
                prediction_hits += 1;
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
    }

    let user_count = match &args.user_list {
        Some(path) => {
            let mut user = io::read_user_list(path)?;
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
            edges.extend(user);
            n
        }
        None => 0,
    };

    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("create output directory {}", parent.display()))?;
        }
    }
    io::write_edges(&out_path, &mut edges)?;

    eprintln!(
        "motifs: {}, TF–motif pairs: {}, promoters: {}, prediction hits: {}, user edges: {}, total edges: {}",
        motif_to_tfs.len(),
        tf_map.len() - missing_motif,
        promoters.len(),
        prediction_hits,
        user_count,
        edges.len()
    );
    eprintln!("wrote {}", out_path.display());
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
