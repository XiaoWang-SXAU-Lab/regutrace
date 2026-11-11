# ReguTrace

**Local, offline TF–target evidence backend.** ReguTrace builds candidate
TF–target relationships from genome/GFF/promoter sequences and PWM motifs, then
merges predictions with user-provided evidence (DAP-seq, ChIP-seq, EMSA,
reporter, curated lists) into one evidence ledger with explicit provenance.

ReguTrace is not a binding predictor and does not claim causality. Every row is
annotated with an evidence level and source; low-confidence predictions stay
separate from experimental evidence.

## Status

Phase 1: CLI + SQLite evidence store + synthetic/Biostrings benchmark.
The wgcna-dashboard toolkit (`/toolkits/regutrace`) bundles the CLI and JASPAR 2024 CORE motif libraries for seven taxa,
with local TF–target scans, free-text queries, and promoter binding-site
visualization.

## Evidence levels

| level | meaning | examples |
|---|---|---|
| A | point-to-point validation | EMSA, reporter, ChIP-qPCR, mutant |
| B | broad in vivo binding | ChIP-seq, ChIP-chip |
| C | broad in vitro binding | DAP-seq, PBM, SELEX |
| D | computational prediction | promoter PWM / motif |
| E | functional/associative evidence | perturbation, co-expression |

DAP-seq is in vitro; ChIP-seq is in vivo. Do not conflate them.

## Build

```bash
cargo build --release
```

Output: `target/release/regutrace` (`.exe` on Windows).

## Quick start

Fetch a minimal Arabidopsis data pack (genome/GFF/protein + JASPAR plant PFMs):

```bash
bash scripts/fetch-data.sh
```

Scan promoters (1 kb upstream) with an optional TF→motif map:

```bash
regutrace scan \
  --genome data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
  --gff data/Arabidopsis_thaliana.TAIR10.60.gff3.gz \
  --motifs data/JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt \
  --tf-map tf_map.tsv \
  --promoter-len 1000 --pvalue 1e-4 \
  --out evidence.tsv
```

Store and query evidence:

```bash
regutrace import --db reguledger.db --evidence evidence.tsv
regutrace stats  --db reguledger.db
regutrace export --db reguledger.db --out evidence_export.tsv
```

`--tf-map` is optional. When it is omitted, `tf_id`/`tf_symbol` are taken from
the motif header name (for example JASPAR `>MA1758.2 HSFA1E`); provide a tf-map
when you need gene IDs or family annotations.

Add `--downstream 100` to include 100 bp past the TSS; PlantRegMap uses
TSS −500/+100, i.e. `--promoter-len 500 --downstream 100`.

Or scan a promoter FASTA directly (useful for benchmarks/other tools):

```bash
regutrace scan --promoters promoters.fa --motifs ... --tf-map ... --out evidence.tsv
```

## Input validation

Review genome/GFF/motif/TF-map inputs without scanning and get a
machine-readable report:

```bash
regutrace validate \
  --genome data/genome.fa --gff data/annotation.gff3 \
  --motifs data/jaspar.txt --tf-map tf_map.tsv \
  --report-json report.json
```

`scan` also accepts `--report-json`. Reports use stable codes, severities,
context and suggestions for frontends and RAG assistants; see
[docs/DIAGNOSTICS.md](docs/DIAGNOSTICS.md).

## Large scans, streaming and progress

A whole-genome scan with hundreds of motifs can emit millions of rows. Add
`--stream` to write the evidence TSV incrementally instead of retaining every
row in memory (rows are in scan order, not sorted), and `--progress-json` to
publish stage/progress for a frontend. Scanning is parallel by default, using
`CPU cores - 1` worker threads; set `--threads N` explicitly or `--threads 1`
for serial execution:

```bash
regutrace scan \
  --genome genome.fa --gff annotation.gff3 --motifs jaspar_plants.txt \
  --promoter-len 1000 --pvalue 1e-4 \
  --out evidence.tsv --stream --threads 8 \
  --progress-json progress.json --report-json report.json
```

`progress.json` is updated atomically with fields `stage`, `status`, `done`,
`total`, `progress` and `message`. Stages include `reading_motifs`,
`reading_genome`, `reading_gff`, `extracting_promoters`, `scanning`,
`writing_evidence`, `done` and `failed`. The dashboard toolkit consumes this
file, indexes the evidence TSV on disk by TF/target/motif, and keeps queries
interactive without loading millions of rows into memory.

On a Tartary buckwheat benchmark (2,000 promoters × 805 JASPAR motifs), 7
worker threads reduced the scan from 59 s to 17.2 s (about 3.4×) with identical
528,201 hits. The default is `available_parallelism() - 1`; use `--threads 1`
for the old serial behaviour.

## Evidence schema

Unified TSV columns:

```text
tf_id  tf_symbol  target_id  target_symbol  evidence_type  evidence_level
method  source  pmid  condition  promoter_start  promoter_end  motif_id
hit_start  hit_end  strand  score  p_value  notes
```

Coordinates are 1-based inclusive. `hit_start/hit_end` are relative to the
promoter oriented 5'→3' in the gene's direction. See
[docs/EVIDENCE.md](docs/EVIDENCE.md).

## Benchmark

`scripts/benchmark.R` compares ReguTrace with Biostrings::matchPWM at a matched
false-positive rate. Current synthetic result (1,200 promoters × HSFA1E motif):
comparable accuracy (AP 0.970 vs 0.970) and ~80× faster (0.05 s vs 4.19 s).
See [docs/BENCHMARK.md](docs/BENCHMARK.md).

## Repository layout

```text
crates/regutrace-core/   FASTA/GFF/promoter/PWM/scan/evidence/SQLite
crates/regutrace-cli/    `regutrace` CLI
scripts/                 data-pack fetch + benchmark
docs/                    plan, evidence schema, benchmark, API plan
```

## Representative transcript

Each gene uses one representative transcript: the longest spliced transcript
(sum of exon lengths; CDS + UTR lengths when exon features are absent; genomic
span as a last resort). Ties prefer the most upstream TSS; a gene without any
transcript feature falls back to its gene coordinates.

## Limitations

- PWM prediction has high false-positive rates; many plant TF families share
  degenerate motifs. Predictions are family-level candidates.
- Promoter definition (length, TSS, UTR, transcript choice) strongly affects
  results; coordinates are always reported.
- Without ATAC/ChIP, prediction cannot distinguish accessible from bound sites.
- A `tf-map` is optional; without it, TF IDs come from motif header names.
  Mapping those names to gene IDs and Pfam HMM TF classification remain Phase 2.
- GFF support covers `gene` + `mRNA`/`transcript`, `mRNA` without an explicit
  gene, gene-only models, and CDS/exon-only fallbacks. Unusual feature layouts
  are reported by `validate` rather than failing silently.
- MIT licensed. Reference data licenses are documented by `fetch-data.sh`.
