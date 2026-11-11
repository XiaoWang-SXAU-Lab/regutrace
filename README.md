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
HTTP API and dashboard integration are Phase 3.

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

Scan promoters (1 kb upstream) with a TF→motif map:

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

Or scan a promoter FASTA directly (useful for benchmarks/other tools):

```bash
regutrace scan --promoters promoters.fa --motifs ... --tf-map ... --out evidence.tsv
```

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

## Limitations

- PWM prediction has high false-positive rates; many plant TF families share
  degenerate motifs. Predictions are family-level candidates.
- Promoter definition (length, TSS, UTR, transcript choice) strongly affects
  results; coordinates are always reported.
- Without ATAC/ChIP, prediction cannot distinguish accessible from bound sites.
- The current `tf-map` is user-provided; automatic JASPAR-name/annotation
  mapping and Pfam HMM TF classification are Phase 2.
- MIT licensed. Reference data licenses are documented by `fetch-data.sh`.
