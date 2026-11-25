# Phase 3 scan architecture and statistics

Run from the repository root. Paths below are examples; replace `RX_HOME` with
the local build/data root.

## Build

```bash
cargo build --release -p regutrace-cli
```

## Optimized scan

```bash
regutrace scan \
  --genome HF1.fa --gff HF1.gff \
  --motifs JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt \
  --tf-map tf_map.tsv --tf-genes tf_genes.tsv \
  --promoter-len 1000 --pvalue 1e-4 --threads 8 \
  --reuse-exact-pwm --scan-cache "$RX_HOME/cache" \
  --sites-out sites.tsv \
  --motif-testability-out testability.tsv \
  --promoter-summary-out promoter-summary.tsv \
  --scan-manifest scan-manifest.json \
  --out evidence.tsv
```

`--sites-out`, `--motif-testability-out`, `--promoter-summary-out`,
`--scan-manifest`, `--scan-cache`, `--reuse-exact-pwm` or a non-uniform
`--background-mode` switches to the Phase 3 path. Without them the legacy path is
used unchanged.

## Empirical background

```bash
regutrace scan ... --background-mode empirical ...
```

Do not pass `--background` with `--background-mode empirical`; the values are
estimated from the supplied promoters and recorded in `scan-manifest.json`.

## Pipeline

```bash
regutrace pipeline \
  --reference-pack pack --genome HF1.fa --gff HF1.gff --proteins HF1.pep \
  --out-dir run --threads 8 \
  --reuse-exact-pwm --phase3-artifacts --scan-cache "$RX_HOME/cache"
```

`--phase3-artifacts` writes `motif-site-hits.tsv`, `motif-testability.tsv`,
`promoter-summary.tsv` and `scan-manifest.json` alongside the usual outputs.
`--background-mode empirical` is also supported.

## HF1 benchmark

```bash
python scripts/phase3/benchmark_hf1.py \
  --run run1=.../hf1-run1-2024-legacy \
  --run run2=.../hf1-run2-2024-optimized-v2 \
  --run run3=.../hf1-run3-2026-uniform-v2 \
  --run run4=.../hf1-run4-2026-empirical-v2 \
  --run run5=.../hf1-run5-2026-calibrated-v3 \
  --run1-run2-equality run1 run2 \
  --run5-strict-subset run5 run3 run3/tf_map.tsv \
  --out benchmark
```

The script reads the evidence databases and, when present, the site,
testability and manifest tables. It does not rescan.

## Peak memory

```bash
python scripts/phase3/measure_peak_memory.py \
  --label scan --out memory.json --log scan.log -- \
  regutrace scan ...
```

On Windows the wrapper samples `tasklist`; on POSIX it reads `VmHWM`.

## Independent synthetic validation

```bash
python scripts/phase3/synthetic_validation.py \
  --binary ./target/release/regutrace \
  --workdir synthetic/work --out synthetic/out
python scripts/phase3/test_phase3.py
```

The validation reimplements the scanner in Python and compares p-values to the
CLI on fixed windows. It is independent of the Rust implementation.

## Figures

```bash
python scripts/phase3/make_figures.py \
  --benchmark benchmark/hf1-benchmark.tsv \
  --legacy-memory memory/legacy.json \
  --optimized-memory memory/optimized.json \
  --testability run3/motif-testability.tsv \
  --synthetic synthetic/out/synthetic-validation.json \
  --out docs/development/figures/phase3
```

## Enrichment

```bash
regutrace enrich \
  --sites run/motif-site-hits.tsv \
  --promoter-summary run/promoter-summary.tsv \
  --gene-sets sets.gmt \
  --out enrichment.tsv --summary-out enrichment-summary.tsv
```

This is exploratory Fisher/BH on a binary promoter-hit matrix. It is not a
calibrated regulatory enrichment model. See
`docs/development/phase3-enrichment-design.md`.

## Tests

```bash
cargo test --workspace
python scripts/phase3/test_phase3.py
```
