# Benchmark

## Design

A single relative threshold (e.g. 80% of max score) is not comparable across
implementations because PWM score scaling and background/pseudocount choices
differ. `scripts/benchmark.R` instead calibrates each tool on negative promoters
to a matched false-positive rate (target FPR = 1%), then compares precision,
recall, F1 and average precision (AP).

Synthetic data:

- 1,200 promoters × 1 kb
- 600 positive promoters: 300 consensus sites + 300 weak sites (50–85% of max
  PWM score), plus/minus strand
- 600 negative promoters with no planted site
- motif: HSFA1E `MA1758.2`
- both strands scanned

Reference tool: Biostrings::matchPWM (Bioconductor).

## Result (2026-11-11)

| tool | threshold | TP | FP | FN | precision | recall | F1 | AP | time |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| regutrace | 21.08 (raw log2 odds) | 300 | 4 | 300 | 0.987 | 0.5 | 0.664 | 0.970 | **0.05 s** |
| Biostrings::matchPWM | 0.975 (unit scale) | 300 | 2 | 300 | 0.993 | 0.5 | 0.665 | 0.970 | 4.19 s |

- Both tools detected all 300 consensus sites; both missed the 300 weak sites
  at 1% FPR (threshold sensitivity, not a tool difference).
- Accuracy is comparable at matched FPR (AP 0.970 vs 0.970); regutrace is
  ~80× faster end-to-end (1,200 promoters × 1 motif, including process startup
  and file parsing).
- Real 1,000 Arabidopsis promoters: 4 vs 2 hits, Jaccard 0.5; 0.05 s vs 2.06 s.
  Both tools are very stringent at 1% FPR.

## Reproduce

```bash
Rscript scripts/benchmark.R target/release/regutrace.exe \
  data/JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt \
  MA1758.2 benchmark_out 1200 1000 \
  benchmark_out/real_promoters.fa
```

## Limitations

- Biostrings timing includes R loop overhead; a vectorized/compiled reference
  may be faster.
- No FIMO/HOMER comparison yet; add container-based runs in Phase 2.
- No DAP-seq/ChIP-seq biological benchmark yet; this is the next accuracy check.
- Synthetic consensus/weak sites test sensitivity and FPR, not biological
  relevance.
