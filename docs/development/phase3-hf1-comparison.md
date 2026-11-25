# Phase 3 HF1 comparison

Date: 2026-11-25  
Inputs: HF1 T2T genome/GFF (`a.Fagopyrum_tataricum.t2t.fa`,
`d.Fagopyrum_tataricum.gene.as.gff`), 37,415 genes, 1,000 bp promoters,
p = 1e-4, both strands, 8 threads  
Machine-readable tables:
`analysis/regutrace-yeast-validation/phase3-scan-2026-11-25/benchmark/`

## Runs

| run | pack | mapping | background | architecture |
|---|---|---|---|---|
| 1 | JASPAR 2024 | strict RBH, 260 TFs | uniform | legacy in-memory sorted |
| 2 | JASPAR 2024 | same map as Run 1 | uniform | exact PWM reuse, streamed |
| 3 | JASPAR 2026 | strict RBH, 271 TFs | uniform | exact reuse, streamed |
| 4 | JASPAR 2026 | strict RBH, 271 TFs | empirical promoter | exact reuse, streamed |
| 5 | JASPAR 2026 | strict 271 + 9 calibrated HD-ZIP | uniform | exact reuse, warm cache |

Run 5 uses the provenance-preserving combined map
`hf1-strict-plus-calibrated-with-provenance.tsv`: the 271 rows from the Phase 3
strict map plus the 9 calibrated HD-ZIP rows that are not already strict. The
Phase 2.5 calibrated configs were not modified.

## Main counts

| metric | Run 1 | Run 2 | Run 3 | Run 4 | Run 5 |
|---|---:|---:|---:|---:|---:|
| identified TFs | 1,907 | 1,907 | 1,907 | 1,907 | 1,907 |
| mapped TFs | 260 | 260 | 271 | 271 | 280 |
| TF–motif rows | 260 | 260 | 271 | 271 | 280 |
| mapped motif IDs | 257 | 257 | 268 | 268 | 275 |
| PWM groups | — | 257 | 268 | 268 | 275 |
| unique PWMs | — | 257 | 268 | 268 | 275 |
| duplicate PWM groups | — | 0 | 0 | 0 | 0 |
| unique PWM sites | — | 3,068,556 | 3,369,011 | 2,798,756 | 3,674,864 |
| motif-ID site hits | 3,068,556 | 3,068,556 | 3,369,011 | 2,798,756 | 3,674,864 |
| TF evidence rows | 3,090,137 | 3,090,137 | 3,390,592 | 2,818,188 | 3,757,212 |
| motif–target pairs | 1,767,857 | 1,767,857 | 1,942,819 | 1,777,195 | 2,081,982 |
| TF–target pairs | 1,785,471 | 1,785,471 | 1,960,433 | 1,793,597 | 2,133,069 |
| target union | 37,407 | 37,407 | 37,413 | 37,415 | 37,413 |
| testable mapped motifs | 229 | 229 | 239 | 249 | 246 |
| unattainable motifs | 28 | 28 | 29 | 19 | 29 |

`motif-ID site hits` counts each motif ID once per hit. If two motif IDs shared
an exact matrix, the unique site count would be smaller; no mapped exact
duplicate groups occurred in these runs, so the two columns are equal.

The 2024 pack has 257 scanned motif IDs but 260 TF–motif rows: three TFs share a
motif with another TF. This is why 3,068,556 site hits become 3,090,137 evidence
rows (21,581 extra rows). The 2026 pack has 268 scanned motif IDs and 271 rows
(3 shared rows), giving 3,369,011 sites and 3,390,592 evidence rows.

### Targets per mapped TF

| run | min | median | mean | max |
|---|---:|---:|---:|---:|
| 1 | 146 | 5,327.5 | 7,696.0 | 30,694 |
| 2 | 146 | 5,327.5 | 7,696.0 | 30,694 |
| 3 | 146 | 5,618 | 8,101.0 | 30,694 |
| 4 | 1,886 | 6,666.5 | 7,117.4 | 17,189 |
| 5 | 146 | 5,847 | 8,498.3 | 30,694 |

The empirical background raises the minimum (weaker TFs lose their only sites)
and lowers the maximum (the strongest TFs lose many redundant sites), while the
target union changes by only two genes.

### Output sizes

| run | evidence.tsv | evidence.db | site table |
|---|---:|---:|---:|
| 1 | 1,166,668,577 B | 1,645,625,344 B | — |
| 2 | 1,166,668,577 B | 1,645,617,152 B | 741,221,483 B |
| 3 | 1,281,685,499 B | 1,805,856,768 B | 813,787,596 B |
| 4 | 1,063,343,432 B | 1,523,843,072 B | 676,384,454 B |
| 5 | 1,424,729,883 B | 1,995,755,520 B | 887,695,921 B |

DB sizes differ slightly between apparently identical evidence runs because the
insert order differs (legacy sorted vs optimized group order), so SQLite packs
pages differently; row counts and full-row sets are identical. The evidence TSV
byte sizes for Runs 1/2 and the pre-Phase-3 2026 run are identical because the
rows and notes are identical.

## Exactness checks

Full-row set comparisons include all 19 evidence columns, including `notes`.

| comparison | left-only | right-only | result |
|---|---:|---:|---|
| Run 1 vs Run 2 | 0 | 0 | identical |
| Run 3 optimized vs pre-Phase-3 JASPAR 2026 DB | 0 | 0 | identical |
| Run 5 strict subset vs Run 3 | 0 | 0 | identical |

The Run 3 `evidence.tsv` is 1,281,685,499 bytes, exactly the size of the
pre-Phase-3 HF1 2026 evidence table. The Phase 3 scan does not change mapping,
evidence rows or notes.

## Runtime and memory

| measurement | value |
|---|---:|
| legacy full scan, wall | 69.3 s |
| optimized cold scan, wall | 61.6 s (2024) / 66.5 s (2026) |
| optimized warm cache scan, wall | 10.2 s (2026) |
| Run 5 warm cache scan | 13.0 s |
| legacy peak working set | 3,308 MiB |
| optimized peak working set | 454 MiB |

The optimized path is 7.3× smaller in peak working set and slightly faster cold.
The large speed gain appears when the cache is warm. Memory is lower because
edges are streamed per PWM group instead of being accumulated and sorted in RAM.

Pipeline wall time (including TF identification, mapping, scan, SQLite
conversion, trace and report) was 151 s for Run 1, 152 s for Run 3 and 142 s for
Run 4. Run 2 and Run 5 were direct scans using existing maps, so their wall times
(63 s and 14 s) are not pipeline totals.

## Cache behavior

| run | cache hits | cache misses | note |
|---|---:|---:|---|
| Run 2 | 0 | 229 | cold 2024 cache |
| Run 3 | 0 | 239 | cold 2026 uniform cache |
| Run 4 | 0 | 249 | cold 2026 empirical cache |
| Run 5 | 246 | 0 | warm, reuses Run 3 cache |

A cold Run 5 with the same cache had 239 hits and 7 misses; the 7 misses were the
new HD-ZIP matrices. The warm Run 5 then had 246 hits and 0 misses. The cache key
includes the p-value, background, promoter set, strand mode and PWM hash; changing
the p-value or background produced a miss in the integration test.

## Background effect (Run 3 vs Run 4)

Same mapping, same 271 TFs, only the background changed.

| quantity | uniform | empirical | change |
|---|---:|---:|---:|
| unique sites | 3,369,011 | 2,798,756 | −16.9% |
| TF evidence rows | 3,390,592 | 2,818,188 | −16.9% |
| motif–target pairs | 1,942,819 | 1,777,195 | −8.5% |
| TF–target pairs | 1,960,433 | 1,793,597 | −8.5% |
| target union | 37,413 | 37,415 | +2 |
| testable motifs | 239 | 249 | +10 |
| unattainable motifs | 29 | 19 | −10 |
| median targets per TF | 5,618 | 6,667 | +18.7% |

The empirical background is 65.1% AT (A/T 0.32555, C/G 0.17445). It reduces
sites for GC-rich motifs and increases testability for short AT-rich motifs. It
does not change the mapping. The target-level effect is smaller than the site-
level effect because many TFs retain at least one site per target.

This is not evidence that the empirical model is "more correct" in a biological
sense. The synthetic GC calibration shows that it removes the first-order
base-composition bias of a uniform background; it remains a zero-order model.

## Calibrated DBD addition (Run 5 vs Run 3)

Run 5 adds 9 HD-ZIP TFs that are not strict 2026 mappings. All 9 are
`dbd_calibrated`; no `dbd_ambiguous` or `dbd_exploratory` row was added.

| quantity | Run 3 strict | Run 5 + calibrated | difference |
|---|---:|---:|---:|
| mapped TFs | 271 | 280 | +9 |
| PWM groups | 268 | 275 | +7 |
| unique sites | 3,369,011 | 3,674,864 | +305,853 |
| TF evidence rows | 3,390,592 | 3,757,212 | +366,620 |
| TF–target pairs | 1,960,433 | 2,133,069 | +172,636 |

The 9 new TFs have 14,203–25,595 targets each; their target counts sum to
172,636, exactly the Run 5 − Run 3 TF–target difference. The strict 271-TF
subset is byte-for-byte identical to Run 3, so the calibrated addition does not
disturb strict evidence.

The 9 TFs remain evidence level D. `dbd_calibrated` is a mapping confidence, not
experimental binding evidence. No Phase 2.5 threshold was changed.

## Interpretation

1. **Architecture is behavior-preserving.** Runs 1/2 and Run 3/original have
   identical full-row evidence sets and identical file sizes. The new site,
   testability and manifest outputs are additive.
2. **Exact reuse is real but rare in these packs.** No mapped exact-duplicate
   PWM groups occurred in the 2024 or 2026 strict maps. The value of the grouping
   layer is correctness and future packs, not a large current speedup. The site
   layer and streaming are the main operational gains.
3. **Caching is the main repeat-run speedup.** Warm-cache scans are 5–6.5× faster
   than cold scans while producing the same evidence.
4. **Background choice matters more than scan architecture.** The empirical
   model changes 16.9% of sites; the architecture changes none. This is a
   scientific parameter, not an optimization.
5. **Testability is not optional.** 29/268 (10.8%) of mapped 2026 motifs cannot
   reach p=1e-4 under uniform; without the table they look like motifs with no
   biology. The empirical background reduces this to 19/268 (7.1%).
6. **Calibrated DBD coverage is still small and additive.** Nine HD-ZIP TFs add
   366,620 evidence rows but remain level D. This is consistent with Phase 2.5's
   conservative calibration.

## Reproducibility

Commands and full logs are in `scripts/phase3/README.md` and the bundle
`benchmark/`, `manifests/`, `testability/`, `synthetic/` directories. No installer
was rebuilt, no commit or push was performed, and no Phase 2.5 config was
overwritten.

## Limitations

- One genome, one promoter length (1,000 bp), one p-value (1e-4), one thread
  count (8). Runtime scales with these choices.
- Peak working set was sampled externally at 250 ms; short transients may be
  missed.
- The empirical background is zero-order and estimated from the scanned
  promoters.
- No site/promoter FDR is reported.
- No external DAP-seq/ChIP-seq benchmark was used.
- The 9 calibrated HD-ZIP TFs are a mapping-coverage addition, not a validated
  target-level improvement.
