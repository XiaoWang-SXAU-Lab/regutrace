# Phase 3 report — motif scan architecture and statistical calibration

Date: 2026-11-25  
Repository: `regutrace` (uncommitted working tree based on `9ac4837` plus Phase
0–2.6 changes)  
Status: Phase 3 implementation complete; stopped before Phase 4

## Scope

Phase 3 separates TF–motif mapping, unique PWM scanning, motif–promoter sites,
promoter aggregation and TF–target evidence; adds exact scan reuse and a
persistent cache; adds an empirical promoter background and an attainable
p-value audit; validates the scanner independently; and runs five HF1 benchmarks.

It does **not** start Phase 2.7, tune family thresholds on the Phase 2.5
calibration set, promote exploratory mappings, change the default reference pack,
or rebuild an installer.

## Constraints preserved

- Strict RBH output is unchanged.
- JASPAR 2024 remains the default; JASPAR 2026 remains a candidate.
- `dbd_calibrated`, `dbd_ambiguous` and `dbd_exploratory` remain distinct from
  binding evidence.
- No non-reciprocal mapping was promoted.
- No commit, push, merge, branch switch or installer rebuild.
- CPU-only, offline, cross-platform code paths.

## Task status

| task | status | evidence |
|---|---|---|
| A. Scanner audit | done | [phase3-scanner-audit.md](phase3-scanner-audit.md) |
| B. Mapping / scan / site / aggregation separation | done for TSV/JSON layers; SQLite normalized projection deferred | [phase3-scan-data-model.md](phase3-scan-data-model.md) |
| C. Exact scan reuse and cache keys | done | `Pwm::digest`, `--reuse-exact-pwm`, `--scan-cache`, manifest |
| D. Empirical promoter background | done | `--background-mode empirical`, manifest raw counts |
| E. Attainable p-value audit and testability | done | `motif-testability.tsv`, statuses |
| F. Site/promoter/gene-set statistical layers | site p-values done; promoter descriptive; gene-set exploratory; q-values deferred | [phase3-statistical-methodology.md](phase3-statistical-methodology.md) |
| G. Independent statistical validation | done | `synthetic_validation.py`, 6 Python tests |
| H. Five HF1 runs | done | [phase3-hf1-comparison.md](phase3-hf1-comparison.md) |
| I. Motif enrichment interface | implemented and tested as exploratory | [phase3-enrichment-design.md](phase3-enrichment-design.md) |
| J. Deferred tasks preserved | done | [phase3-deferred-tasks.md](phase3-deferred-tasks.md) |

## What was implemented

### Exact PWM identity and reuse

`Pwm::digest()` is a SHA-256 over the normalized probability matrix bits and
matrix length. It excludes motif ID, name, family and source. Exact duplicates
share a digest; near-duplicates do not.

`--reuse-exact-pwm` scans one representative per digest and emits evidence for
every motif ID in the group. On the HF1 2024 and 2026 strict maps there were no
mapped exact-duplicate groups, so this did not reduce the current job count; it
fixes the correctness and future-pack behavior. The site layer nevertheless makes
the sharing explicit.

### Persistent cache

`--scan-cache <dir>` stores per-PWM hit lists keyed by SHA-256 over:

```text
algorithm version, scale, background bits, p-value, min-score,
strand mode, promoter-set hash, PWM digest
```

A changed threshold, background or promoter set produces a miss. The cache is
not keyed by motif accession, family or similarity. Cache hits and misses are
recorded in the manifest. The cache has no eviction or migration policy; that is
deferred.

### Site layer

`--sites-out` writes one row per unique PWM group and site:

```text
pwm_hash, motif_ids, gene_id, chrom, promoter_start/end/strand,
hit_start/end/strand, genomic_hit_start/end/strand,
score, p_value, tf_count, scan_key
```

For genome/GFF input, genomic coordinates are derived for both promoter
orientations. For `--promoters` FASTA input they are empty, because the FASTA
header is the only identifier. HF1 JASPAR 2026 had zero site rows missing
genomic coordinates.

### Testability

`--motif-testability-out` writes, for every library motif:

```text
motif_id, pwm_hash, matrix_length, information_content, max_score,
min_pvalue, requested_pvalue, attainable, status, mapped,
motif_id_hit_count, tf_count, background_mode
```

Statuses:

```text
significant_hits
no_significant_hits
threshold_unattainable
not_scanned
```

`invalid_motif` is not emitted: `read_jaspar` rejects malformed or duplicate-ID
motifs before scanning. The requested threshold is never silently lowered.
Motifs whose exact floor exceeds the threshold are skipped and reported.

### Empirical background

`--background-mode empirical` counts A/C/G/T in the supplied promoters, ignores
other bases, symmetrizes A/T and C/G, and adds one pseudocount per base. On HF1:

```text
A/T = 0.32555, C/G = 0.17445 (65.1% AT)
raw A,C,G,T = 12,146,728 / 6,536,535 / 6,517,360 / 12,214,377
```

The model is zero-order and is estimated from the same promoters that are
scanned. P-values are conditional on the estimate; no uncertainty correction is
applied. Uniform remains the default.

### Promoter summary and manifest

`--promoter-summary-out` writes per-promoter site count, best site p-value and
best PWM. The minimum p-value is descriptive, not a promoter-level p-value.

`--scan-manifest` records inputs, promoter hash, background provenance, thresholds,
PWM groups, duplicate groups, cache statistics, site/evidence counts, runtime and
outputs. It is the audit trail for a run.

### Enrichment

`regutrace enrich` implements a one-sided Fisher exact test on a binary
motif-hit matrix over the scanned promoter universe, with BH q-values across
motifs within each gene set. It is explicitly exploratory; it does not match
promoter length, GC, motif similarity or gene-set overlap. See
[phase3-enrichment-design.md](phase3-enrichment-design.md).

## Statistical calibration

### What is claimed

- Site p-values are exact for the discretized score and the supplied zero-order
  background for one pre-specified window.
- The minimum attainable p-value is exact for the same model.
- The empirical background removes a first-order base-composition bias in the
  synthetic GC test.

### What is not claimed

- No genome-wide site FDR.
- No promoter q-value.
- No calibrated regulatory enrichment.
- No dinucleotide/GC-window/chromatin null.
- No FIMO agreement (FIMO was not available).

### Independent validation highlights

- Python reimplementation vs CLI: maximum absolute p-value difference `0` to
  `3.3e-16` across 8 motif/background combinations; minimum p-values matched
  exactly.
- Fixed-window null calibration, uniform: fraction p≤0.05 was 0.035 (AT6) and
  0.036 (GC6); binomial tails 0.99.
- GC-rich promoters, GC motif: at 80% GC, uniform background gave 16.9%
  p≤0.05 (binomial tail 3.7e-43); empirical gave 4.7% (tail 0.69).
- Overlap dependence: lag-1 score correlation was -0.278 for ACGTAC and 0.840
  for AAAAAA; adjacent windows are not independent.
- Discretization: exhaustive enumeration of 4–8 bp motifs gave 0 difference
  between discretized and continuous tails for the tested matrices.
- Attainable floor: 4 bp/3.14 bits cannot reach 1e-4; 8 bp/6.29 bits can.

Full tables: [phase3-statistical-methodology.md](phase3-statistical-methodology.md).

## HF1 benchmark

Five runs, 37,415 genes, 1,000 bp promoters, p=1e-4, both strands, 8 threads.

| metric | R1 2024 legacy | R2 2024 reuse | R3 2026 uniform | R4 2026 empirical | R5 2026 + cal |
|---|---:|---:|---:|---:|---:|
| mapped TFs | 260 | 260 | 271 | 271 | 280 |
| PWM groups | — | 257 | 268 | 268 | 275 |
| unique sites | — | 3,068,556 | 3,369,011 | 2,798,756 | 3,674,864 |
| TF evidence rows | 3,090,137 | 3,090,137 | 3,390,592 | 2,818,188 | 3,757,212 |
| TF–target pairs | 1,785,471 | 1,785,471 | 1,960,433 | 1,793,597 | 2,133,069 |
| target union | 37,407 | 37,407 | 37,413 | 37,415 | 37,413 |
| testable motifs | 229 | 229 | 239 | 249 | 246 |
| unattainable motifs | 28 | 28 | 29 | 19 | 29 |
| cache hits/misses | — | 0/229 | 0/239 | 0/249 | 246/0 |
| scan wall (s) | 69.3 | 61.6 | 66.5 | 65.7 | 13.0 |

Full-row evidence equality including `notes`:

| comparison | left-only | right-only |
|---|---:|---:|
| Run 1 vs Run 2 | 0 | 0 |
| Run 3 vs pre-Phase-3 2026 DB | 0 | 0 |
| Run 5 strict subset vs Run 3 | 0 | 0 |

Peak working set: legacy 3,308 MiB; optimized 454 MiB (7.3× lower). The optimized
path streams edges per PWM group instead of accumulating and sorting all edges.

### Background effect

Run 3 → Run 4, same mapping:

- unique sites −16.9%;
- TF–target pairs −8.5%;
- testable motifs 239 → 249;
- unattainable motifs 29 → 19;
- target union 37,413 → 37,415.

This is a scientific parameter change, not an optimization. Uniform remains the
default.

### Calibrated DBD addition

Run 5 adds 9 HD-ZIP TFs (`dbd_calibrated`, not ambiguous/exploratory), 7 new PWM
groups, 305,853 sites and 172,636 TF–target pairs. The strict 271-TF subset is
byte-for-byte identical to Run 3. No Phase 2.5 config changed.

## Performance

| measurement | legacy | optimized |
|---|---:|---:|
| peak working set | 3,308 MiB | 454 MiB |
| cold scan wall | 69.3 s | 61.6–66.5 s |
| warm cache scan wall | — | 10.2 s (2026) / 13.0 s (Run 5) |
| evidence.tsv (2026) | 1,281,685,499 B | 1,281,685,499 B |
| evidence.db (2026) | 1,805,856,768 B | 1,805,856,768 B |
| site table (2026) | — | 813,787,596 B |

The evidence file size is identical because the evidence rows and notes are
unchanged. The site/testability/manifest outputs are additive.

## Tests

Rust:

```text
46 CLI unit tests
7 existing CLI integration tests
2 Phase 3 CLI integration tests
69 core unit tests
0 failed
```

Phase 3 tests cover:

- PWM digest equality/inequality;
- empirical background symmetry and normalization;
- plus/minus genomic coordinate conversion;
- cache roundtrip;
- exact reuse, site counting, cache invalidation, empirical background and
  unattainable status end-to-end;
- enrichment Fisher/BH and CLI integration.

Python:

```text
6 tests passed (scripts/phase3/test_phase3.py)
```

The Python tests cover probability normalization, null convolution against brute
force, tail monotonicity, binomial tails, GC frequencies and the KS diagnostic.

## Figures

- [scan-architecture.svg](figures/phase3/scan-architecture.svg)
- [background-effect.svg](figures/phase3/background-effect.svg)
- [attainable-floor.svg](figures/phase3/attainable-floor.svg)
- [pvalue-calibration.svg](figures/phase3/pvalue-calibration.svg)

## Artifacts

Code:

```text
crates/regutrace-core/src/pwm.rs          Pwm::digest
crates/regutrace-core/src/io.rs           EdgeStream::append
crates/regutrace-cli/src/scan_phase3.rs   Phase 3 scan path
crates/regutrace-cli/src/enrich.rs        exploratory enrichment
crates/regutrace-cli/src/main.rs          CLI options and dispatch
crates/regutrace-cli/src/pipeline.rs      pipeline integration
crates/regutrace-cli/tests/phase3.rs      integration tests
scripts/phase3/                           validation, benchmark, figures, tests
```

Local machine-readable outputs:

```text
benchmark-final/hf1-benchmark.tsv
benchmark-final/hf1-benchmark.json
benchmark-final/hf1-targets-per-tf.tsv
benchmark-final/hf1-testability-by-length.tsv
benchmark-final/hf1-run1-run2-equality.json
benchmark-final/hf1-run5-strict-subset-equality.json
hf1-run{1..5}/scan-manifest.json
hf1-run{2..5}/motif-testability.tsv or testability.tsv
hf1-run{2..5}/sites.tsv or motif-site-hits.tsv
synthetic/out/synthetic-validation.{tsv,json}
mem/{legacy,optimized}.json
```

The sync workspace receives a data-only bundle at
`analysis/regutrace-yeast-validation/phase3-scan-2026-11-25/`. It does not contain
`evidence.tsv`, `evidence.db`, the site tables or the cache; those remain local
because of size.

## Limitations

1. One genome, one promoter length, one p-value and one thread count.
2. Empirical background is zero-order and estimated from scanned promoters.
3. No site/promoter q-values; no dinucleotide or chromatin null.
4. Enrichment is exploratory only.
5. FIMO was not available; no comparison.
6. Peak memory was sampled externally at 250 ms.
7. No external DAP-seq/ChIP-seq benchmark; calibrated DBD remains level D.
8. Cache eviction and SQLite normalized schema are deferred.
9. `--promoters` input cannot produce genomic coordinates without GFF.

## Conclusion

Phase 3 delivers a behavior-preserving scan architecture: legacy and optimized
evidence are identical, memory is 7.3× lower, repeat runs are 5–6.5× faster with a
warm cache, unique sites are separated from TF evidence, and short motifs that
cannot reach the threshold are reported instead of appearing as biological
negatives. The empirical background is a material scientific option (16.9% fewer
sites), not a default. No mapping or evidence tier was promoted.

Per the Phase 3 instruction, work stops here. Phase 4 has not been started.
