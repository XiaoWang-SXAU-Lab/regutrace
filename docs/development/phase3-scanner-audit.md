# Phase 3 scanner audit

Date: 2026-11-25  
Code base: `regutrace` at working tree `9ac4837` plus Phase 0–2.6 changes  
Scope: audit of the existing promoter scanner before Phase 3 changes

## Existing scan path

`cmd_scan` (`crates/regutrace-cli/src/main.rs`) builds:

1. a motif library (`pwm::read_jaspar`);
2. a TF–motif map (`io::read_tf_map`);
3. a `motif_to_tfs: BTreeMap<String, Vec<&TfMotif>>` index;
4. one `scan::Scanner` per motif ID;
5. promoter chunks scanned in parallel.

`scan_chunk` iterates promoters, then motif IDs, then hits, then TFs that use
that motif. It emits one `EdgeRecord` per TF, motif and hit. The scanner is
therefore **already motif-level, not TF-level**: a motif shared by several TFs is
scanned once per promoter. What it did not do was reuse an identical matrix
stored under two motif IDs.

## Scanner mathematics

- `Pwm::log_odds` uses `log2(p / background)` and clamps zero probabilities to
  `1e-9`.
- `Scanner::new(pwm, background, scale=10.0)` builds:
  - a quantized per-position log-odds table;
  - an exact null distribution by dynamic programming over quantized scores;
  - an upper-tail table `P(null >= score)`.
- `Scanner::scan` skips any window containing a non-ACGT base.
- Both strands are scanned; reverse-strand hit coordinates are converted back to
  oriented-promoter coordinates.
- `Scanner::min_pvalue()` is the exact minimum attainable p-value for that PWM
  and background.

The p-value is exact for the **discretized** score and the supplied **zero-order
iid background**. It is not a p-value for a dinucleotide, GC-matched or
length-matched null.

## Problems found

### P1. Exact duplicate matrices were rescanned

Two motif IDs with bit-identical matrices each received a scanner and each was
scanned across every promoter. The duplication is not visible in the motif
library count but multiplies work and scan-cache size.

### P2. TF evidence rows duplicate a shared site

A single motif site used by several TFs produces one evidence row per TF. On
HF1 JASPAR 2024, 3,068,556 motif-level site hits became 3,090,137 evidence rows
(21,581 extra rows). The legacy report calls the former `prediction hits`, which
is easy to mistake for independent target evidence.

### P3. No site layer

`evidence.tsv` is a TF–motif–site table. There was no one-row-per-PWM-site table,
so the number of unique sites could not be recovered without post-hoc
deduplication, and exact duplicate PWMs could not be identified without rehashing
the library.

### P4. No persistent scan cache

A repeated run with the same promoters, motifs, background and threshold
recomputed every scan. There was no key that would make cache invalidation
auditable.

### P5. Background was numeric only

`--background` accepted four probabilities. There was no empirical promoter
background, no provenance of how it was estimated, and no cache key binding the
background values.

### P6. Testability was partially reported and hard-coded to uniform

`pipeline.rs` computed `Scanner::min_pvalue()` with `[0.25; 4]` even when the
scan used a different background. It reported `scanned_unreachable` and a
warning, but there was no per-motif testability table, and no status distinction
between "no significant hits" and "this motif could never reach the threshold".

### P7. No explicit mapping/scan/site/aggregation boundary

Mapping diagnostics existed (Phase 0–2.5), but the scan output did not carry a
stable PWM hash, a scan key, or a separate site table. Dashboard and export code
therefore had to infer these relationships.

## HF1 audit counts (before Phase 3)

| item | value |
|---|---:|
| genes / promoters | 37,415 / 37,415 |
| identified TFs | 1,907 |
| strict mapped TFs (JASPAR 2024) | 260 |
| mapped motif IDs | 257 |
| motif-level site hits | 3,068,556 |
| TF evidence rows | 3,090,137 |
| shared-motif extra rows | 21,581 |
| unmapped rows | 3 |
| motifs with p-value floor above 1e-4 | 28 |
| peak working set, legacy full scan | 3,308 MiB |
| legacy scan wall time | 69.3 s |

JASPAR 2026 strict had 271 TFs, 268 scanned motif IDs, 3,369,011 unique sites,
3,390,592 evidence rows and 29 motifs above the 1e-4 floor. Its original
`evidence.tsv` was 1,281,685,499 bytes; the Phase 3 optimized run reproduced the
same row set and byte size.

## What Phase 3 changed

- bit-identical normalized matrix hashing (`Pwm::digest`, SHA-256);
- opt-in `--reuse-exact-pwm` grouping;
- opt-in persistent `--scan-cache` keyed by algorithm, background, promoters,
  threshold, strand mode, min-score and PWM hash;
- `--sites-out` one-row-per-unique-PWM-site table with genomic coordinates when
  genome/GFF input is available;
- `--motif-testability-out` with attainable p-value, status and hit count;
- `--promoter-summary-out` with per-promoter site counts and best site;
- `--scan-manifest` with background, promoter hash, cache and reuse provenance;
- `--background-mode empirical` strand-symmetric zero-order promoter model;
- pipeline flags `--reuse-exact-pwm`, `--phase3-artifacts`,
  `--background-mode`, `--scan-cache`;
- fixed pipeline testability to use the effective background.

The default `scan` path remains the legacy in-memory sorted path. Phase 3 is
entered only when at least one Phase 3 option is supplied.

## Residual risks

- Site p-values are conditional on a zero-order background estimated globally
  (or supplied as uniform). They do not model dinucleotides, GC windows,
  chromatin accessibility or TF cooperation.
- Overlapping windows are dependent. No site-level or promoter-level FDR is
  claimed.
- `--promoters` FASTA input has no genomic coordinates; site rows then carry
  empty genomic coordinates. Use genome/GFF for coordinate-bearing output.
- The cache is content-keyed but has no eviction policy or cross-version
  migration. Changing the scanner algorithm requires a key-version bump.
