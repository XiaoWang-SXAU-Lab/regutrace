# Phase 2.5 HF1 comparison

Date: 2026-11-24  
Full methodology and calibration results:
[phase2.5-report.md](phase2.5-report.md)

## Mapping methods

| method | mapping rows | mapped TFs | calibrated TFs | ambiguous TFs | unique motifs | overlap with strict 2026 | exact donor | exact motif | motif similarity ≥0.80 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| JASPAR 2024 strict RBH | 260 | 260 | 0 | 0 | 257 | — | — | — | — |
| JASPAR 2026 strict RBH | 271 | 271 | 0 | 0 | 268 | — | — | — | — |
| uncalibrated DBD | 32,765 | 1,125 | 0 | 0 | 657 | 253 | 251 | 251 | 253 |
| calibrated DBD | 57 | 19 | 19 | 11 | 57 | 10 | 9 | 9 | 10 |
| strict 2026 + calibrated DBD | 280 | 280 | 19 | 11 | 275 | 271 | 271 | 271 | 271 |

The calibrated method maps only HD-ZIP and TCP. It adds 9 HD-ZIP TFs beyond
strict 2026; calibrated TCP adds no TF beyond strict 2026. The uncalibrated
method adds 872 TFs relative to strict 2026, but those are diagnostic candidates
and are not promoted.

## Family coverage for calibrated families

| family | catalogue TFs | strict 2024 | strict 2026 | uncalibrated DBD | calibrated DBD | ambiguous DBD | combined |
|---|---:|---:|---:|---:|---:|---:|---:|
| HD-ZIP | 52 | 11 | 13 | 49 | 16 | 5 | 22 |
| TCP | 30 | 7 | 8 | 19 | 3 | 6 | 8 |

No other family has a calibrated threshold, so calibrated mapping returns
`NO_CALIBRATED_THRESHOLD` for them.

## Calibrated diagnostics

Over 1,981 HF1 catalogue TFs:

| reason | count |
|---|---:|
| NO_CALIBRATED_THRESHOLD | 1,484 |
| NO_DBD_SEQUENCE | 341 |
| NOT_IDENTIFIED | 74 |
| NOT_RECIPROCAL_DBD_BEST | 29 |
| BELOW_DBD_THRESHOLD | 23 |
| DBD_CALIBRATED | 19 |
| DBD_AMBIGUOUS | 11 |

## JASPAR 2024 → 2026 strict changes

54 TFs differ:

- 21 newly mapped;
- 10 lost;
- 23 changed donor or motif.

Lost mappings:

| reason | count |
|---|---:|
| donor absent from 2026 pack | 4 |
| motif absent from 2026 library | 1 |
| competition/threshold change with donor and motif present | 5 |

Changed mappings: 20 donor changes, 2 motif version updates, 1 representative
motif change. Neither loss nor change is automatically an error or an
improvement; these are candidates for Phase 5 target-level benchmarking.

## Runtime

Mapping only, same machine:

| method | wall time |
|---|---:|
| strict 2024 | ~44 s |
| strict 2026 | ~44 s |
| uncalibrated DBD | 34.4 s |
| calibrated DBD | ~19 s |

Peak sampled working set for the calibrated DBD run was 78.8 MiB
(parent plus `blastp` sampling). Peak memory for the other mapping runs was not
instrumented.

## Interpretation

- The uncalibrated DBD method has high coverage but is not a release method.
- The calibrated method has low coverage but its family-specific validation
  passed the pre-registered criteria.
- Agreement with strict RBH is consistency, not accuracy.
- The 11 ambiguous calibrated-family TFs remain visible but are excluded from
  calibrated and combined counts.
- No exploratory row set was genome-wide scanned in Phase 2.5.
