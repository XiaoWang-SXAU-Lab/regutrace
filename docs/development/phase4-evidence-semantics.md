# Phase 4 evidence semantics

**Date:** 2026-11-25

## 1. Four independent coverage layers

The UI and exports must not collapse these into one “coverage” number:

| Layer | Question | Failure state |
|---|---|---|
| TF identification | Was a TF gene identified in this genome? | not in catalogue |
| Motif mapping | Is there a supported TF–motif assignment? | `no_supported_motif` |
| Motif testability | Can the PWM attain the requested p-value under the background? | `threshold_unattainable` / `not_scanned` |
| Target prediction | Does a testable motif produce at least one promoter hit? | zero targets at the requested cutoff |

A TF with a motif but no target is not the same as a TF with no motif. A TF with a motif that cannot attain the cutoff is not the same as a TF whose scan found nothing. The catalogue exposes a diagnostic reason for each case.

## 2. Mapping statuses

| Status | Meaning | Default strict | Default balanced | Default exploratory |
|---|---|:---:|:---:|:---:|
| `strict_rbh` | Full-length protein reciprocal best hit; unchanged Phase 1 semantics | yes | yes | yes |
| `dbd_calibrated` | Phase 2.5 calibrated DBD transfer | no | yes | yes |
| `dbd_ambiguous` | DBD evidence does not resolve a specificity group | no | no | yes (visible; excluded from coverage/prior by default) |
| `dbd_exploratory` | Sequence-supported DBD candidate; uncalibrated | no | no | yes |
| `no_supported_motif` | Identified TF with no supported assignment | visible | visible | visible |
| `not_computed` | No mapping table was supplied | visible | visible | visible |

`strict_calibrated` is the prior-export policy name for strict + calibrated; it is not a fourth mapping mode.

## 3. Testability statuses

| Status | Meaning |
|---|---|
| `significant_hits` | At least one site passes the requested cutoff |
| `no_significant_hits` | Scan completed; no site passes |
| `threshold_unattainable` | The minimum attainable p-value for this PWM/background/length is above the requested cutoff |
| `not_scanned` | Motif was not part of the scan |
| `not_computed` | No testability table |

A motif that is `threshold_unattainable` must not be reported as “no target”; it was never testable at that cutoff.

## 4. Evidence levels

- **Level D** — computational PWM promoter hit. All Phase 3/4 target associations are level D.
- The UI never upgrades a PWM hit to experimental validation.
- Mapping status describes how the TF–motif assignment was obtained; it does not validate the target.
- Shared effective PWMs mean multiple TFs can explain the same site; the UI shows a `shared PWM × N` warning and does not treat the TFs as independent evidence.

## 5. Coverage definitions

- **Identified TFs:** unique TF gene IDs in the identified-TF catalogue after identity mapping.
- **TFs with candidate motifs:** identified TFs with at least one mapping row in the current mode.
- **TFs with testable motifs:** mapped TFs with at least one motif whose status is `significant_hits`.
- **TFs with completed scanning:** mapped/testable TFs whose scan completed.
- **TFs with ≥1 candidate target:** TFs with at least one promoter hit at the requested cutoff.
- **Actual target coverage:** TFs with ≥1 candidate target divided by identified TFs.
- **Motif–target associations:** distinct motif–target pairs with at least one site.
- **TF–target candidate associations:** distinct TF–target pairs after expanding motifs to mapped TFs.
- **Median targets per TF:** median over TFs with ≥1 target.

Exploratory coverage counts ambiguous TFs in “candidate motifs” but excludes them from “testable”, “≥1 target”, and prior export by default. The HF1 table reports 11 ambiguous TFs separately.

## 6. Mode policy

### Strict

- Existing full-length protein RBH only.
- No DBD inference.
- Default mode; preserves legacy results.

### Balanced

- Strict plus Phase 2.5 calibrated DBD transfer.
- Calibration thresholds are not reused for non-reciprocal inference.
- Conservative by design; HF1 adds only 9 TFs over strict.

### Exploratory

- Sequence-supported DBD candidates, non-reciprocal donors, alternative specificity groups, uncalibrated candidates.
- Never labeled calibrated.
- Candidate targets are viewable, but the UI labels them exploratory.
- Not automatically promoted to balanced or strict.

## 7. Background and p-values

- Default background is uniform A/C/G/T = 0.25.
- Empirical background is a zero-order model estimated from the run's promoters; it does not model dinucleotides.
- The requested p-value is not relaxed automatically.
- If a motif cannot attain the cutoff, the testability table and UI say so.
- The same p-value has different attainable floors for different motif lengths and base compositions.

## 8. Enrichment semantics

The enrichment panel is exploratory:

- one-sided Fisher exact test on binary promoter hits;
- BH across the requested motifs only;
- background can be user-supplied or all scanned promoters;
- no GC/length matching, gene-set overlap correction, or motif-similarity correction;
- shared PWMs create dependence across TFs;
- a small q-value is not evidence of experimental regulatory enrichment.

## 9. Export semantics

- `strict`: strict RBH only; all motifs must be testable.
- `strict_calibrated`: strict + calibrated.
- `extended`: strict + calibrated + exploratory; ambiguous excluded.
- `selected`: explicit TF IDs; strict/calibrated included; exploratory only with `include_exploratory=true`.
- All exported prior rows remain level D.
- Site-only exploratory jobs have no materialized TF–target prior.
