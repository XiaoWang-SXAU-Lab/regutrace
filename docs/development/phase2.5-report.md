# Phase 2.5 report — DBD motif-transfer calibration

Date: 2026-11-24  
Repository: `regutrace`  
Base: Phase 0–2 working tree (base HEAD `9ac4837`; not committed)

## Outcome

Phase 2.5 is implemented with explicit abstention rather than a universal
threshold.

- Calibration dataset: **596 donors / 619 experimental motifs / 32 families**.
- **11 families** met the sample-size floor and were tested.
- **2 families** met the nested validation criteria:
  - `HD-ZIP`: identity ≥75%, coverage ≥90%, nested OOF precision 1.00
    (95% CI lower 0.70), calibrated coverage 26.5% (9/34).
  - `TCP`: identity ≥70%, coverage ≥50%, nested OOF precision 1.00
    (95% CI lower 0.87), calibrated coverage 71.4% (25/35).
- **30 families abstain**, either before threshold search (donor/species floor)
  or because nested validation did not support a calibrated call.
- HF1 strict + calibrated transfer maps **280 TFs** versus 271 for strict 2026:
  **+9 HD-ZIP TFs**; calibrated TCP adds no TF beyond strict 2026.
- No exploratory mapping is promoted by coverage alone.

## Calibration dataset

Outputs:

```text
calibration-donors.tsv          596 donors
calibration-motifs.tsv          619 donor-associated experimental motifs
calibration-dbd.fasta           596 DBD sequences
calibration-family-summary.tsv  32 families
calibration-exclusions.tsv      351 excluded donor/motif records
exact-pwm-groups.tsv            947 exact PWM groups over the 2024/2026 union
```

Exclusion counts:

| reason | records |
|---|---:|
| JASPAR 2026 motif not experimental (empty/Transfac) | 97 |
| JASPAR 2024 motif not experimental (empty/Transfac) | 96 |
| JASPAR 2024 family unresolved | 34 |
| JASPAR 2026 family unresolved | 29 |
| no DBD sequence (abstained family) | 33 |
| JASPAR 2026 donor not classified | 26 |
| JASPAR 2024 accession not retrieved | 13 |
| JASPAR 2024 taxonomy unresolved | 12 |
| JASPAR 2026 taxonomy unresolved | 11 |

Each retained motif has one donor; exact motif ID recovery is zero by
construction. The calibration target is therefore PWM similarity, not exact ID
recovery.

Family-level sample sizes are shown in
[family-sample-sizes.svg](figures/phase2_5/family-sample-sizes.svg). The 11
supported families are Dof, ERF, HD-ZIP, MYB, MYB-related, NAC, SBP, TCP, WRKY,
bHLH and bZIP. The rest abstain because they have fewer than 15 donors, fewer
than 3 species, or fewer than 2 exact PWM groups.

## Pairwise DBD and motif similarity

Within-family donor pairs: **13,848**; **11,784** have a motif similarity with
≥50% overlap and both identity scores.

Predictor comparison for the label `motif similarity ≥ 0.80`:

| predictor | rank AUC | bootstrap 95% CI |
|---|---:|---:|
| DBD identity | **0.696** | 0.687–0.705 (pair-level; see Phase 2.6) |
| full-length protein identity | 0.608 | 0.598–0.618 (pair-level; see Phase 2.6) |

Phase 2.6 showed that pair-level bootstrap is pseudoreplicated. The
homology-cluster bootstrap intervals are 0.677–0.776 for DBD identity and
0.572–0.704 for full-length protein identity; use those for uncertainty
reporting.

DBD identity is better overall, but not universally. Family-specific AUCs range
from 0.35 (MIKC) to 0.99 (AP2); ARF, Dof, SBP and MIKC have AUC below 0.5,
meaning higher DBD identity does not monotonically imply higher PWM similarity
within those families. This is why no universal identity threshold was used.

Figures:

- [DBD identity vs PWM similarity](figures/phase2_5/dbd-identity-vs-motif-similarity.svg)
- [Family-specific panels](figures/phase2_5/family-dbd-motif-similarity.svg)
- [DBD vs protein predictor AUC](figures/phase2_5/predictor-auc.svg)

## Motif specificity groups

Using TOMTOM-style Pearson correlation with ≥50% overlap:

| threshold | groups | motifs | singleton groups |
|---|---:|---:|---:|
| 0.80 | 127 | 619 | 29 |
| 0.90 | 174 | 619 | 61 |

Groups are operational, not biological equivalence classes. Exact PWM hashes
are always merged. The largest 0.90 groups are WRKY (29 motifs), ERF (17), TCP
(16), bHLH (12) and bZIP (12). Full tables:

```text
motif-specificity-groups.tsv
specificity-group-summary.tsv
motif-pairwise-similarity.tsv
```

## Threshold calibration

Primary evaluation is nested species-grouped cross-validation. The final
deployment threshold is selected on all donors with the same leave-one-species-out
rule. Promotion requires nested OOF precision ≥0.90, Wilson lower bound ≥0.70,
and at least `max(5, 10% of donors)` calibrated calls.

| family | identity | coverage | donors | calibrated OOF | OOF precision | OOF coverage | OOF precision CI |
|---|---:|---:|---:|---:|---:|---:|---:|
| HD-ZIP | 75% | 0.90 | 34 | 9 | 1.00 | 0.265 | 0.70–1.00 |
| TCP | 70% | 0.50 | 35 | 25 | 1.00 | 0.714 | 0.87–1.00 |

Precision–coverage tradeoff:
[precision-coverage.svg](figures/phase2_5/precision-coverage.svg).

Abstentions among the 11 sample-size-supported families:

| family | reason |
|---|---|
| Dof | nested validation precision below target |
| ERF | no nested species calls met the criterion |
| MYB | no nested species calls met the criterion |
| MYB-related | insufficient nested calls |
| NAC | no nested species calls met the criterion |
| SBP | no nested species calls met the criterion |
| WRKY | nested validation precision below target |
| bHLH | insufficient nested calls |
| bZIP | insufficient nested calls |

The remaining 21 families abstain before threshold search due to the
donor/species/exact-group floor.

## HF1 reassessment

Methods:

1. JASPAR 2024 strict RBH;
2. JASPAR 2026 strict RBH;
3. uncalibrated DBD transfer (2026 donors);
4. calibrated DBD transfer;
5. strict 2026 plus calibrated DBD transfer.

| method | mapping rows | mapped TFs | calibrated TFs | ambiguous TFs | unique motifs | strict-2026 overlap | exact donor | exact motif | motif similar ≥0.80 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| strict 2024 | 260 | 260 | 0 | 0 | 257 | — | — | — | — |
| strict 2026 | 271 | 271 | 0 | 0 | 268 | — | — | — | — |
| DBD uncalibrated | 32,765 | 1,125 | 0 | 0 | 657 | 253 | 251 | 251 | 253 |
| DBD calibrated | 57 | 19 | 19 | 11 | 57 | 10 | 9 | 9 | 10 |
| strict 2026 + calibrated | 280 | 280 | 19 | 11 | 275 | 271 | 271 | 271 | 271 |

The calibrated DBD method maps only HD-ZIP and TCP. It recovers 10 strict-2026
TFs; 9 have the exact strict donor and motif and all 10 are within the ≥0.80
similarity label. Its unique contribution is 9 HD-ZIP TFs not present in strict
2026. The uncalibrated method adds 872 TFs relative to strict 2026, but those are
diagnostic candidates and are not promoted.

Calibrated diagnostics over 1,981 catalogue TFs:

| reason | count |
|---|---:|
| NO_CALIBRATED_THRESHOLD | 1,484 |
| NO_DBD_SEQUENCE | 341 |
| NOT_IDENTIFIED | 74 |
| NOT_RECIPROCAL_DBD_BEST | 29 |
| BELOW_DBD_THRESHOLD | 23 |
| DBD_CALIBRATED | 19 |
| DBD_AMBIGUOUS | 11 |

Runtime (mapping only, wall time, same machine):

| method | wall time | peak sampled working set |
|---|---:|---:|
| strict 2024 / 2026 | ~44 s each | not instrumented |
| uncalibrated DBD | 34.4 s | not instrumented |
| calibrated DBD | ~19 s | 78.8 MiB (parent + blastp sampling) |

The calibrated map is not a genome-wide target release. It is a mapping
reassessment; Phase 3 will not scan all 32,765 exploratory rows.

### JASPAR 2026 strict changes

The 2024 → 2026 strict comparison has 54 changed TFs:

- 21 new mappings;
- 10 lost mappings;
- 23 changed donor/motif pairs (20 donor changes, 2 version updates, 1
  representative change).

Lost mappings split into:

- 4 donors absent from the 2026 pack;
- 1 motif absent from the 2026 library;
- 5 competition/threshold changes where both donor and motif are present.

Full table: `hf1-strict-2024-vs-2026-changes.tsv`. The combined
strict+calibrated row set is `hf1-strict-plus-calibrated.tsv`. Agreement with strict RBH is
consistency, not ground truth; neither changes nor losses are automatically
errors.

## Short-motif p-value floor

The 2026 HF1 scan scanned 268 motifs. **29** have a minimum attainable p-value
above the default `1e-4` under the uniform background:

| item | value |
|---|---|
| motifs | 29 |
| length 4 | 5 |
| length 5 | 1 |
| length 6 | 22 |
| length 7 | 1 |
| information content | 3.996–11.49 bits |
| minimum p-value | 1.22e-4 – 3.91e-3 |
| floor / 1e-4 | 2.44–39.06 |

Full table: `phase3-motif-threshold-audit.tsv`. These motifs are short, so the
discretized null cannot produce a p-value below the cutoff. The scan warns and
continues; the cutoff is not silently relaxed. Phase 3 must expose the attainable
floor per motif and choose an explicit policy.

## Code and tests

New/updated code:

- `dbd-extract` and `map-dbd-motifs` extensions;
- calibrated threshold config support;
- specificity groups and reciprocal DBD-best filter;
- mapping statuses `dbd_calibrated`, `dbd_ambiguous`, `dbd_exploratory`,
  `no_supported_motif`;
- calibration scripts in `scripts/phase2_5/`.

Rust tests:

```text
40 CLI unit tests passed
7 CLI integration tests passed
69 core unit tests passed
0 failed
```

Python synthetic tests:

```text
7 tests passed (scripts/phase2_5/test_phase2_5.py)
```

No strict RBH code path or output semantics changed.

## Calibration configs

```text
config/dbd-models/plant-calibrated-v1.json
config/dbd-thresholds/plant-calibrated-v1.json
```

Only HD-ZIP and TCP appear. All other families are absent by design, so a
calibrated run returns `NO_CALIBRATED_THRESHOLD` rather than transferring an
unvalidated motif.

## Reproducibility and hashes

Calibration inputs and commands are documented in
[dbd-calibration-methodology.md](dbd-calibration-methodology.md). Key output
hashes:

```text
calibration-donors.tsv          2c3f233f03992aa0d35a92bb1e54abffabcf4bf4c8312056749730b1d119194b
calibration-motifs.tsv          3b0a8503c163f01f11d5743e8f8d84d06077a37a0a219c886f0f2f8716e06b10
calibration-family-summary.tsv  b569a2dc0f89ddabe35048261fca66665c3b2fc6871304a4f840ad2b3b031381
calibration-thresholds.tsv      c083a3979ddf7dcf4f00f901d753a22ea89672f15b0feccb1a700fa5de93790d
hf1-method-summary.tsv          b9843dc1426c7b8b1cce2aa956197c3f81db20d63bb1135090addc48c7ee25a8
```

The machine-readable calibration bundle is copied to the sync workspace at:

```text
analysis/regutrace-yeast-validation/phase2.5-calibration-2026-11-24/
```

## Limitations and blockers

1. **Small independent test sets.** JASPAR plant CORE is Arabidopsis-dominated;
   species holdouts for several families contain only 1–5 donors. This is why
   most families abstain.
2. **Motif similarity is operational.** 0.80 is not a biological equivalence
   threshold; 0.90 sensitivity is reported. TOMTOM E-values were not computed.
3. **Method metadata is incomplete.** 95 representative 2026 records have no
   method; they are excluded from ground truth. This reduces calibration data.
4. **DBD boundaries are curated.** C2H2/C3H array extraction remains abstained.
5. **No DNA-contact-residue model.** The task suggested considering recognition
   residues where supported; no structure/contact annotation was available, so
   this is not implemented.
6. **No biological benchmark.** HF1 overlap is consistency with strict RBH, not
   validation against DAP-seq/ChIP-seq targets.
7. **No test-set leakage.** Donors and species are held out; exact releases are
   merged; final thresholds use all donors but the reported estimate is nested
   OOF.
8. **Reciprocity panels differ.** Calibration checks a donor's reverse best
   within the donor panel plus the query; HF1 checks against all local query
   DBDs. The HF1 criterion is stricter and may reduce calibrated coverage.

## Next steps

1. Phase 3: separate mapping, motif–site hits, motif–target aggregation, and
   TF–target associations; exact PWM scan reuse by hash/config.
2. Phase 3: empirical promoter background and per-motif attainable p-value
   policy; do not silently relax `1e-4`.
3. Phase 4: dashboard shows `dbd_calibrated`, `dbd_ambiguous`,
   `dbd_exploratory` separately; no ambiguous edge enters prior export.
4. Phase 5: benchmark calibrated HD-ZIP/TCP against DAP-seq/ChIP-seq and test
   whether adding full-length protein or DNA-contact features improves precision.
