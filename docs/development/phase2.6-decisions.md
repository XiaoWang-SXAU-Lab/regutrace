# Phase 2.6 decisions

Date: 2026-11-24

## D1. Phase 2.5 calibrated models are preserved

`config/dbd-models/plant-calibrated-v1.json` and
`config/dbd-thresholds/plant-calibrated-v1.json` are unchanged. Phase 2.6 is an
ablation and audit; it does not overwrite the HD-ZIP/TCP calibrated configs.

## D2. Reciprocity is a precision filter, not a CIS-BP requirement

CIS-BP 2014 describes nearest-neighbor transfer and does not require a reciprocal
best hit. ReguTrace's reciprocal DBD method (B) is therefore an additional
confidence filter, not a reproduction of CIS-BP.

The ablation shows the tradeoff under species-held-out evaluation:

| method | precision | coverage |
|---|---:|---:|
| B DBD reciprocal | 0.933 | 0.151 |
| D DBD non-reciprocal, family thresholds | 0.846 | 0.272 |

Reciprocity raises precision but roughly halves coverage. This is a real
precision–coverage tradeoff, not an implementation defect.

## D3. No non-reciprocal method is promoted

Methods C, D and E are not promoted to a production/default config because:

- C's single global threshold met the precision target in only one outer fold
  (species-held-out coverage 0.072, precision 0.767);
- D has higher coverage but species-held-out OOF precision 0.846 and
  homology-cluster-held-out precision 0.727;
- E's unique-group rule has species-held-out precision 0.864, coverage 0.074,
  and 7.0% ambiguity.

The default remains strict RBH; HD-ZIP/TCP calibrated DBD remains the only
calibrated mode.

## D4. Validation settings are not interchangeable

Donor-held-out, species-held-out and homology-cluster-held-out give materially
different estimates:

| method | donor precision / coverage | species | homology cluster |
|---|---:|---:|---:|
| A | 0.937 / 0.455 | 0.878 / 0.233 | 0.870 / 0.077 |
| B | 0.945 / 0.369 | 0.933 / 0.151 | 0.825 / 0.096 |
| C | 0.911 / 0.681 | 0.767 / 0.072 | no calls |
| D | 0.922 / 0.622 | 0.846 / 0.272 | 0.727 / 0.258 |
| E | 0.978 / 0.225 | 0.864 / 0.074 | 0.750 / 0.054 |

Donor-held-out folds retain same-species and homologous donors and are
optimistic. Homology-cluster-held-out is the strictest setting. Species-held-out
is the primary setting for the Phase 2.6 comparison, but it is not equivalent to
an external dataset.

## D5. Pairwise confidence intervals must be cluster-aware

Phase 2.5 bootstrapped donor pairs as independent observations. Phase 2.6 shows
that this understates uncertainty:

| predictor | pair bootstrap CI width | donor bootstrap | homology-cluster bootstrap |
|---|---:|---:|---:|
| DBD identity | 0.018 | 0.062 | 0.099 |
| full-length protein identity | 0.021 | 0.064 | 0.132 |

53% of motif-similarity pairs are within a 70% DBD homology cluster and 56% are
same-species pairs. Pair-level CIs are therefore pseudoreplicated. Phase 2.6
reports cluster bootstrap intervals for AUC and method precision/coverage.

## D6. 0.80 is the primary motif-similarity label; 0.90 is sensitivity

At motif similarity ≥0.90, precision drops for every method (for example B:
0.933 → 0.856; A: 0.881 → 0.713). The 0.80 label is retained as the primary
operational label, with 0.90 reported separately. Neither is biological
equivalence.

## D7. Exact motif recovery is not a useful calibration target

Each calibration motif has one donor. Held-out donors cannot recover an exact
motif ID from another donor. Exact recovery is zero by construction; PWM
similarity is the primary label.

## D8. Family-support floor is donor-count driven

Across the tested grid, the minimum exact-group count (1/2/3) never changed the
supported-family set, and the minimum species count changed it only at the
10-donor floor. The donor count is the binding constraint. The Phase 2.5 floor
(15 donors / 3 species / 2 exact groups) is retained for comparability; Phase
2.6 does not lower it to increase coverage.

## D9. Representative motif selection is not a major factor

Using all donor motifs versus the maximum-information-content representative
changed the correctness label for 0/143 (A), 0/90 (B), 0/162 (D) and 1/44 (E)
calls. This does not justify collapsing motifs; it only shows that the
representative choice is not driving the ablation.

## D10. No UI or genome-wide scan in Phase 2.6

The ablation is donor-level calibration on experimental reference data. It does
not scan the 32,765 exploratory HF1 rows and does not change dashboard behavior.
Phase 3 will define the scan/aggregation interface and Phase 5 will benchmark
target-level performance.
