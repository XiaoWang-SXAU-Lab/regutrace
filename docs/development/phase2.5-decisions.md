# Phase 2.5 decisions

Date: 2026-11-24

This file records the decisions that materially affect calibration results. It
is intended to prevent later sessions from silently changing the criteria or
reading exploratory results as validated results.

## D1. Calibration unit is donor-held-out motif transfer

A query donor is removed from the reference panel; its experimentally supported
motif set is the label. The prediction comes from same-family reference donors
only. This avoids testing a donor against itself and avoids treating overlapping
JASPAR releases as independent datasets.

## D2. Exact duplicates are not independent

JASPAR 2024 and JASPAR 2026 share 781 matrix IDs with identical PFMs. Donor–
motif records are merged by `(donor_id, motif_id)`. Exact duplicate PWMs are
grouped by canonical normalized matrix hash. They remain visible as provenance
but do not become independent positive evidence.

## D3. Only experimental methods enter calibration

DAP-seq, ampDAP-seq, PBM, ChIP-seq, ChIP-chip, ChIP-exo, SELEX, SELEX-seq,
HT-SELEX, CAP-SELEX, bacterial 1-hybrid and the PBM/CSA/DIP-chip composite are
accepted. `Transfac`, `COMPILED` and empty method fields are excluded from
ground truth and listed in `calibration-exclusions.tsv`. This is why the
calibration set is smaller than the full JASPAR CORE pack.

## D4. DBD families with ambiguous architecture abstain

The following families do not receive a DBD sequence or a threshold:

```text
C2H2, C3H, CO-like, CPP, GRAS, LSD, NF-YA, NF-YB, NF-YC, ZF-HD
```

Reasons: variable zinc-finger arrays, composite DNA binding, or a DBD whose
boundary is not encoded by the pack rule. No fallback boundary was invented.

## D5. TOMTOM correlation is used, not TOMTOM significance

The offline build has no TOMTOM binary. We use the Pearson-correlation component
of the TOMTOM approach across offsets and strands. It is an operational
similarity measure, not an E-value and not proof of biological equivalence.
Overlap must be at least 50% of the shorter motif.

## D6. Primary motif-similarity label is 0.80; 0.90 is sensitivity

The 0.80 threshold is a pre-registered operational label chosen from the
within-family distribution. It is not a universal biological constant. Results
are also reported at 0.90; no result was selected by maximizing coverage at an
arbitrary similarity cutoff.

## D7. Calibrated transfer requires a reciprocal DBD top hit

For a calibrated call, the unique top donor must also have the query as its best
DBD hit within the reference panel plus the query. This is a confidence filter,
not a claim that DBD RBH is biological truth. It reduced coverage but removed
many non-reciprocal top-donor errors in the calibration folds.

## D8. Ambiguity is not promoted

If candidate motifs fall into more than one specificity group, or a candidate
motif lacks a group annotation, the TF is `dbd_ambiguous`. Ambiguous rows remain
visible with all donors/motifs but are excluded from calibrated coverage and
from the combined strict+calibrated set.

## D9. Family support minimum is 15 donors / 3 species / 2 exact groups

This is a sample-size and diversity floor, not a biological threshold. Families
below it abstain before threshold search. The floor prevents a family with one
species or one motif from producing a seemingly precise threshold.

## D10. Nested species-grouped CV is the primary evaluation

Threshold selection is nested inside outer species folds. Promotion requires
nested out-of-fold precision ≥ 0.90, Wilson lower bound ≥ 0.70, and at least
`max(5, 10% of donors)` calibrated calls. This is deliberately conservative
because JASPAR plant CORE is Arabidopsis-dominated and species holdouts are
small.

## D11. Only two families met the promotion criteria

`HD-ZIP` and `TCP` are calibrated. The other supported families abstain because
their nested validation precision or number of calibrated calls was
insufficient. This is an intentional outcome; there is no minimum coverage
requirement.

## D12. The strict map is not ground truth

Agreement with strict RBH is reported as consistency, not accuracy. Disagreements
are not automatically errors or improvements. The 10 lost and 23 changed
JASPAR 2026 strict mappings are listed for review.

## D13. The 29 short-motif p-value floors are not relaxed

For 29 scanned motifs, the minimum attainable p-value under the configured
uniform background is above `1e-4`. The scan warns and continues; it does not
silently lower the cutoff. Phase 3 must address this with a per-motif attainable
threshold report and an explicit statistical decision.

## D14. No genome-wide scan of exploratory mappings

The 32,765 uncalibrated HF1 candidate rows are a diagnostic set, not a release
target set. Phase 3 will separate mapping, motif–site hits, motif–target
aggregation and TF–target associations, and will not scan all exploratory rows.

## D15. Calibration reciprocity is donor-panel based; HF1 is all-local based

In calibration, the reciprocal DBD check compares a donor against the other
calibration donors plus the query. In HF1 mapping, `donor-dbd` is aligned against
all local query DBDs, so a non-donor local TF can take the donor's reverse top
hit. The HF1 criterion is therefore stricter than the calibration criterion and
may reduce calibrated coverage. This is recorded as a limitation; it was not
changed after seeing the HF1 result.

## D16. Pair-level bootstrap intervals are superseded

Phase 2.6 showed that 53% of motif-similarity donor pairs are within a 70% DBD
homology cluster and 56% are same-species pairs. Treating pairs as independent
understates the DBD AUC confidence interval by roughly 5–6×. Use the Phase 2.6
homology-cluster bootstrap intervals (DBD AUC 0.677–0.776; full-length protein
AUC 0.572–0.704) for uncertainty reporting.

