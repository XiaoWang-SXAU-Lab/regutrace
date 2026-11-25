# Phase 5 decisions

**Date:** 2026-11-25

## D1. Denominator is the PEP catalogue (1,907), not the CDS catalogue (1,908)

**Decision:** keep 1,907 for the current PEP-based HF1 background.

**Why:** the denominator must come from the TF identification catalogue, not from mapping rows. The PEP and CDS catalogues differ by exactly one gene, `Fta_6G0009940` (PEP ambiguous AP2;ERF; CDS identified AP2). No other identified gene differs. The current production pipeline used PEP.

**Consequence:** Phase 5 reports the reconciliation at gene-ID level; a CDS-based run may legitimately use 1,908. The denominator is a catalogue property, not a threshold choice.

## D2. Coordinate grouping is not the HF1 redundancy mechanism

**Decision:** do not reduce motifs by splitting query DBD coordinates on HF1.

**Why:** all coordinate-covered TFs have one query DBD region group at 0.5–0.9 overlap. The redundancy is donor/PWM alternatives for one region.

**Consequence:** best-supported and specificity-group policies are the useful levers; coordinate grouping is retained for future multi-domain genomes.

## D3. Best-supported means one motif per DBD group, not one per TF

**Decision:** select one motif per functional DBD architecture group; composite architectures remain one group.

**Why:** a multidomain TF may genuinely have more than one DNA-binding region. Collapsing to one motif per TF would discard a real architecture.

**Consequence:** HF1 best-supported has median 1 motif because most TFs have one group; future multidomain TFs can retain more.

## D4. Representative specificity groups are an intermediate policy

**Decision:** offer `reps@0.8` in addition to `best` and `all`.

**Why:** alternative specificity hypotheses may be biologically real; a single best motif can discard them. Specificity-group representatives keep alternatives without retaining every near-duplicate donor PWM.

**Consequence:** HF1 reps keeps 6 motifs/TF and 57.89% coverage, between all and best.

## D5. Preserve all candidate evidence

**Decision:** best/reps are query/export policies, never destructive filtering.

**Why:** downstream rsWGCNA methods may model motif dependence, and reproducibility requires the original mappings.

**Consequence:** the background DB stores `tf_motif` and `motif_target`; the best-policy `pairs` is an additional import contract, not a replacement for the full mapping.

## D6. No target-count or p-value ranking of motif support

**Decision:** representative ranking uses architecture, coverage, identity, bitscore, and deterministic tie-breaks only.

**Why:** ranking by the number of predicted targets or best promoter p-value would optimize the result being evaluated.

**Consequence:** ties are recorded; the top donor is not claimed to be biologically correct.

## D7. Multiple-testing correction is an optional interpretation layer

**Decision:** keep p=1e-4 as the scanner default; report site/promoter FDR separately; do not silently replace p-values.

**Why:** the complete 74M-window universe makes q≤0.05 extremely stringent; the user needs both the broad candidate background and a statistically interpretable layer.

**Consequence:** original p-values remain; FDR tables and scripts are additive. No scanner redesign.

## D8. Conservative BH, not an exact dependence model

**Decision:** use hit ranks plus total tested windows for site BH; Bonferroni windows plus BH across promoters for promoter-level.

**Why:** all windows with p≤cutoff are in the site table, so hit ranks are exact; non-hit terms can only lower q, making the hit-only result conservative. A full dependence model is not available.

**Consequence:** q-values are upper bounds, explicitly labelled. Do not present them as exact biological FDR.

## D9. No multiple-motif Bonferroni factor

**Decision:** do not multiply per-motif corrected p-values by the number of motifs assigned to a TF.

**Why:** motifs are correlated and a defensible dependence model is unresolved. Best/reps policies reduce the hypothesis count; the all layer remains for downstream modeling.

**Consequence:** per-motif q-values and TF-level policy comparisons are reported separately.

## D10. Standalone SQLite background with a policy-selected prior

**Decision:** build `tf_catalogue`, `tf_motif`, `motif_target`, `motif_registry`, empty canonical `hits`, and best-policy `pairs`.

**Why:** the existing rsWGCNA Data Import contract requires `hits`/`pairs`; materializing all 25M exploratory TF–target rows is unnecessary for import.

**Consequence:** the DB is 1.28 GiB and imports in ~61 s; `motif_target` preserves the 4.6M motif-level layer. Full exploratory pairs remain available through the all-policy analysis/export, not the default DB.

## D11. DAP-seq benchmark is promoter-overlap, not a negative set

**Decision:** use promoter overlap for positives and treat unpeaked genes as unknown.

**Why:** distal DAP peaks and unobserved genes cannot be assigned as negatives without additional data.

**Consequence:** PPV_assay/AUPRC are comparative diagnostics, not biological precision. The report labels this explicitly.

## D12. No threshold optimization on the benchmark

**Decision:** do not tune p-value, background, or motif-selection thresholds on the DAP results.

**Why:** that would overfit the benchmark and violate Phase 2.5/2.6 validation discipline.

**Consequence:** the benchmark reports the fixed Phase 3/4 configuration and its limitations.

## D13. Job registry remains in-memory; exported background is the handoff

**Decision:** do not rewrite ReguTrace job persistence in Phase 5.

**Why:** the acceptance requirement is a reusable background independent of the service/job ID. The exported SQLite plus rsWGCNA project state satisfies it; job-registry persistence is a separate engineering change.

**Consequence:** after a dashboard restart, the project must be reopened and old ReguTrace job IDs are not listed; the imported prior and upstream queries still work. This is documented, not hidden.

## D14. No Phase 6 work

**Decision:** stop after Phase 5.

**Why:** the task book ends at Phase 5. No WGCNA-specific regulatory-network algorithm was added to ReguTrace.

## D15. Freeze background schema v3 and exact policy associations

**Decision:** freeze the reusable background at schema **v3** / policy version **phase5-v1**. Store the exact all/best/reps TF–motif selections in `tf_motif_policy`, the coordinate rerun in `tf_motif_coordinate`, and policy parameters/provenance in `tf_motif_policy_config`; reconstruct target sets by joining to `motif_target` rather than materializing per-policy TF–target tables.

**Why:** the earlier v2 DB could not exactly reproduce representative-group target sets. Exact per-TF verification and the actual v3 import/reopen regression remove that ambiguity while keeping the artifact compact.

**Consequence:** downstream consumers must use the documented join. Any future change to selection rules requires a new policy version; any schema change requires a new schema version. Core prediction algorithms and thresholds remain frozen.
