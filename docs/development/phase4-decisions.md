# Phase 4 decisions

**Date:** 2026-11-25

## D1. Keep the existing dashboard and extend it

**Decision:** add Phase 4 views to the existing ReguTrace toolkit page and backend rather than build a separate dashboard.

**Why:** the page already had upload/run/progress/query/promoter/diagnostics and a disk-backed evidence index. A rewrite would have risked regressions and split the user workflow.

**Consequence:** Phase 4 additions share the existing API/client types and visual language.

## D2. Site-only table instead of TF × target × site materialization

**Decision:** exploratory scans write one site row per unique PWM group/site. TF identity is resolved at query time.

**Why:** materializing all TF–target–site rows for 1,143 TFs would create hundreds of millions of rows and a multi-GB-to-tens-of-GB table, with no scientific gain. The user explicitly forbade this default.

**Consequence:** site-only jobs cannot export a TF–target prior directly; the UI explains this and offers the site table plus TF–motif mapping.

## D3. Strict/balanced/exploratory share one scan

**Decision:** run one exploratory site scan and intersect with mode-specific TF–motif maps for strict/balanced/exploratory coverage.

**Why:** scanning three times would produce identical sites for shared motifs and waste hours. The site table is mode-independent; the mapping is not.

**Consequence:** the coverage table reports the same 8,191,211 sites for all three modes, with different TF/motif intersections. This is stated explicitly.

## D4. Denominator is 1,907, not 1,908

**Decision:** use the actual identified-TF count (1,907) in all coverage percentages.

**Why:** the pipeline's identity summary maps 32,792 TF–motif rows to 1,907 genes. Using 1,908 from the task book would make the table internally inconsistent and slightly understate coverage.

**Consequence:** the report notes the 1-gene discrepancy; it is <0.01 percentage points and does not change any conclusion.

## D5. Exploratory coverage is not accuracy

**Decision:** report 58.94% as *actual target coverage*, never as precision or validation.

**Why:** exploratory mappings are sequence-supported DBD candidates, not calibrated assignments. Phase 2.5's ~90% precision applies to high-confidence calibrated models, not to this layer.

**Consequence:** all exploratory targets remain level D; the UI labels the mode and never auto-promotes.

## D6. No threshold re-optimization

**Decision:** keep the requested p-value and report unattainable motifs instead of lowering thresholds per family.

**Why:** re-tuning on the calibration dataset would invalidate the Phase 2.5/2.6 validation and inflate coverage.

**Consequence:** 65 motifs are `threshold_unattainable` at p=1e-4 on HF1; they are visible, not hidden.

## D7. Bounded TF target cache

**Decision:** cache up to four aggregated TF target results per job, keyed by TF/status/max-p/mode.

**Why:** without caching, every page of a 36k-target TF took ~6.3 s. Caching makes pagination ~1–2 ms with bounded memory.

**Consequence:** first query still costs seconds; cache is per job and lost on restart.

## D8. Separate background files, not one archive

**Decision:** export manifest, catalogue, TF–motif, site table, and prior separately.

**Why:** the site table is ~2 GB; bundling it would duplicate data and slow every download. The manifest points to the pieces and carries checksums.

**Consequence:** users must keep the set together; the manifest makes the relationship explicit.

## D9. Prior policies

**Decision:** implement `strict`, `strict_calibrated`, `extended`, and `selected`, with ambiguous excluded and exploratory opt-in only for `selected`.

**Why:** rsWGCNA users need a conservative default and an explicit exploratory option. Silent inclusion of exploratory rows would violate the evidence boundary.

**Consequence:** `extended` is available for non-site-only exploratory jobs; site-only jobs have no prior.

## D10. Shared-PWM warning rather than deduplication

**Decision:** show how many TFs share an effective PWM; do not merge them or drop duplicates.

**Why:** merging would hide ambiguity and change the biological interpretation. The user must see that multiple TFs can explain one site.

**Consequence:** target/regulator tables can show several TFs with the same p-value; a warning explains why.

## D11. Keep JASPAR 2024 as default

**Decision:** do not switch the default to JASPAR 2026 in Phase 4.

**Why:** Phase 3 validated the 2024 baseline; switching would confound Phase 4 coverage and require a new audit. 2026 remains a candidate for Phase 5.

## D12. Keep the legacy free-query path

**Decision:** do not remove `/query`; add dedicated `/targets`, `/regulators`, and motif endpoints.

**Why:** existing users and links depend on it; the dedicated endpoints provide the scalable path for new views.

**Consequence:** the legacy path remains slower for large TFs and is documented as a bottleneck.

## D13. In-memory job registry

**Decision:** accept the existing manager design in Phase 4; do not persist/reload job metadata.

**Why:** changing job persistence is outside Phase 4 scope and risks the run pipeline.

**Consequence:** restarting the dashboard loses job IDs even though files remain. This is documented as an operational limitation; a future phase should add job metadata persistence.

## D14. No Phase 5 work started

**Decision:** stop after Phase 4.

**Why:** the user requested a Phase 4 report before Phase 5. External benchmarks, non-reciprocal evaluation, FIMO, dependence-aware statistics, and final JASPAR decisions remain Phase 5.
