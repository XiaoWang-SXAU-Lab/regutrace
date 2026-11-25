# ReguTrace Phase 5 — Scientific validation and reusable background finalization

**Date:** 2026-11-25
**Status:** Phase 5 implementation and validation completed for the stated scope. Phase 6 not started.
**Scope:** standalone ReguTrace; no expression input, no WGCNA-specific filtering, no new DBD threshold optimization.

Supporting records:

- [phase5-background-roundtrip.md](phase5-background-roundtrip.md)
- [phase5-motif-aggregation.md](phase5-motif-aggregation.md)
- [phase5-target-benchmark.md](phase5-target-benchmark.md)
- [phase5-statistics.md](phase5-statistics.md)
- [phase5-decisions.md](phase5-decisions.md)

Machine-readable results are under `%LOCALAPPDATA%\rsWGCNA-build\phase5\` and in the sync package described below.

## 1. Headline results

### HF1 denominator reconciliation

The 1,908 vs 1,907 difference is **one gene**:

| Catalogue | Identified | Ambiguous | Difference |
|---|---:|---:|---|
| CDS-derived `hf1-pipeline-2026/tf_genes.tsv` | 1,908 | 73 | `Fta_6G0009940` identified AP2 |
| PEP-derived current catalogue | 1,907 | 74 | `Fta_6G0009940` ambiguous AP2;ERF |

No other identified gene differs. The denominator used for Phase 5 remains 1,907 **because it comes from the PEP TF identification catalogue, not from motif-mapping rows**. The CDS run's 1,908 is a valid alternative catalogue, not a lost target. See `denominator-reconciliation.json`.

### HF1 motif policy comparison

All values use 1,907 identified TFs; 1,143 have a candidate motif in exploratory mode. Target counts are evidence level D.

| Policy | Candidate motifs | ≥1 target | Coverage | Median motifs/TF | Median targets/TF | Target Q1/Q3 | Total TF–target pairs |
|---|---:|---:|---:|---:|---:|---:|---:|
| All candidate motifs | 1,143 | 1,124 | 58.94% | 17 | 21,715 | 13,891 / 34,348 | 25,481,360 |
| Best-supported per DBD group | 1,143 | 991 | 51.97% | 1 | 4,091 | 2,290 / 8,385 | 7,394,435 |
| Representative per specificity group | 1,143 | 1,104 | 57.89% | 6 | 14,657 | 7,644 / 28,711 | 19,848,826 |

The best-supported policy cuts the candidate pair total by **71.0%** (median per-TF target count to 25.2% of the union) while retaining 51.97% actual coverage. The representative policy keeps alternative specificity hypotheses and retains 57.89% coverage with a 22.1% pair reduction.

### Coordinate-aware audit

The hypothesis that dozens of motifs come from different query DBD regions is **false for HF1**. The coordinate rerun (28,465 mapping rows with query/donor coordinates) shows:

- every coordinate-covered DBD-mapped TF has **one biological DBD region group** at overlap thresholds 0.5, 0.7, 0.8 and 0.9;
- all 19 `Fta_1G0005130` motifs align to the same query interval 152–198 (bHLH), with 19 distinct PWMs but only 7 specificity groups;
- the inflation is therefore **many donor-derived motif hypotheses for one DBD region**, not multiple independent DBDs;
- target overlap between those 19 motifs is low: median pairwise Jaccard 0.0018 (some pairs are identical, max 1.0, but most target sets differ);
- the all-candidate union is 8,709 targets for `Fta_1G0005130`; best-supported is 2,290; representative groups is 4,112.

This directly answers the addendum: grouping by DBD coordinates does not reduce motif count in HF1 because the redundancy is donor/PWM redundancy, not coordinate redundancy. Motif redundancy removal and specificity grouping are the effective levers.

### Multiple testing

Using the complete tested-window universe (37,415 promoters, both strands, 1,000 bp) gives very different target sets:

| Policy | Original p≤1e-4 | Site-level BH q≤0.05 | Promoter-level BH q≤0.05 |
|---|---:|---:|---:|
| All candidate motifs | 1,124 TFs / 58.94% | 429 / 22.50% | 417 / 21.87% |
| Best-supported | 991 / 51.97% | 86 / 4.51% | 81 / 4.25% |
| Representative groups | 1,104 / 57.89% | 326 / 17.09% | 287 / 15.05% |

At p=1e-4, a 10 bp motif is tested in ~74 million windows; ~7,400 false positives are expected per motif before correction. The raw cutoff is therefore a **candidate prefilter**, not a genome-wide significance claim. Site/promoter FDR is reported as an optional interpretation layer; the scanner default is unchanged.

A dinucleotide-preserving shuffle of 200 promoters (3 benchmark motifs, 200 kb) produced **61.5% promoters with a hit**, versus **58.5% observed** in the same genes. At this scale, raw promoter hit presence is not enriched over composition-matched sequence; this reinforces the need for promoter-level null/FDR and against treating hit counts as evidence. Full details in [phase5-statistics.md](phase5-statistics.md).

### Independent DAP-seq benchmark

O'Malley et al. 2016 Arabidopsis DAP-seq peaks (GEO GSE60141) were used, with TAIR10 genome/annotation and 1,000 bp promoters. Positive genes are promoters overlapping a peak; unpeaked genes are unknown, not negatives.

| Setting | TFs evaluated | Median AUPRC | Median enrichment | Median recall | Median PPV_assay |
|---|---:|---:|---:|---:|---:|
| Direct motif, uniform background | 7 | 0.330 | 4.89× | 0.489 | 0.376 |
| Direct motif, empirical background | 7 | 0.346 | 3.20× | 0.685 | 0.236 |
| Held-out DBD-transferred best motif | 7 | 0.261 | 2.81× | 0.335 | 0.255 |

On the 6 TFs common to direct and transferred settings, direct motif AUPRC median 0.314 vs transfer 0.279; direct enrichment 4.61× vs 3.75×. Direct motif lookup is better on average, but transfer still recovers substantial signal for several TFs. `TGA1` transferred better than direct (AUPRC 0.449 vs 0.255), showing that donor choice matters.

Limitations: promoter overlap only (distal peaks ignored); unobserved genes treated as non-positive for AUPRC/PPV, not as experimental negatives; DAP-seq binding is not functional regulation; one assay series; 7 TFs, not a genome-wide benchmark; RAV1 direct motif was unattainable at p=1e-4; GBF3 had no supported DBD transfer; PDF2 lacked a direct JASPAR motif in the benchmark map.

### Incremental DBD quality

- Strict 271 TFs, balanced 280, exploratory 1,143.
- Exploratory adds **872 TFs over strict**, **863 over balanced**; 857 of the 872 have at least one target.
- 11 ambiguous TFs (29 rows) are visible but excluded from best/representative target policies by default.
- The additional TFs are concentrated in ERF (91), MYB (145), bHLH (113), bZIP (69), WRKY (62), NAC (51), HD-ZIP (38), Dof (37), LBD (34).
- Phase 2.6 donor-held-out method D precision was 0.922 (95% CI low 0.890); homology-cluster-held-out method D precision was 0.727 (CI low 0.652). These are **method- and setting-specific** estimates and are **not** applied to all 1,143 exploratory HF1 TFs. No universal threshold was introduced.

Full table in `incremental-quality.json`.

## 2. Task A — export/import round trip

A standalone SQLite background was built from the HF1 site-only job and the coordinate-aware mapping:

| Item | Value |
|---|---:|
| File | `hf1-background.db` |
| Size | 1,382,932,480 bytes (1.29 GiB) |
| SHA-256 | `f515861759b481343177fa6bb9656a315f781a743d9313bb2a0513cfbde4e622` |
| TF catalogue | 1,981 rows (1,907 identified + 74 ambiguous) |
| TF–motif rows | 32,792 (27,720 with query DBD coordinates; 1,089 TFs) |
| Motif–target associations | 4,644,497 |
| TF–target pairs (best-supported policy) | 7,394,435 |
| Distinct TF with pairs | 991 |
| Distinct motifs in mapping / with targets | 679 / 614 |
| Target union | 37,415 |
| Normalized layers | motif_testability 927; promoter 37,415; scan_config 13 |
| Build wall time | ~180 s |
| Observed build working set | ~244 MB |
| Dashboard peak working set during import | ~22 MB (validator child excluded) |
| `validate-db` time | 62.7 s |
| Project Data Import time | 61.3 s |
| Standalone prior query (3 genes) | 1.03 s |

Round trip actually performed:

1. Built the SQLite with `scripts/phase5/build_background.py`.
2. `regutrace validate-db` passed.
3. `regutrace prior` queried 3 targets and returned 321 upstream rows without the ReguTrace service.
4. Created a dashboard project, imported the 1.3 GB `.db` as `tf_target`; import succeeded.
5. Queried upstream TFs through the project API (200; 52 KB).
6. Stopped and restarted the dashboard, reopened the project, and repeated the upstream query successfully (200). This verifies the background remains usable after application restart.
7. The imported project file's SHA-256 equals the exported DB SHA-256.

The background is independent of the original service, job ID, genome/GFF/CDS/protein uploads, and re-scanning. The job registry itself is still in-memory; the exported DB and project state are the reusable artifacts. The project must be opened after a dashboard restart, but it no longer depends on the ReguTrace job.

Schema: `tf_catalogue`, `tf_motif` (with coordinates when available), `motif_target`, `motif_registry`, canonical `hits` (empty by design for site-only), and `pairs` (best-supported policy). `validate-db` requires `hits`/`pairs`, so the DB is importable; `motif_target` preserves the broader 4.6M motif-level layer. The full 25M exploratory TF–target union is deliberately not materialized.

See [phase5-background-roundtrip.md](phase5-background-roundtrip.md).

## 3. Task B — why exploratory target sets are large

Median exploratory TF: **17 motifs, 17 effective PWMs, 6 specificity groups, 21,715 targets**.

Evidence that the union of alternative donor motifs dominates:

- motif count vs target count Spearman correlation **0.704**;
- best-supported reduces median targets to 4,091 and total pairs to 25.2% of the union;
- representative specificity groups reduce pairs to 77.9% of the union while keeping 6 motifs;
- the regression TF `Fta_1G0005130`: 19 motifs, 7 specificity groups, 8,709 targets; best motif 2,290; representative groups 4,112; median pairwise motif Jaccard 0.0018.

The 19 motifs are not independent experimental evidence. They are donor/PWM hypotheses for the same biological bHLH DBD region. The low target overlap means the union grows almost linearly with the number of donor motifs.

Five biological/technical causes are distinguishable:

1. 615 donors / 679 motif IDs cover many plant families, with several donors per family.
2. Many donor PWMs are similar but not identical, producing different promoter hits.
3. Short motifs have low information content and high hit rates.
4. At p=1e-4 across ~74M windows per motif, random hits are expected.
5. Overlapping/shared PWMs let several TFs explain the same site.

The correct response is not arbitrary target removal. The background preserves all candidate evidence; best/representative policies are analysis and export options.

## 4. Addendum 1 — coordinate-aware selection

Implementation: `scripts/phase5/aggregate_audit.py`.

- Query DBD coordinates come from a coordinated rerun of `map-dbd-motifs` with full alignment columns (`query_full_dbd_start/end`, donor coordinates, identity/coverage, architecture).
- Grouping is by compatible extraction family/domain mode and query interval overlap; composite architectures (RAV) are not split.
- Exact PWM duplicates are removed by PWM hash (specificity-group table).
- Representative policies:
  - **best@0.8** — one best-supported motif per DBD group;
  - **reps@0.8** — one best motif per specificity group (group_080) within each DBD group;
  - **all** — unchanged union.
- Ranking uses architecture compatibility, query/donor coverage, identity, bitscore, and deterministic tie-breaking; it never uses target count or promoter p-value.
- Overlap thresholds 0.5/0.7/0.8/0.9 produced the same HF1 groups because every TF's candidate alignments overlap one query region. The threshold is retained as a configurable safeguard for future multi-domain cases.
- 14 strict/calibrated mappings in the layered map are preserved unchanged; only exploratory/ambiguous rows are subject to grouping. Best/reps exclude ambiguous rows from target sets but the original evidence remains in `all`.

## 5. Addendum 2 — multiple testing and promoter-level statistics

Implementation: same audit script, using exact score p-values from the Phase 3 scanner and the actual tested-window count.

- **FIMO audit:** FIMO computes exact score p-values per motif occurrence and BH q-values; both strands are scored by default; q-values are approximate if the stored-score cap is reached. FIMO's q-value is a sequence-level statistic, not a biological interaction FDR. The local scanner's p-values are exact for its discrete score/zero-order background, but the default 1e-4 is not a genome-wide FDR.
- **Site-level BH:** number of tests per motif = `37,415 × 2 × (1,000 − motif_length + 1)` (all promoters, both strands). Rank among hits is exact for p≤cutoff; non-hit windows have larger p and cannot reduce a hit's q, so the hit-only implementation is a conservative upper bound. This avoids ordinary BH on the prefiltered table alone.
- **Promoter-level:** best site p per promoter, Bonferroni over the actual windows in that promoter (`p_best × 2 × (1,000−m+1)`), then BH across 37,415 promoters. This is conservative and ignores dependence between overlapping windows; it is not a biological FDR.
- **Empirical null:** 200 promoters shuffled with exact dinucleotide preservation; 3 motifs; 198 hits in shuffled sequences vs 117/200 observed promoters with a hit. Shuffled hit rate 61.5% vs observed 58.5%; median best p 4.51e-5 vs 3.05e-5. Raw promoter hit presence is not enriched over composition-matched sequence at this scale.

Important: site/promoter FDR collapses coverage at q≤0.05 because the testing universe is enormous. This is a real statistical property, not a bug; the default cutoff remains p=1e-4 and original p-values are preserved. Adjusted fields are optional analysis layers, not replacements.

## 6. Task D and final policy

Existing evidence tiers are preserved:

- strict only;
- strict + calibrated DBD;
- extended inferred motifs;
- selected exploratory candidates;
- complete evidence-aware background.

All PWM associations remain level D. No expression-based filtering was added. No threshold was re-optimized using the DAP benchmark.

## 7. Tests and reproducibility

Actually run:

| Check | Result |
|---|---|
| Rust `cargo test --workspace` | 124 passed (46 CLI unit + 7 CLI + 2 Phase 3 + 69 core) |
| Python Phase 3 tests | 6 passed |
| Phase 5 unit tests (`scripts/phase5/test_phase5.py`) | 5 passed |
| Go `internal/regutrace` | 46 passed |
| Go `internal/api` | 13 passed |
| `go vet` changed packages | pass |
| `vue-tsc -b && vite build` | pass |
| `regutrace validate-db` on HF1 background | pass |
| Dashboard project import of 1.3 GB background | pass |
| Project upstream query after dashboard restart | pass |
| Invalid schema rejection | pass (missing columns reported) |

New Phase 5 scripts are reproducible from the documented inputs. The analysis does not require re-scanning promoters; it reuses the Phase 3 site table and cache.

## 8. Limitations

- DAP-seq benchmark is promoter-overlap only, 7 TFs, one assay series; distal peaks and unknown negatives are not resolved.
- No FIMO execution comparison was performed (FIMO is not installed); the audit is documentation/source-level plus our own BH implementation.
- The HF1 background DB uses best-supported pairs for the importable prior; the full exploratory union is represented only by `tf_motif`/`motif_target`, not by 25M `pairs`.
- The dashboard job registry remains in-memory; exports and projects survive restart, but old ReguTrace job IDs do not.
- The empirical null is 200 promoters, 3 motifs; it is illustrative, not a genome-wide null model.
- No external ChIP-seq cross-species benchmark beyond the O'Malley DAP-seq subset.
- No commit/push/installer rebuild was performed.

## 9. Deliverables

- this report and the four supporting Phase 5 documents;
- machine-readable JSON/TSV in `%LOCALAPPDATA%\rsWGCNA-build\phase5\`;
- HF1 SQLite background (local; sync copy in `analysis/regutrace-hf1-phase5-background-2026-11-25/`);
- scripts:
  - `scripts/phase5/aggregate_audit.py`
  - `scripts/phase5/build_background.py`
  - `scripts/phase5/add_background_layers.py`
  - `scripts/phase5/benchmark_dap.py`
  - `scripts/phase5/null_shuffle.py`

**Stop point:** Phase 5 ends here. Phase 6 has not been started.
