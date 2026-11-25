# Phase 5 — Regulatory background export and round-trip

**Date:** 2026-11-25

## 1. What was built

`scripts/phase5/build_background_v2.py` builds the final standalone SQLite database in one command: it runs `build_background.py` (canonical v1 hits/pairs), `add_background_layers.py` (testability/promoter/scan layers), `rebuild_motif_registry.py` (registry from the exact scan motif library), and `add_policy_layers.py` (coordinate rerun layer plus exact all/best/reps TF–motif selections), then records the final hash and counts. The database does not depend on the ReguTrace service or job ID.

HF1 artifact:

```
%LOCALAPPDATA%\rsWGCNA-build\phase5\background\hf1-background.db
bytes    1,406,148,608
sha256   7a33d3fc0c00b1b78caf32437842cb02269fba73ddbfe5fb5a44b757cb1bbee4
user_version 3
```

## 2. Schema

| Table | Rows | Purpose |
|---|---:|---|
| `tf_catalogue` | 1,981 | all identified + ambiguous TF genes, family/status, motif/target counts |
| `tf_motif` | 32,792 | strict/calibrated/exploratory mapping; coordinates and provenance when available |
| `tf_motif_coordinate` | 28,465 | exact coordinate rerun layer used by Phase 5 selection |
| `tf_motif_selected` | 1,273 | best-supported motif per TF/DBD group at build time (1,143 mapped TFs, including TFs with no target) |
| `tf_motif_policy` | 44,295 | exact selected TF–motif associations for all/best/reps |
| `tf_motif_policy_config` | 3 | policy version, parameters, source files and SHA-256 |
| `motif_target` | 4,644,497 | aggregated motif–target associations (n_sites, best p/score) |
| `motif_registry` | 927 | motif name, length, PWM hash, specificity group, source |
| `motif_testability` | 927 | attainable p, requested p, status, hit/mapped counts |
| `promoter` | 37,415 | per-promoter scan/testability/site aggregation |
| `scan_config` | 13 | scan parameters, background, cache and runtime |
| `pairs` | 7,394,435 | policy-selected TF–target prior (best-supported per DBD group) |
| `hits` | 0 | canonical schema required by `validate-db`; site-only runs do not materialize evidence rows |
| `meta` | — | schema v3, policy version, counts, checksums, limitations |

The `motif_registry` was rebuilt from the exact 927-row motif library used by the scan, so all 679 mapped motif IDs have registry metadata. The canonical `pairs` table is the rsWGCNA import contract. `tf_motif_selected` records the best-supported selection even for TFs whose selected motifs have no significant target; `pairs` contains only the 991 TFs with at least one target. `tf_motif_policy` stores the exact selection for all three supported policies (32,792 / 1,273 / 10,230 rows), with `policy_version`, parameters and source hashes in `tf_motif_policy_config`. Target sets are reconstructed by joining `tf_motif_policy` to `motif_target`; no per-policy multi-million-row TF–target table is materialized.

## 3. Reproducible build

```bash
python scripts/phase5/build_background_v2.py \
  --coord-map hf1-exploratory-coord-map.tsv \
  --base-map hf1-exploratory-map.tsv \
  --tf-genes tf_genes.tsv \
  --sites <sites.tsv> \
  --testability <motif-testability.tsv> \
  --motifs <job motifs.jaspar> \
  --pwm-hash <motif-pwm-hash.tsv> \
  --specificity-groups <motif-specificity-groups.tsv> \
  --promoter-summary <promoter-summary.tsv> \
  --scan-manifest <scan-manifest.json> \
  --out <background/hf1-background.db> \
  --policy best@0.8
```

The wrapper runs the four builders in the only supported order and writes the final SHA-256 and counts to `background-build.json`. The v3 artifact in this package was produced from the v2 database plus `add_policy_layers.py`; the policy table was verified against `phase5-policy-targets.tsv` and `phase5-policy-summary.json`.

## 4. Round-trip steps actually performed on schema v3

The final v3 database (SHA-256 `7a33…`) was used for every step below; no v1/v2 import is counted for the v3 regression.

1. Built DB from the HF1 site-only job.
2. `regutrace validate-db --db hf1-background.db` → `valid TF-target prior database` (63.0 s).
3. `regutrace prior --db ... --targets 3 genes` → 321 rows (1.03 s), no service running.
4. Created dashboard project `v3-verify` (`a72d242e43af5b60e52fa65a9edb6273`).
5. Imported the 1.4 GB v3 DB as `tf_target` via `POST /api/projects/{id}/inputs/tf_target`; 74.4 s, 200, revision 2.
6. Queried upstream TFs for three target genes through the project API; 200, 753 rows, 2.4 s.
7. Stopped the dashboard, restarted it, reopened the project with `POST /api/projects/open`, repeated the upstream query; 200, 753 rows. The reopened project retained the `7a33…` input hash.
8. The earlier v2 regression (`v2-verify`, `v2-reuse-2`) remains valid for the v2 schema; v3 adds only compact policy/coordinate tables and was re-imported as above.
9. Compared the imported project file and the exported DB with SHA-256: both `7a33d3fc0c00b1b78caf32437842cb02269fba73ddbfe5fb5a44b757cb1bbee4`.

Query genes: `Fta_6G0020080`, `Fta_1G0002980`, `Fta_2G0028680`. All three are present in `promoter`, all returned upstream TFs are in `pairs`, and no `pairs.target_id` falls outside the 37,415-promoter universe.

## 5. Resulting counts

| Item | Value |
|---|---:|
| TF catalogue | 1,981 |
| Identified TFs | 1,907 |
| TFs with pairs (best policy) | 991 |
| Motifs in mapping | 679 |
| Motifs with targets | 614 |
| Best-supported selected motifs | 1,273 rows / 1,143 mapped TFs |
| Coordinate rerun rows | 28,465 |
| Exact policy rows (all / best / reps) | 32,792 / 1,273 / 10,230 |
| Motif–target associations | 4,644,497 |
| TF–target pairs | 7,394,435 |
| Target union | 37,415 |
| Pairs' targets outside `promoter` | 0 |
| Mapped motifs missing `motif_registry` | 0 |
| Gene ID matching | 100% (all targets come from the 37,415-promoter universe) |

## 6. Export timing and memory

| Step | Result |
|---|---:|
| v2 DB build (site read + aggregation + pairs + indexes) | ~196 s |
| v3 policy/coordinate layers | 14.8 s including vacuum |
| Observed build working set | ~244 MB |
| Dashboard peak working set during import | ~22 MB (streamed upload; validator child not included) |
| `validate-db` | 63 s (v2); v3 PASS |
| Project import (upload + validation) | 74.4 s (v3) |
| Standalone prior query (3 genes) | 1.03 s |
| Project API upstream query (3 genes) | 2.4 s |
| Dashboard restart + reopen + repeat query | PASS, 753 rows |

The v2 build reads the 1.9 GB site TSV once. The final v3 database is 1.31 GiB; the WAL and temporary files require additional free space.

## 7. Backward compatibility and failure handling

- `validate-db` passes the canonical `hits`/`pairs` schema and required indexes.
- A SQLite with only a malformed `hits` table fails with an explicit missing-column list.
- The dashboard Data Import accepts `.db` and `.sqlite` and validates before copying.
- TSV/CSV upload remains available through `regutrace convert` for legacy consumers.
- The background does not require genome/GFF/CDS/protein FASTA, promoter rescanning, or the original job directory.

## 8. Downstream access to alternative motif policies

All three policies are exact and require only joins to `tf_motif_policy`:

| View | Source tables | Verified exact result |
|---|---|---|
| All candidates | `tf_motif_policy` (policy=`all`) ⨝ `motif_target` | 1,124 TFs / 25,481,360 pairs; 0 per-TF target-set mismatches |
| Best-supported | `tf_motif_policy` (policy=`best`) ⨝ `motif_target`; also `pairs` | 991 TFs / 7,394,435 pairs; 0 mismatches vs `pairs` |
| Representative specificity groups | `tf_motif_policy` (policy=`reps`) ⨝ `motif_target` | 1,104 TFs / 19,848,826 pairs; 0 per-TF target-set mismatches |

Verification is recorded in `policies/phase5-policy-verification.json`: selection TF/motif sets are compared exactly with `phase5-policy-targets.tsv`, and target sets are compared per TF as sets (not only aggregate counts). Pair-set SHA-256: all `97315d02…`, best `4bba3330…`, reps `3fb1cbd0…`.

Minimal exact query:

```sql
SELECT p.tf_id, mt.target_id, MIN(mt.best_p) AS best_p,
       COUNT(DISTINCT p.motif_id) AS n_motifs
FROM tf_motif_policy p
JOIN motif_target mt ON mt.motif_id = p.motif_id
WHERE p.policy = 'reps'
GROUP BY p.tf_id, mt.target_id;
```

`tf_motif_coordinate` preserves all 28,465 coordinate rerun rows for provenance and recomputation. `tf_motif_policy_config` stores policy version `phase5-v1`, parameters, and SHA-256 for the policy targets and coordinate sources. `motif_registry.specificity_group_080` supplies the representative key; `motif_target` has indexes on both `motif_id` and `target_id`. The policy tables are compact (44,295 rows total); no policy-specific TF–target table is materialized.

Phase 6 should consume the prior through this documented interface; it should not add WGCNA-specific filtering or ranking inside ReguTrace.

## 9. What remains in-memory

The ReguTrace job registry is still in-memory. After a dashboard restart, old ReguTrace job IDs are not listed; this is why the export is the handoff artifact. The imported rsWGCNA project survives restart and reopens from its project directory, as tested. A future change could persist and lazily reload ReguTrace job metadata, but it is not required for background reuse.

## 10. Limitations

- `pairs` is the best-supported policy (7.39M rows), not the all-candidate union.
- `hits` is empty; site-level evidence remains in the original site table or the `motif_target` aggregate.
- The DB contains evidence level D only; no experimental validation is implied.
- The database is not compressed. A future release could use gzip or a page-compressed SQLite build.
- The v3 round trip was tested through the HTTP API and headless dashboard, not through the desktop WebView UI.
- Coordinate status is **1,089 successfully joined + 2 key mismatches + 52 unavailable**. The 52 unavailable base TFs use the documented opaque-group fallback in the original Phase 5 selection; the v3 exact policy table preserves that selection. The DB cannot add coordinates that were absent.
- The standard Phase 4 job `tf_map.tsv` has no `query_full_dbd_*` columns. The dashboard policy UI therefore reports `coordinate_available=false` and uses the opaque fallback for current jobs; the coordinate rerun exists only in the Phase 5 build artifacts.
- Representative groups are stored as exact TF–motif associations, not as a TF–target table; downstream must perform the join. The current Phase 4 `tf_map.tsv` still lacks coordinate columns, so the live UI uses its documented fallback even though the background itself is exact.
