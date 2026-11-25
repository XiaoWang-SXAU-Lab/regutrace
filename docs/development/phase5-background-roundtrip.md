# Phase 5 — Regulatory background export and round-trip

**Date:** 2026-11-25

## 1. What was built

`scripts/phase5/build_background.py` builds a standalone SQLite database from a completed site-only ReguTrace run and its TF–motif mapping. The database does not depend on the ReguTrace service or job ID.

HF1 artifact:

```
%LOCALAPPDATA%\rsWGCNA-build\phase5\background\hf1-background.db
bytes    1,382,932,480
sha256   f515861759b481343177fa6bb9656a315f781a743d9313bb2a0513cfbde4e622
```

## 2. Schema

| Table | Rows | Purpose |
|---|---:|---|
| `tf_catalogue` | 1,981 | all identified + ambiguous TF genes, family/status, motif/target counts |
| `tf_motif` | 32,792 | strict/calibrated/exploratory mapping; coordinates and provenance when available |
| `motif_target` | 4,644,497 | aggregated motif–target associations (n_sites, best p/score) |
| `motif_registry` | 805 | motif name, length, PWM hash, specificity group, source |
| `motif_testability` | 927 | attainable p, requested p, status, hit/mapped counts |
| `promoter` | 37,415 | per-promoter scan/testability/site aggregation |
| `scan_config` | 13 | scan parameters, background, cache and runtime |
| `pairs` | 7,394,435 | policy-selected TF–target prior (best-supported per DBD group) |
| `hits` | 0 | canonical schema required by `validate-db`; site-only runs do not materialize evidence rows |
| `meta` | — | schema v2, policy, counts, checksums, limitations |

The canonical `pairs` table is the rsWGCNA import contract. The normalized testability/promoter/scan layers were added without rebuilding the site table. The 25M exploratory TF–target union is **not** materialized; `motif_target` + `tf_motif` preserve the broader evidence for advanced joins.

## 3. Round-trip steps actually performed

1. Built DB from the HF1 site-only job.
2. `regutrace validate-db --db hf1-background.db` → `valid TF-target prior database` (62.7 s).
3. `regutrace prior --db ... --targets 3 genes` → 321 rows (1.03 s), no service running.
4. Created dashboard project `phase5-roundtrip` (`4bd3e43e046cc62f321b14ce804b37a9`).
5. Imported the 1.3 GB DB as `tf_target` via `POST /api/projects/{id}/inputs/tf_target`; 61.3 s, 200.
6. Queried upstream TFs for 3 genes; 200, 52,465 bytes.
7. Stopped the dashboard, restarted it, reopened the project with `POST /api/projects/open`, repeated the upstream query; 200, 21,014 bytes for one gene.
8. Imported the same DB into a second independent project (`a2f10d5bc8e9cbd370692530a6ee7382`); 200; upstream query returned 19,510 bytes.
9. Compared the imported project file SHA-256 with the exported DB SHA-256: equal.

## 4. Resulting counts

| Item | Value |
|---|---:|
| TF catalogue | 1,981 |
| Identified TFs | 1,907 |
| TFs with pairs (best policy) | 991 |
| Motifs in mapping | 679 |
| Motifs with targets | 614 |
| Motif–target associations | 4,644,497 |
| TF–target pairs | 7,394,435 |
| Target union | 37,415 |
| Pairs' targets outside motif_target | 0 |
| Gene ID matching | 100% (all targets come from the 37,415-promoter universe) |

## 5. Export timing and memory

| Step | Result |
|---|---:|
| DB build (site read + aggregation + pairs + indexes) | ~180 s |
| Observed build working set | ~244 MB |
| Dashboard peak working set during import | ~22 MB (streamed upload; validator child not included) |
| `validate-db` | 62.7 s |
| Project import (upload + validation) | 61.3 s |
| Standalone prior query | 1.03 s |

The build reads the 1.9 GB site TSV once. The database is 1.29 GiB; the WAL and temporary files require additional free space.

## 6. Backward compatibility and failure handling

- `validate-db` passes the canonical `hits`/`pairs` schema and required indexes.
- A SQLite with only a malformed `hits` table fails with an explicit missing-column list.
- The dashboard Data Import accepts `.db` and `.sqlite` and validates before copying.
- TSV/CSV upload remains available through `regutrace convert` for legacy consumers.
- The background does not require genome/GFF/CDS/protein FASTA, promoter rescanning, or the original job directory.

## 7. What remains in-memory

The ReguTrace job registry is still in-memory. After a dashboard restart, old ReguTrace job IDs are not listed; this is why the export is the handoff artifact. The imported rsWGCNA project survives restart and reopens from its project directory, as tested. A future change could persist and lazily reload ReguTrace job metadata, but it is not required for background reuse.

## 8. Limitations

- `pairs` is the best-supported policy (7.39M rows), not the all-candidate union.
- `hits` is empty; site-level evidence remains in the original site table or the `motif_target` aggregate.
- The DB contains evidence level D only; no experimental validation is implied.
- The database is not compressed. A future release could use gzip or a page-compressed SQLite build.
