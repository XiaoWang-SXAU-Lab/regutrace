# Phase 4 query architecture

**Date:** 2026-11-25
**Scope:** dashboard backend (`wgcna-dashboard/backend/internal/regutrace`) and the Phase 3 site table.

## 1. Two dataset modes

The backend keeps the existing evidence-row dataset for ordinary runs and adds a site-only dataset for Phase 3 `--sites-only` runs.

### Evidence mode (existing)

- `evidence.tsv` is indexed by TF and target with on-disk byte offsets.
- `Query`, `PromoterView`, `TFView`, and export operate on evidence rows.
- Memory is proportional to the number of TF/target keys, not to all evidence rows.

### Site mode (new)

`buildSiteDataset` scans the unique site TSV once and builds:

- `targetOffsets[target_id] -> []byte offset` for each site row;
- `motifOffsets[motif_id] -> []byte offset` for each site row;
- `motifTargets[motif_id] -> bitset` of target indexes;
- `tfTargetCounts[tf_id]` by OR-ing the bitsets of the TF's mapped motifs;
- summary counts: unique sites, motif–target pairs, TF–target pairs, unique targets.

It never creates TF × target × site rows. A site row stores `motif_ids` as a semicolon list for PWM groups that share a matrix; TF identity is resolved from the mapping table at query time.

## 2. Query paths

### TF → targets

`siteTFTargets(tfID)`:

1. resolve the TF ID (case-insensitive, catalogue fallback);
2. collect the TF's mapped motif IDs, excluding single-domain fallback rows unless enabled;
3. read only the site offsets for those motifs;
4. aggregate per target: best p/score, site count (distinct hit coordinates), motif IDs, mapping statuses, shared-PWM count;
5. filter by mapping mode, explicit status, and max p;
6. sort by best p, then target ID;
7. paginate.

The aggregation result is cached per job for up to 4 `(tf, status, max_p, mapping_mode)` keys. The cache is bounded; eviction is FIFO.

### Target → regulators

`siteTargetRegulators(target)`:

1. resolve the target through the promoter catalogue/site index;
2. read only the target's site offsets;
3. for every motif at each site, expand to mapped TFs;
4. aggregate per TF: best p, distinct sites, motifs, statuses, shared-PWM count;
5. sort by best p, then TF ID; paginate.

This path is naturally bounded by the number of sites in one promoter (typically hundreds), so no cache is used.

### Promoter view

`sitePromoterHits(target)` groups sites by `(motif_ids, hit_start, hit_end, hit_strand)`, attaches mapped TFs, sorts by position, and adds sequence slices. The promoter sequence itself is loaded from the promoter FASTA/Promoter map.

### Exact site evidence

`siteEvidenceFromSites` reads the target's offsets and matches motif, hit coordinates, and strand. It returns one row per mapped TF, preserving shared-PWM ties.

### Motif list/detail

Motif metadata and testability are loaded from the motif library and Phase 3 testability table. Detail adds:

- normalized matrix and PWM hash;
- registry entries (donor/source/method/species/license);
- mapped TFs and specificity groups;
- target count from the site dataset (site mode) or evidence rows (evidence mode).

## 3. Pagination and limits

| Endpoint | Default | Maximum |
|---|---:|---:|
| `/targets` | 50 | 500 |
| `/regulators` | 50 | 500 |
| `/query` | 200 | 1,000 |
| `/motifs` | 50 | 200 |
| enrichment motifs | — | 200 |

All offsets are non-negative; over-range offsets clamp to the total. Rows are deterministically ordered, so pagination is stable.

## 4. Mapping-mode filtering

`modeAllows`:

- `strict` → `strict_rbh`
- `balanced` → `strict_rbh`, `dbd_calibrated`
- `exploratory` → all statuses, including `dbd_ambiguous` and `dbd_exploratory`

The UI sends no mode parameter: the job's configured mode is authoritative. Per-row status filters can only narrow within the mode. This prevents a user from silently viewing exploratory rows in a strict job.

## 5. Stale and concurrent requests

- The main TF/promoter/site-evidence views use monotonic sequence counters and discard older responses.
- New Phase 4 components use the same pattern for target, regulator, motif, enrichment, and provenance loads.
- HTTP requests themselves are not aborted; stale responses are discarded client-side.
- The backend serves reads from immutable completed-job artifacts; no partial generation is exposed.

## 6. Memory characteristics

The site dataset index is the dominant memory cost:

- 8,191,211 target offsets and 8,191,211 motif offsets (int64 each);
- per-key slices and map overhead;
- 677 motif target bitsets (about 37,415 bits each, ~3.2 MB total);
- TF target counts and summary structures.

Measured server working set after indexing and heavy exploratory queries reached ~2.0 GB; after light use it was ~0.6–0.9 GB. See [phase4-performance.md](phase4-performance.md). The index is rebuilt when a job is loaded after a server restart; it is not persisted as a binary index.

## 7. What was intentionally not done

- No TF × target × site materialization.
- No loading of the site TSV into browser memory.
- No client-side full-table sorting.
- No silent p-value relaxation.
- No automatic promotion of exploratory mappings.
- No WGCNA-specific filtering or expression-based ranking inside ReguTrace.

## 8. Known bottlenecks

- First TF target aggregation reads all sites for the TF's motifs; 36,513 targets took ~6.3 s.
- The legacy free-query endpoint still builds all matching rows before pagination; a large TF query took ~4.2 s and returned 483 kB for 200 rows.
- Promoter views can be ~0.5 MB per target because they include sequence and every site/TF mapping.
- The job summary embeds the full catalogue (~1.08 MB).
