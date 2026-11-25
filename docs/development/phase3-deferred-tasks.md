# Phase 3 deferred tasks

Date: 2026-11-25

These items are explicitly **not completed**. They must not be reported as
Phase 3 results or used to justify a promotion.

## Statistical

- Site-level and promoter-level q-values with a dependence-aware null.
- A sequence-based promoter null preserving length, GC and dinucleotide
  composition (for example, dinucleotide-preserving shuffles or a Markov
  background).
- Parameter-uncertainty correction for the empirical background, which is
  estimated from the same promoters that are scanned.
- Block/bootstrap calibration at the promoter or gene level; pair-level or
  site-level resampling is pseudoreplicated.
- Motif-family and shared-PWM dependence in multiple testing.
- A calibrated motif enrichment model with length/GC matching and gene-set
  overlap handling.
- FIMO comparison. FIMO is not installed locally; no agreement claim is made.

## Data model and operations

- Normalized SQLite tables/views for scan config, site, promoter aggregation and
  testability; the current DB remains the legacy `hits`/`pairs` projection.
- Cache eviction, size limits and cross-version migration.
- Genomic coordinates for `--promoters` FASTA input without a GFF. Current site
  rows correctly leave genomic coordinates empty in that mode.
- A stable scan-manifest schema version and migration policy.
- Dashboard integration for site/testability/manifest panels.

## Mapping and validation (Phase 5)

- External Nature Plants 2025 cross-species binding dataset audit.
- Plant Physiology 2025 conserved motif dataset audit.
- Independent DAP-seq/ChIP-seq target benchmark.
- Family/subfamily cis-element validation.
- Non-reciprocal DBD candidate strategy evaluation under species- and
  homology-cluster-held-out validation.
- Incremental mapping precision among TFs not recovered by strict RBH.
- Combined/cascading mapping strategies (strict first, calibrated/exploratory
  fallback) with explicit evidence tiers.
- Final JASPAR 2026 default pack decision; JASPAR 2024 remains the default,
  JASPAR 2026 remains a candidate.
- Phase 2.7 threshold refinement; not started, and not to be tuned on the Phase
  2.5 calibration set without new external data.

## Scan performance

- Multi-motif batching per promoter or SIMD/GPU acceleration; Phase 3 stays
  CPU-only and offline.
- Threshold-free scan cache storing complete score distributions (would allow
  removing p-value from the cache key).
- Parallel group-level scanning with deterministic merge if memory permits.
- Approximate score distributions or branch-and-bound scanning; not implemented
  because exactness and p-value auditability were prioritized.

## Enrichment

- GMT collections for HF1; none were supplied, so the interface was tested
  synthetically.
- Motif aggregation by exact PWM hash or specificity group with explicit
  reporting.
- Expression, chromatin and orthology evidence integration; belongs to Phase 6–7.
