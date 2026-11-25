# Phase 3 scan and evidence data model

Date: 2026-11-25  
Status: implemented for scan/site/testability layers; normalized SQLite
projection deferred

## Goal

Separate objects that were adjacent in the legacy evidence table:

1. TF–motif mapping;
2. motif registry;
3. unique PWM scan groups;
4. motif–promoter site hits;
5. motif–target aggregation;
6. TF–target candidate associations.

The separation prevents a shared PWM site from being counted as independent
evidence for every TF that carries the motif, and allows exact scan reuse without
merging similar-but-different matrices.

## Implementation map

| layer | implementation | status |
|---|---|---|
| TF–motif mapping | `tf_map.tsv`, mapping diagnostics (Phase 0–2.5) | existing |
| motif registry | `registry.tsv`, `Pwm::digest` | existing + Phase 3 hash |
| unique PWM groups | in-memory `PwmGroup`, `scan-manifest.json` | new |
| motif–promoter site | `motif-site-hits.tsv` (`--sites-out`) | new |
| promoter aggregation | `promoter-summary.tsv` (`--promoter-summary-out`) | new |
| TF–target association | legacy `evidence.tsv` / SQLite `hits`+`pairs` | compatibility projection |
| testability | `motif-testability.tsv` (`--motif-testability-out`) | new |
| enrichment | `regutrace enrich` over sites + gene sets | new, exploratory |

The SQLite database remains the legacy `hits`/`pairs` compatibility projection.
Normalized SQLite tables for site/scan/aggregation are deferred to Phase 5 or
later so the dashboard and `prior` exporter keep working unchanged.

## Layer 1 — TF–motif mapping

`tf_map.tsv` is unchanged. Mapping status remains one of:

```text
strict_rbh
dbd_calibrated
dbd_ambiguous
dbd_exploratory
no_supported_motif
```

A `dbd_calibrated` mapping still produces evidence level D when a promoter site
is found. Mapping confidence never becomes binding evidence.

## Layer 2 — motif registry and exact PWM identity

`Pwm::digest()` hashes:

- a format string (`regutrace-pwm-matrix-v1`);
- the number of positions;
- each normalized probability as `f64::to_bits()` in A,C,G,T order.

The motif accession and name are excluded. Two accessions with a bit-identical
normalized matrix share a digest; matrices differing by one floating-point bit do
not. There is no similarity threshold in this layer.

`motif-testability.tsv` also exposes the digest, matrix length, information
content and maximum score per motif.

## Layer 3 — unique PWM groups

`--reuse-exact-pwm` groups motif IDs by digest. One representative matrix is
scanned per group. The group keeps the sorted list of motif IDs, so evidence
rows are still emitted for every TF–motif mapping.

`scan-manifest.json` records:

- mapped motif IDs;
- unique PWMs;
- PWM groups;
- duplicate PWM groups (groups with more than one motif ID);
- whether exact reuse was enabled.

The cache key is a SHA-256 over:

```text
algorithm version
score scale (10)
background A,C,G,T bit patterns
requested p-value
min-score (or its absence)
strand mode
promoter-set hash (ordered gene IDs + sequences)
PWM digest
```

The p-value cutoff remains in the key because the scanner stores only hits, not
the complete score distribution. A future threshold-free cache could remove it.

## Layer 4 — motif–promoter site hits

`motif-site-hits.tsv` columns:

```text
pwm_hash
motif_ids                 semicolon-separated exact group
gene_id
chrom
promoter_start
promoter_end
promoter_strand
hit_start                 oriented-promoter coordinate, 1-based inclusive
hit_end
hit_strand
genomic_hit_start         empty when only --promoters FASTA was supplied
genomic_hit_end
genomic_hit_strand
score                     log2-odds
p_value
tf_count                  distinct TFs mapped to the group
scan_key                  persistent-cache key
```

One row is written per unique PWM group and site. If two motif IDs share the
matrix, they share one site row and `motif_ids` contains both. The site row is
not duplicated per TF.

### Coordinate convention

Promoter `seq` is oriented 5′→3′ in the gene direction. Hit coordinates refer to
that oriented sequence. For a minus-strand gene, genomic coordinates are:

```text
genomic_start = promoter_end - hit_end + 1
genomic_end   = promoter_end - hit_start + 1
genomic_strand = opposite(hit_strand)
```

The unit test `genomic_coordinates_follow_promoter_orientation` covers both
orientations. HF1 JASPAR 2026 had 0 site rows missing genomic coordinates when
genome/GFF input was used.

## Layer 5 — promoter aggregation

`promoter-summary.tsv` columns:

```text
gene_id, chrom, promoter_start, promoter_end, promoter_strand,
promoter_length, scanned_pwms, testable_pwms, pwms_with_hits,
sites, best_pvalue, best_score, best_pwm_hash, best_motif_ids
```

`best_pvalue` is the smallest **site** p-value in that promoter. It is a
descriptive statistic, not a promoter-level p-value. Overlapping windows and
multiple motifs make the minimum non-uniform under the null; no promoter q-value
is reported.

## Layer 6 — TF–target association

The legacy evidence projection remains:

```text
one row per TF × motif ID × site
evidence_type = prediction_motif
evidence_level = D
```

This preserves downstream `pairs`, `prior` and dashboard behavior. The site layer
is the correct place to count unique sites. The evidence layer may contain
several rows for one site when several TFs share the motif; those rows are not
independent biological evidence.

HF1 checks:

- Run 1 (legacy) vs Run 2 (exact reuse + streaming): full-row set difference
  including notes = 0 / 0.
- Run 3 optimized vs the pre-Phase-3 JASPAR 2026 pipeline DB: 0 / 0.
- Run 5 strict subset vs Run 3: 0 / 0.

## Testability and attainable p-values

`motif-testability.tsv` columns:

```text
motif_id, pwm_hash, matrix_length, information_content, max_score,
min_pvalue, requested_pvalue, attainable, status, mapped,
motif_id_hit_count, tf_count, background_mode
```

`status` is one of:

```text
significant_hits
no_significant_hits
threshold_unattainable
not_scanned
```

`invalid_motif` is not emitted because `read_jaspar` rejects malformed or
duplicate-ID motifs before scanning. The requested threshold is always preserved;
the scanner never silently lowers it. Motifs whose exact floor is above the
requested p-value are skipped and reported as `threshold_unattainable`.

## Empirical background

`--background-mode empirical` estimates a strand-symmetric zero-order model from
the supplied promoters:

```text
A = T = (count_A + count_T) / 2
C = G = (count_C + count_G) / 2
p_base = (symmetrized_count + 1) / (2*(A+C) + 4)
```

Non-ACGT bases are counted as ambiguous, not as background. The scan still skips
windows containing them. The manifest records raw counts, ambiguous bases,
pseudocount and symmetrization.

The model is estimated from the same promoters that are scanned. P-values are
conditional on that estimate; no parameter-uncertainty correction is applied.
This is a zero-order model, not a dinucleotide or GC-window model.

## Enrichment interface

`regutrace enrich` reads a site table, a promoter summary and a GMT. It builds a
binary motif × gene hit matrix over the scanned promoter universe and runs a
one-sided Fisher exact test per motif and gene set, with BH q-values across
motifs within each gene set.

The interface is implemented and tested, but the test is exploratory: it does
not match promoter length, GC, dinucleotide composition, gene-set overlap or
motif similarity. See [phase3-enrichment-design.md](phase3-enrichment-design.md).

## Acceptance tests

Implemented:

- exact PWM reuse: one scan for two identical matrices;
- similar but non-identical matrices have different digests and separate scans;
- changed p-value, background or promoter set invalidates the cache;
- shared motif site appears once in the site table and once per TF in evidence;
- plus/minus genomic coordinate conversion;
- testability status for an unattainable short motif;
- legacy vs optimized full-row equality on HF1;
- enrichment Fisher/BH and CLI integration;
- UTF-8/ASCII deterministic output.

Deferred:

- normalized SQLite tables/views for scan/site/aggregation;
- dependence-aware site/promoter q-values;
- cache eviction and cross-version migration;
- genomic coordinates for `--promoters` FASTA input without a GFF.
