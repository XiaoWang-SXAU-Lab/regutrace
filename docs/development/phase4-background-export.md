# Phase 4 regulatory-background export

**Date:** 2026-11-25

## 1. Purpose

ReguTrace builds a genome-wide regulatory background once and reuses it across rsWGCNA projects. The background is a set of files plus a manifest, not a single archive: the site table is ~2 GB and should not be forced through every project import.

## 2. Manifest

Schema: `regutrace-regulatory-background-v1`.

Top-level fields:

| Field | Contents |
|---|---|
| `schema_version` | `regutrace-regulatory-background-v1` |
| `job_id` | ReguTrace job identifier |
| `software` | software name/version |
| `mapping_mode` | strict / balanced / exploratory |
| `gene_id_namespace` | GFF gene ID namespace for targets |
| `reference_pack` | key, version, species, assembly, checksums |
| `inputs` | role, path/name, size, SHA-256 |
| `promoter_definition` | upstream bp, downstream bp, strand mode |
| `scan` | p-value, background mode, exact-PWM reuse, site-only |
| `mapping_methods` | human-readable method list |
| `tables` | pointers to catalogue, TF–motif, sites, prior |
| `statistics` | summary counts (catalogue excluded) |
| `limitations` | evidence-level and promotion caveats |

The manifest deliberately excludes the full TF catalogue (`statistics.tf_catalogue`) to stay small (3,660 B); the catalogue is a separate table.

## 3. TF catalogue

`kind=tf_catalogue`

Header:

```
tf_id  tf_symbol  family  identification_status  mapping_status  mapping_method
reference_tf_id  motif_count  target_count  testable_motifs  unattainable_motifs
domain_architecture  protein_length  description  diagnostic_reason
```

Properties:

- includes identified TFs with no motif and no target;
- `diagnostic_reason` explains the failure state;
- `target_count` is the mode-filtered candidate count;
- no expression or WGCNA fields.

## 4. TF–motif table

`kind=tf_motif&mode=strict|balanced|exploratory`

Header:

```
tf_id  motif_id  pwm_hash  mapping_status  mapping_method  mapping_confidence
reference_tf_id  donor_species  dbd_identity  query_dbd_coverage  donor_dbd_coverage
specificity_groups  alternative_donor_count  domain_mode  extraction_family
reference_pack_version
```

Properties:

- `pwm_hash` is the SHA-256 of the normalized effective PWM matrix; it identifies exact PWM reuse across motif IDs;
- `mapping_status` and `mapping_confidence` are separate from target evidence;
- `donor_species` is resolved from the motif registry;
- `domain_mode` exposes full-architecture vs single-domain fallback;
- the mode filter uses the same status policy as the UI;
- exploratory mappings are uncalibrated and labeled as such.

## 5. Site table

`kind=sites` (optional, large)

Header:

```
pwm_hash  motif_ids  gene_id  chrom  promoter_start  promoter_end  promoter_strand
hit_start  hit_end  hit_strand  genomic_hit_start  genomic_hit_end  genomic_hit_strand
score  p_value  tf_count  scan_key
```

Properties:

- one row per unique PWM site per PWM group;
- `motif_ids` is a semicolon list for identical matrices;
- no TF identity column; expand with the TF–motif table;
- `scan_key` supports cache/provenance checks;
- this is the optional site-level attachment, not the rsWGCNA prior.

Motif–target associations can be derived by joining `motif_ids` to the TF–motif table and grouping by `gene_id`; the dashboard motif views do this without materializing a second wide table.

## 6. TF–target prior

`kind=prior&policy=...`

| Policy | Included mappings | Testability |
|---|---|---|
| `strict` | `strict_rbh` | testable motifs only |
| `strict_calibrated` | strict + `dbd_calibrated` | testable only |
| `extended` | strict + calibrated + `dbd_exploratory` | testable only |
| `selected` | explicit TF IDs; strict/calibrated by default; exploratory with `include_exploratory=true` | testable only |

Ambiguous mappings are excluded. The prior is produced from evidence rows, so it is only available for non-site-only jobs. Site-only exploratory jobs return an explicit error and the UI explains why.

Output format: canonical SQLite database produced by `regutrace convert`.

### SQLite contract

```sql
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE hits (
  tf_id TEXT, tf_symbol TEXT, target_id TEXT, target_symbol TEXT,
  evidence_type TEXT, evidence_level TEXT, method TEXT, source TEXT,
  pmid TEXT, condition TEXT, promoter_start INTEGER, promoter_end INTEGER,
  motif_id TEXT, hit_start INTEGER, hit_end INTEGER, strand TEXT,
  score REAL, p_value REAL, notes TEXT
);
CREATE TABLE pairs AS
  SELECT tf_id, tf_symbol, target_id, target_symbol,
         COUNT(*) AS n_sites,
         COUNT(DISTINCT motif_id) AS n_motifs,
         MIN(p_value) AS best_p, MAX(score) AS best_score,
         GROUP_CONCAT(DISTINCT motif_id) AS motifs,
         GROUP_CONCAT(DISTINCT evidence_level) AS evidence_levels
  FROM hits GROUP BY tf_id, target_id;
CREATE INDEX hits_tf ON hits(tf_id);
CREATE INDEX hits_target ON hits(target_id);
```

The dashboard Data Import path accepts this canonical `.db` as `tf_target`. The Phase 4 integration test verifies the SQLite header and export path; a full rsWGCNA network run with a Phase 4 prior was not performed.

## 7. Legacy TSV support

- Evidence export remains TSV (`kind=evidence`).
- `regutrace convert --input <tsv|csv|tsv.gz|csv.gz> --out <db|tsv|csv|jsonl> --format ...` supports round-tripping.
- CSV, TSV, gzip, JSONL, and SQLite are supported by the CLI.
- Legacy consumers that expect TSV can use the evidence export; the prior contract is SQLite.

## 8. What the background does not contain

- expression matrices, traits, samples, or WGCNA modules;
- project-specific filtering or ranking;
- experimental validation claims;
- automatic promotion of exploratory mappings;
- a calibrated accuracy estimate for exploratory coverage.
