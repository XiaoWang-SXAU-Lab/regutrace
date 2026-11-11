# Evidence model

ReguTrace stores one row per evidence record. A TF–target pair may have several
rows (prediction + DAP-seq + ChIP-seq + EMSA), each with its own provenance.

## Levels

| level | meaning | examples |
|---|---|---|
| A | point-to-point validation | EMSA, reporter assay, ChIP-qPCR, mutant/genetic validation |
| B | broad in vivo binding | ChIP-seq, ChIP-chip |
| C | broad in vitro binding | DAP-seq, PBM, SELEX/HT-SELEX |
| D | computational prediction | promoter PWM, motif, homology |
| E | functional/associative evidence | perturbation, expression correlation, module membership |

DAP-seq is **in vitro**; ChIP-seq is **in vivo**. A prediction alone is D and
must not be reported as validated binding.

## TSV / SQLite fields

| field | meaning |
|---|---|
| `tf_id` | TF gene ID |
| `tf_symbol` | TF symbol/name |
| `target_id` | target gene ID |
| `target_symbol` | target symbol/name |
| `evidence_type` | controlled method type (e.g. `prediction_motif`, `in_vitro_dapseq`, `in_vivo_chipseq`, `point_validation_emsa`, `user_curated`) |
| `evidence_level` | A/B/C/D/E |
| `method` | method detail (e.g. `jaspar_pwm`, `DAP-seq`, `EMSA`) |
| `source` | database/paper/import source |
| `pmid` | PubMed ID, if known |
| `condition` | condition/time/genotype label |
| `promoter_start` / `promoter_end` | 1-based inclusive genomic promoter window |
| `motif_id` | JASPAR/PWM ID |
| `hit_start` / `hit_end` | 1-based inclusive motif hit relative to the oriented promoter |
| `strand` | `+`/`-` relative to the oriented promoter |
| `score` | raw PWM log2-odds score (prediction) or source score |
| `p_value` | discretized PWM null p-value (prediction; approximate) |
| `notes` | free text |

SQLite adds `id`, `record_key` (unique import key) and `imported_at`.
Duplicate imports are ignored via `record_key`.

## Coordinate conventions

- GFF3 coordinates are parsed as 1-based inclusive.
- The promoter is oriented 5'→3' in the gene's direction; for minus-strand genes
  the genomic window is reverse-complemented.
- `hit_start/hit_end` are positions in the oriented promoter sequence.
- `promoter_start/promoter_end` are the genomic window (1-based inclusive).
- With `--downstream > 0` the downstream count includes the TSS base (the window
  is symmetric on both strands); with the default 0 it is upstream-only and
  excludes the TSS base.

## Adding evidence

Use the unified TSV with `--evidence`/`--user-list`, or a user list with columns
`tf,target,evidence_type,evidence_level,method,source,pmid,condition,score,notes`.
Import into SQLite with `regutrace import`. `evidence_level` is required for user
lists and must be one of A/B/C/D/E; missing or invalid values are rejected (lower
case is normalized to upper case). Keep `evidence_level` honest; do not label
predictions as A/B/C.
