# Input diagnostics

`regutrace validate` reviews genome/GFF/motif/TF-map/annotation/user-list
inputs without scanning. `regutrace scan --report-json` writes the same
diagnostics for a real run. The report is JSON so a frontend or a RAG assistant
can display, filter, and explain it without parsing CLI text.

## Commands

```bash
regutrace validate \
  --genome genome.fa --gff annotation.gff3 \
  --motifs jaspar_plants.txt --tf-map tf_map.tsv \
  --promoter-len 1000 --downstream 0 \
  --report-json report.json

regutrace scan ... --report-json scan-report.json
```

Exit status: `0` when there are no errors (warnings are allowed); `1` when any
diagnostic has severity `error`. `validate --strict` also exits `1` when there
are warnings, emitting `STRICT_WARNINGS`.

## JSON schema

```json
{
  "tool": "regutrace",
  "version": "0.1.0",
  "command": "validate",
  "inputs": {
    "genome": "/data/genome.fa",
    "gff": "/data/annotation.gff3"
  },
  "summary": { "errors": 0, "warnings": 1, "infos": 7 },
  "diagnostics": [
    {
      "severity": "warning",
      "code": "GFF_GENES_WITHOUT_TRANSCRIPT",
      "message": "12 gene(s) have no mRNA/transcript record; TSS uses gene bounds",
      "file": "/data/annotation.gff3",
      "context": { "count": "12" },
      "suggestion": "check whether the GFF uses CDS/exon features without mRNA records"
    }
  ]
}
```

- `severity`: `error`, `warning`, or `info`.
- `code`: stable machine-readable identifier; use it for UI grouping, filtering,
  and assistant prompts.
- `file` / `line`: present when the diagnostic points at an input location.
- `context`: small string key/value map (counts, names, percentages).
- `suggestion`: actionable next step, when one is known.

## Codes

### FASTA

| code | severity | meaning |
|---|---|---|
| `FASTA_SUMMARY` | info | record count, total bases, length range |
| `FASTA_HIGH_N` | warning | more than 5% N bases (assembly gaps) |
| `FASTA_AMBIGUOUS_BASES` | warning | bases outside A/C/G/T/N (IUPAC codes) |
| `FASTA_CHROMOSOMES` | info | sequence names in the FASTA |

### GFF3

| code | severity | meaning |
|---|---|---|
| `GFF_SUMMARY` | info | gene/transcript counts and feature-type histogram |
| `GFF_NO_GENES` | error | no usable gene/transcript models |
| `GFF_INFERRED_MODELS` | info | genes/transcripts inferred from CDS/exon records |
| `GFF_GENES_WITHOUT_TRANSCRIPT` | warning | gene has no mRNA/transcript; TSS uses gene bounds |
| `GFF_TRANSCRIPTS_WITHOUT_FEATURES` | warning | no CDS/exon/UTR; length falls back to genomic span |
| `GFF_MULTI_ISOFORM` | info | multiple isoforms; one representative transcript is selected |
| `GFF_UNUSUAL_STRAND` | warning | strand values other than `+`/`-` |
| `GFF_CHROMOSOMES` | info | chromosome names in the GFF |

### Chromosome matching

| code | severity | meaning |
|---|---|---|
| `CHROMOSOME_NAMES_NO_MATCH` | error | no GFF chromosome matches the FASTA |
| `CHROMOSOME_NAMES_PARTIAL` | warning | some GFF chromosomes are missing from the FASTA |
| `CHROMOSOME_ALIASES_RESOLVED` | info | names matched through `chr`/`chromosome`, `scaffold`, `contig`, case, or zero-padding aliases |

### Motifs and TF map

| code | severity | meaning |
|---|---|---|
| `MOTIF_SUMMARY` | info | motif count and length range |
| `MOTIF_DUPLICATE_IDS` | warning | duplicate motif IDs; only the first matrix is used |
| `MOTIF_UNNAMED` | info | motif header has no name; auto tf-map uses the ID |
| `MOTIF_SHORT` | warning | shortest motif is below 4 bp |
| `MOTIF_LONG` | info | longest motif is above 30 bp |
| `TFMAP_SUMMARY` | info | TF–motif row count and duplicates |
| `TFMAP_MISSING_MOTIFS` | warning | some tf-map motif IDs are not in the motif library |
| `TFMAP_NO_MATCHING_MOTIFS` | error | no tf-map motif ID matches the motif library |
| `TFMAP_DUPLICATE_ROWS` | warning | duplicate `(tf_id, motif_id)` rows |
| `TFMAP_EMPTY_SYMBOLS` | info | empty `tf_symbol`; `tf_id` is used |
| `TFMAP_WITHOUT_MOTIFS` | warning | tf-map given without motifs; matching not checked |
| `TFMAP_AUTO` | info | tf-map derived from motif header names |

### Annotation and user evidence

| code | severity | meaning |
|---|---|---|
| `ANNOTATION_COVERAGE` | info | genes with a symbol in the annotation table |
| `ANNOTATION_NOT_CHECKED` | info | annotation coverage needs gene models |
| `USER_EVIDENCE_SUMMARY` | info | rows by evidence level/type |
| `USER_EVIDENCE_EMPTY` | warning | user evidence list is empty |

### Promoters and run status

| code | severity | meaning |
|---|---|---|
| `PROMOTER_SUMMARY` | info | promoter count and length range |
| `PROMOTER_EMPTY` | error | no promoters could be extracted |
| `PROMOTER_TRUNCATED` | warning | promoter shorter than the requested window |
| `PROMOTER_SKIPPED` | warning | gene produced no promoter (missing chromosome or boundary TSS) |
| `NO_INPUTS` | error | validate received no input files |
| `VALIDATE_OK` | info | input review completed without errors |
| `VALIDATE_FAILED` | error | an input parse or validation step failed |
| `STRICT_WARNINGS` | error | `--strict` was set and warnings were found |
| `SCAN_SUMMARY` | info | genes, promoters, hits, user rows, total evidence rows |
| `SCAN_FAILED` | error | a scan step failed |

## Frontend / assistant use

1. Run `validate --report-json` before a scan (or `scan --report-json` for a
   real run).
2. Show `summary` as a badge and group `diagnostics` by `severity`/`code`.
3. For an assistant, pass the whole JSON plus the user's question. The `code`,
   `context`, and `suggestion` fields are designed to be quoted directly.
4. Treat `error` as blocking; `warning` as a user decision; `info` as context.
