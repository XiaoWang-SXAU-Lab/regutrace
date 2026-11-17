# Traceability: TF → motif source → target → site

ReguTrace keeps prediction candidates separate from experimental evidence.
The pipeline also writes three small derived tables so a user can follow a
candidate back to its motif provenance without manually joining `tf_map.tsv`,
the reference map, and the JASPAR library.

## Files

| file | one row per | purpose |
|---|---|---|
| `tf_motif_trace.tsv` | local TF × motif | identification status, motif name/library, reference TF donor, mapping method, scan status, ambiguity counts, evidence summary, one example site |
| `motif_ambiguity.tsv` | ambiguous motif | motifs shared by multiple local TFs or mapped through multiple/tied reference donors, with all candidates listed |
| `target_candidate_ties.tsv` | tied site | target + motif + promoter window + hit coordinates + strand supported by more than one TF candidate |

`evidence.tsv` remains the complete site-level table. The trace tables are
derived views; they do not replace or collapse evidence rows.

## Lookup path

1. **TF** — `tf_genes.tsv`: `gene_id`, `tf_symbol`, `family`,
   `identification_status`, `source`.
2. **Motif source** — `tf_motif_trace.tsv`: filter `tf_id` or `tf_symbol`.
   - `motif_id` / `motif_name` / `motif_library` identify the matrix.
   - `reference_donors` is a JSON array; each element has `tf_id`,
     `tf_symbol` and `tf_family`. Use it when one TF–motif row has multiple
     reference donors. The `reference_tf_id` / `reference_tf_symbol` /
     `reference_tf_family` columns are semicolon-joined summaries in donor-ID
     order, not a substitute for the JSON pairing.
   - `mapping_methods`, `mapping_ambiguous`, `domain_modes` show how the
     mapping was inferred.
   - `scan_status` is `scanned`, `excluded_single_domain_fallback`,
     `missing_motif` (matrix not in the supplied motif library), or
     `not_scanned`.
3. **Target** — `evidence.tsv`: filter the same `tf_id` + `motif_id`. Each row
   has `target_id`, `target_symbol`, `promoter_start/end`, and
   `hit_start/hit_end/strand`.
4. **Site ambiguity** — `target_candidate_ties.tsv`: filter `target_id` +
   `motif_id` + `promoter_start` + `promoter_end` + `hit_start` + `hit_end` +
   `strand`. The promoter window is part of the key: the same relative hit in
   two different promoter windows is two distinct genomic sites. `tf_candidates`
   lists all TFs tied at that site.

## Worked example (Drosophila beta)

`tf_motif_trace.tsv` row (abridged):

```text
tf_id=gene-Dmel_CG10021
tf_symbol=gene-Dmel_CG10021
tf_family=zf-C2H2
identification_status=identified
motif_id=MA0454.1
motif_name=odd
motif_library=JASPAR2024_CORE_insects.txt
reference_pack=drosophila-beta-2026-11-15
reference_tf_id=gene-Dmel_CG3851
reference_tf_symbol=odd
reference_donors=[{"tf_id":"gene-Dmel_CG3851","tf_symbol":"odd","tf_family":"zf-C2H2"}]
mapping_methods=reciprocal_dbd_tied_candidate
mapping_ambiguous=true
scan_status=scanned
shared_motif_tf_count=3
shared_motif_tfs=gene-Dmel_CG10021;gene-Dmel_CG3242;gene-Dmel_CG3851
evidence_rows=2409
target_genes=2028
example_target_id=gene-Dmel_CG11887
example_site=694-704(+)
```

This says: the local zf-C2H2 gene `gene-Dmel_CG10021` was assigned the `odd`
matrix `MA0454.1` through the reference donor `gene-Dmel_CG3851`; the
reciprocal DBD hit was tied, so the mapping is ambiguous. Two other local TFs
share the same matrix, and all three are kept as candidates.

`target_candidate_ties.tsv` row (abridged):

```text
target_id=gene-Dmel_CG10001
target_symbol=AstA-R2
motif_id=MA0454.1
motif_name=odd
promoter_start=28744783
promoter_end=28745782
hit_start=896
hit_end=906
strand=-
tf_count=3
tf_candidates=gene-Dmel_CG10021;gene-Dmel_CG3242;gene-Dmel_CG3851
evidence_level=D
```

Filtering `evidence.tsv` on these three TFs + `MA0454.1` + the same promoter
window and hit coordinates returns the three separate prediction rows. The
product rule is to keep all three candidates; do not pick one arbitrarily.

## Ambiguity semantics

- **`mapping_ambiguous=true`** — the same local TF–motif pair has more than one
  mapping row (for example, multiple reference donors), or the reciprocal DBD
  hit was tied. The local TF's matrix assignment is not unique.
- **`shared_motif_tf_count>1`** — multiple local TFs use the same matrix.
  Each produces its own evidence rows; the same genomic site can therefore
  support multiple TF candidates. `shared_motif_tf_count` counts all mappings;
  `shared_motif_scanned_tf_count` counts only mappings actually scanned.
- **`motif_ambiguity.tsv`** groups both situations at the motif level and also
  carries `reference_donors` as a JSON array.
- **`target_candidate_ties.tsv`** is the site-level view. It is the explicit
  table for the rule “when sequence evidence cannot distinguish TFs, preserve
  all candidates and show the tie”.
- `scanned_tf_count` excludes mappings skipped by the fallback filter;
  `tf_count` includes all mappings in the input map. Use `scan_status` in the
  trace table for per-mapping detail. A fallback row is excluded under the
  default full-architecture filter even when the same TF–motif pair also has a
  full-architecture row; `excluded_fallback_rows` counts it separately.

## Structured donor field

When a local TF–motif pair has more than one reference donor, do not pair the
semicolon-joined summary columns by position in a UI. Parse
`reference_donors` instead:

```json
[
  {"tf_id": "D1", "tf_symbol": "BETA", "tf_family": "F2"},
  {"tf_id": "D2", "tf_symbol": "ALPHA", "tf_family": "F1"},
  {"tf_id": "D3", "tf_symbol": "", "tf_family": ""}
]
```

Empty `tf_symbol`/`tf_family` means the reference map had no value for that
donor; the `tf_id` is still authoritative. JSON object field order is not
significant. The same field is present in `motif_ambiguity.tsv`.

These are candidate ambiguities, not evidence of co-binding, cooperativity, or
a validated target. All rows remain evidence level D unless independently
supported.

## Standalone command

The pipeline writes the three tables automatically. For an existing scan or
pipeline output, run:

```bash
regutrace trace \
  --tf-genes run/tf_genes.tsv \
  --tf-map run/tf_map.tsv \
  --evidence run/evidence.tsv \
  --motifs reference-pack/motifs/JASPAR2024_CORE_insects.txt \
  --reference-map reference-pack/reference-map.tsv \
  --pack-version drosophila-beta-2026-11-15 \
  --out-dir run/trace
```

Add `--include-single-domain-fallback` when the scan that produced
`evidence.tsv` used that flag; the default trace behavior mirrors the default
scan (full-architecture mappings only). The standalone command overwrites the
three trace files in `--out-dir`; it does not modify `evidence.tsv`.
