# News

## Unreleased (2026-11-17)

- Added the four-step ReguTrace pipeline: identify local TF genes, match
  reference TFs, map motifs, and scan local promoters. Outputs include the full
  catalogue, TF–motif provenance, site-level evidence, and an SQLite prior.
- Added a multi-species JASPAR 2024 plant reference with 615 donors and one
  representative motif per donor. The plant mapper ranks BLASTP hits within
  the TF family, collapses local isoforms at gene level, and keeps exact
  reciprocal best-hit ties.
- The pipeline accepts CDS or protein FASTA, genome, GFF3, and optional
  normalized gene annotation. Zero identified or mapped TFs produce empty
  evidence outputs with a valid report.
- On HF1, CDS found 1,908 identified TFs and PEP found 1,907. Both gave the
  same 260 TF–motif–donor triples and scanned 257 distinct matrices.
