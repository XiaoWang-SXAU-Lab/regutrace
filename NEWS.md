# News

## Unreleased (2026-11-26)

- **The scan null is now calibrated to the species.** `--background-mode`
  defaults to `empirical`, which fits a strand-symmetric zero-order null to the
  A/C/G/T composition of the scanned promoters. A uniform 0.25-per-base null
  mis-calibrates every p-value on genomes whose promoters are not 25 % each
  base: on HF1 (fitted composition ≈ 67 % AT) observed per-motif promoter
  coverage correlates with the composition-matched expectation at r = 0.94 and
  with the uniform expectation at only r = 0.24, and an A/T-rich Dof motif can
  collect 8 bits of score — over half its information content — from the mismatch
  alone. `--background-mode uniform` still works, warns, and is implied by an
  explicit `--background`. Results produced before this change are uniform-null
  results and are not directly comparable to a rescan. See
  `docs/development/phase6-null-model.md`.
- The scan manifest records `background.calibrated_to_promoter_composition` and a
  `null_model` block, including the chance-coverage floor (~17.8 % of promoters
  at p = 1e-4 on 1 kb promoters), so a result can never be read without knowing
  which null produced it.
- **A background's `n_sites` now counts independent sites.** `regutrace
  background` merges overlapping hit windows before counting, so a degenerate
  match that slides one base at a time is one binding event rather than thirty.
  Previously `n_sites` was `COUNT(*)` over hit windows, which inflated one
  measured motif by 5.0× and five others by 1.01–1.36×. `--site-count
  overlapping` restores the old counting exactly so pre-2026-11 artifacts can be
  rebuilt byte-for-byte; the mode is recorded in `meta.site_count_mode`. This
  changes the meaning of `motif_target.n_sites` and `pairs.n_sites`, so any
  rebuilt background needs a new version label and SHA-256. See
  `docs/development/phase6-target-admission-closure.md` §8.

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
