# DBD motif-transfer calibration methodology (Phase 2.5)

Date: 2026-11-24  
Repository: `regutrace`  
Status: implemented; results in [phase2.5-report.md](phase2.5-report.md)

## 1. Question

Can a query TF's experimentally supported motif be transferred from a
same-family donor using DBD sequence similarity alone, and if so, which
identity/coverage thresholds are supported for each family?

The calibration target is **motif specificity**, not TF–target regulation. A
correct transfer means the predicted motif is in the same operational
specificity neighborhood as the query donor's experimentally supported motif.
It does not mean the motif is bound in vivo, nor that a promoter hit regulates a
target.

## 2. Data sources and deduplication

Inputs:

- JASPAR 2024 CORE plant PFMs and the 2024 donor audit;
- JASPAR 2026 CORE plant PFMs and the 2026 donor audit;
- UniProt donor proteins (2026 pack uses UniProt release 2026_04);
- existing PlantTFHunter-compatible HMM/rules and the `plant-dbd-models-v1`
  family-to-DBD architecture config.

Rules:

1. **Donor identity** is the UniProt accession. Records are merged across
   releases by `(donor_id, motif_id)`; a shared record is counted once.
2. **Motif identity** is the JASPAR matrix ID. The 781 shared 2024/2026 matrix
   IDs have identical PFMs; they are not independent observations.
3. **Exact PWM duplicates** are grouped by a canonical SHA-256 over the
   normalized probability matrix (six decimal places), not by name, family or
   similarity. A donor–motif record is retained once; exact duplicate motif IDs
   are kept as a group for labeling.
4. **Experimental provenance**: calibration ground truth uses methods in
   DAP-seq, ampDAP-seq, PBM, ChIP-seq, ChIP-chip, ChIP-exo, SELEX, SELEX-seq,
   HT-SELEX, CAP-SELEX, bacterial 1-hybrid, and the PBM/CSA/DIP-chip composite.
   `Transfac`, `COMPILED` and empty method fields are excluded from calibration
   ground truth and listed in `calibration-exclusions.tsv`.
5. **Donor eligibility**: one DBD sequence, `full_architecture` mode, supported
   family, at least one experimentally supported motif, no family/species
   conflict. Single-domain fallbacks are not calibrated.
6. **Each retained motif has one donor.** JASPAR matrices are donor-specific;
   no exact motif ID is shared by two donors in this calibration set. Exact
   motif recovery is therefore structurally zero for held-out donors, and PWM
   similarity is the primary correctness label.
7. **Multiple motifs per donor** are preserved. A donor's motif set is the union
   of its experimental motifs across both releases.

The calibration dataset is written to:

```text
calibration-donors.tsv
calibration-motifs.tsv
calibration-dbd.fasta
calibration-family-summary.tsv
calibration-exclusions.tsv
exact-pwm-groups.tsv
```

## 3. DBD extraction

DBD extraction reuses `regutrace dbd-extract` and
`config/dbd-models/plant-v1.json`. Families with variable arrays or composite
DNA binding abstain:

```text
C2H2, C3H, CO-like, CPP, GRAS, LSD, NF-YA, NF-YB, NF-YC, ZF-HD
```

For supported families, the architecture is fixed (single model, repeated model,
or ordered pair). The extractor records coordinates, model list, mode and
domain count. Donor DBD FASTA sequences are checked against the full-protein
slice and against the donor registry checksum.

## 4. PWM similarity

TOMTOM was not available in the offline build. We implement the **correlation
component** of the TOMTOM approach:

- normalize each PWM column to probabilities;
- compare all offsets and both orientations (reverse complement);
- compute Pearson correlation over overlapping columns;
- retain the maximum correlation, offset, overlap and orientation.

The implementation does **not** compute a TOMTOM E-value and does not treat
correlation as biological equivalence. To prevent short-overlap artifacts,
grouping and evaluation require at least 50% overlap of the shorter motif.

Operational labels:

| label | rule |
|---|---|
| exact duplicate | identical canonical normalized matrix hash |
| similar (primary) | max Pearson r ≥ 0.80 and overlap ≥ 0.50 |
| high similarity (sensitivity) | max Pearson r ≥ 0.90 and overlap ≥ 0.50 |

The 0.80 primary threshold is an operational decision, not a biological
constant. It was chosen after inspecting the within-family distribution and is
reported with 0.90 sensitivity in the HF1 comparison.

## 5. Pairwise analysis

Within each family, every donor pair is aligned by BLASTP on:

- the extracted DBD sequences; and
- the full-length donor proteins.

For each pair we record DBD identity/coverage, full-length identity/coverage,
best motif pair similarity, exact-duplicate status, species relationship and
provenance. The predictor comparison uses rank AUC for the label
`motif similarity ≥ 0.80` with overlap ≥ 0.50.

This is a descriptive comparison. Pairwise observations are not independent:
donors from the same species and homologous cluster contribute correlated
pairs. The Phase 2.5 pair-level bootstrap intervals were superseded by Phase 2.6
homology-cluster bootstrap intervals; use the Phase 2.6
`pairwise-bootstrap-comparison.tsv` for uncertainty reporting.

## 6. Calibrated transfer rule

For a query donor and a family-specific threshold `(identity, coverage)`:

1. Consider same-family reference donors with a DBD alignment passing identity,
   coverage and `E ≤ 1e-5`.
2. Rank candidates by BLAST bitscore.
3. For calibrated calls, require the **unique top donor** and a reciprocal DBD
   best hit: the donor's best local hit, within the reference panel plus the
   query, must be the query.
4. Transfer the unique top donor's motif set.
5. Map each predicted motif to its operational specificity group.
6. Status:
   - `dbd_calibrated`: calibrated threshold family, reciprocal top donor, and
     all candidate motif groups collapse to one group;
   - `dbd_ambiguous`: calibrated family but multiple specificity groups,
     missing group annotation, or no unique top donor;
   - `dbd_exploratory`: uncalibrated config;
   - `no_supported_motif`: no candidate/reciprocal/motif.

A calibrated call never upgrades the TF–target evidence level. It remains a
motif-mapping confidence label.

## 7. Threshold selection and leakage control

Family support before calibration:

- at least 15 donors;
- at least 3 species;
- at least 2 exact PWM groups.

Threshold grid:

```text
identity: 30, 35, ..., 95%
coverage: 0.5, 0.6, 0.7, 0.8, 0.9
E-value:  1e-5
```

Primary evaluation uses **nested species-grouped cross-validation**:

- outer folds group donors by species;
- inner leave-one-species-out folds select a threshold on the outer training
  donors only;
- the selected threshold is applied to the held-out outer species;
- threshold selection maximizes coverage subject to precision ≥ 0.90 and a
  minimum number of calibrated calls;
- promotion requires nested out-of-fold precision ≥ 0.90, Wilson lower bound
  ≥ 0.70, and at least `max(5, 10% of family donors)` calibrated calls.

Final deployment thresholds are selected with the same leave-one-species-out
procedure on all calibration donors. The nested out-of-fold metrics are the
independent performance estimate; they are not the same observations used to
fit the final threshold.

A secondary homology-cluster sensitivity analysis groups donors at ≥70% DBD
identity / ≥70% coverage. This 70% value is a **leakage-control grouping rule**,
not a motif-transfer threshold.

## 8. What is not claimed

- No family is calibrated from a universal identity threshold.
- No family is calibrated when its nested validation is weak, even if it has
  many donors.
- The 0.80 motif-similarity label is not a claim that two motifs are
  biologically interchangeable.
- The strict RBH map is not treated as biological ground truth.
- Exact duplicate PWMs and shared JASPAR releases are not counted as
  independent experimental evidence.
- The 29 short motifs whose minimum attainable p-value exceeds the default scan
  cutoff are not fixed by lowering the p-value threshold.

## 9. Reproducible command order

```bash
# 1. DBD extraction for donors (2024 and 2026 assets)
regutrace identify-plant-tfs ... --out donor2026-tf_genes.tsv
regutrace dbd-extract ... --out-tsv donor2026-dbd.tsv --out-faa donor2026-dbd.faa

# 2. Calibration dataset
python scripts/phase2_5/build_calibration_dataset.py ...

# 3. All-vs-all BLAST
blastp -query calibration-dbd.fasta -subject calibration-dbd.fasta ...
blastp -query calibration-proteins.faa -subject calibration-proteins.faa ...

# 4. Pairwise similarity and specificity groups
python scripts/phase2_5/build_pairwise.py ...
python scripts/phase2_5/build_specificity_groups.py ...
python scripts/phase2_5/analyze_pairwise.py ...

# 5. Calibration
python scripts/phase2_5/calibrate.py ...
python scripts/phase2_5/make_figures.py ...

# 6. HF1 reassessment
regutrace map-dbd-motifs ... --threshold-config config/dbd-thresholds/plant-calibrated-v1.json \
  --specificity-groups ... --require-reciprocal-dbd true ...
python scripts/phase2_5/hf1_reassess.py ...

# 7. Synthetic helper tests
python scripts/phase2_5/test_phase2_5.py
```

Exact input paths and output hashes are recorded in
[phase2.5-report.md](phase2.5-report.md).
