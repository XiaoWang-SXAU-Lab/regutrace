# CIS-BP methodology review for Phase 2.6

Date: 2026-11-24  
Source: Weirauch MT, et al. *Determination and Inference of Eukaryotic
Transcription Factor Sequence Specificity.* Cell 158(6):1431–1443, 2014.
DOI: 10.1016/j.cell.2014.08.009. PMID: 25215497. PMC4163041.

This review records what the paper actually describes and where ReguTrace
differs. It does not claim that ReguTrace reproduces CIS-BP.

## What CIS-BP did

### DBD classes and alignment

- They characterized more than 1,000 TFs across 54 DBD classes and 131
  eukaryotes, primarily by PBM.
- Within each DBD class, they aligned DBD sequences with Clustal Omega.
- They calculated pairwise amino-acid percent identity (%ID) over the aligned
  DBD sequences.

### Threshold calibration

- They binned construct pairs by %ID using overlapping bins of size 10 from 0
  to 100, step 5.
- DNA-binding preference similarity was measured as the fraction of shared
  high-scoring 8-mers (E-score > 0.45, with several alternative E/Z-score
  schemes).
- A transfer was "correct" when 8-mer similarity exceeded the 25th percentile
  of experimental replicates for the same protein. This corresponds to a 0.75
  precision target.
- For each DBD class, the inference threshold was the last %ID bin before
  precision dropped below 0.75, using the median across eight 8-mer similarity
  measures.
- Insufficient-data classes received a 70% fallback (the mean/median/mode
  threshold across classes); AT-hook received 40%; IRF, CXXC zinc fingers and
  Dof zinc fingers received 100%.

### Transfer rule

- Motifs, 8-mer scores and consensus sequences were transferred by **simple
  transfer from the nearest characterized protein** above the DBD-class
  threshold.
- The paper describes nearest-neighbor transfer. It does **not** require a
  reciprocal best hit between the query and the reference donor.
- The same motif can therefore be transferred to many TFs whose DBDs are close
  enough to the same characterized donor. CIS-BP explicitly catalogues direct
  and inferred motifs for many TFs.

### Validation

- Leave-one-out cross-validation held out each characterized protein and
  compared it with its nearest characterized neighbor (highest DBD %ID).
- They reported ~89% of inferences above threshold as accurate as experimental
  replicates.
- They compared ortholog-only (88%) and paralog-only (90%) inferences, and
  reported that the DBD identity distribution of their characterized set was
  similar to that of sequenced eukaryotes.
- They used ChIP-seq AUROC and promoter positional analyses as downstream
  support, not as the primary threshold calibration.

## Where ReguTrace differs

| item | CIS-BP 2014 | ReguTrace Phase 2.5/2.6 |
|---|---|---|
| Reference unit | 54 DBD classes | PlantTFHunter/PlantTFDB family labels (32 families in the calibration set) |
| DBD boundaries | Construct DBD alignments defined for the PBM constructs | HMMER domain coordinates plus `plant-dbd-models-v1` architecture rules |
| Identity calculation | Clustal Omega DBD alignment %ID | BLASTP DBD alignment identity and coverage |
| Motif similarity | Fraction of shared high-scoring 8-mers | TOMTOM-style PWM Pearson correlation over offsets/strands, ≥50% overlap |
| Threshold target | 0.75 precision against replicate 25th percentile | 0.90 precision against held-out PWM similarity ≥0.80 |
| Threshold granularity | Per DBD class, with 70% fallback | Per family, with explicit abstention; no universal fallback |
| Reciprocity | Not described; nearest-neighbor transfer | Method B requires reciprocal DBD best hit; methods C/D/E do not |
| Multiple TFs per donor | Allowed by design | Allowed in methods C/D/E; ambiguity preserved |
| Validation | Leave-one-out across characterized proteins; ortholog/paralog breakdown | Donor-, species-, and homology-cluster-held-out; external benchmark not available |
| Taxonomic scope | 131 eukaryotes, 290 genomes | JASPAR CORE plant donors, heavily Arabidopsis |
| Downstream use | 8-mer/PWM transfer for genome-wide inference | Candidate TF–motif mapping; promoter scanning remains a separate evidence level |

## Which CIS-BP assumptions transfer

Transferable with care:

1. DBD similarity is a useful proxy for DNA-binding specificity within a DBD
   class/family.
2. Thresholds should be calibrated per DBD class/family, not globally.
3. Nearest-neighbor transfer does not require reciprocity; multiple TFs may
   share one characterized donor.
4. Threshold choice should be based on an explicit precision target and
   validated against held-out proteins.
5. Ortholog/paralog composition affects apparent performance and must be
   audited.

Not transferable without plant-specific validation:

1. The exact 70% fallback. PlantTFHunter families are not CIS-BP DBD classes,
   and the plant JASPAR donor set is small and Arabidopsis-dominated.
2. The 0.75 precision target and 8-mer E-score >0.45 criterion. ReguTrace uses
   PWM correlation and a 0.90 precision target; these are different labels.
3. CIS-BP's global genome-wide inference scale. ReguTrace calibrates 596 donors
   and 619 experimental motifs; most families cannot support a family-specific
   threshold.
4. The assumption that a nearest neighbor is adequate for every family. Phase
   2.6 shows family- and validation-setting-dependent performance, including
   DBD AUC below 0.5 in ARF, Dof, SBP and MIKC.
5. Any claim that a transferred motif is a validated TF–target edge. It remains
   a motif-mapping candidate and, after promoter scanning, evidence level D.

## Implication for Phase 2.6

The original motivation for non-reciprocal DBD transfer is consistent with
CIS-BP: one experimentally characterized donor can support multiple related
TFs. Phase 2.6 therefore tests non-reciprocal methods C/D/E directly. The
results must be judged by held-out precision and coverage, not by whether they
recover more HF1 TFs.

The paper's DBD-class thresholds are not a license to reuse a universal 70%
threshold in ReguTrace. Phase 2.5/2.6 calibrate each family independently and
abstain where the data do not support a threshold.
