# Phase 2.6 report — mapping strategy ablation and cross-species validation

Date: 2026-11-24  
Repository: `regutrace`  
Base: Phase 0–2.5 working tree (base HEAD `9ac4837`; not committed)

## Scope

Phase 2.6 tests whether the reciprocal DBD best-hit requirement in Phase 2.5 is
scientifically necessary for motif transfer, and whether a non-reciprocal rule
can recover more correct mappings without unacceptable precision loss.

It is an ablation and methodology audit, not a new default mapping mode. The
Phase 2.5 HD-ZIP/TCP calibrated configs are unchanged. No UI or genome-wide scan
was performed.

## Methods compared

All methods use the same 596-donor calibration dataset and the same outer
species-grouped folds. B–E are recalibrated independently inside each outer
training fold; their thresholds are never reused across methods.

| method | evidence used | reciprocity | thresholding | prediction |
|---|---|---|---|---|
| A | full-length protein BLAST | yes, legacy strict | fixed legacy 35% / 0.60 | top donor's motif set |
| B | DBD BLAST | yes | family-specific | top donor's motif set |
| C | DBD BLAST | no | one global threshold | top donor's motif set |
| D | DBD BLAST | no | family-specific | top donor's motif set |
| E | DBD BLAST | no | family-specific | unique specificity group from all eligible donors |

Correctness label: predicted and query motif sets share a PWM with TOMTOM-style
Pearson correlation ≥0.80 and ≥50% overlap. Exact motif IDs are unique to each
donor, so exact-ID recovery is zero by construction and is not used for
threshold selection. The ≥0.90 label is reported as sensitivity.

## Main result — species-held-out ablation

| method | accepted | correct | precision (95% Wilson) | coverage | recall | ambiguous | abstained | unique donors | max queries/donor |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A full-length strict RBH | 139 | 122 | 0.878 (0.813–0.922) | 0.233 | 0.205 | 25 | 432 | 98 | 8 |
| B DBD reciprocal | 90 | 84 | **0.933** (0.862–0.969) | **0.151** | 0.141 | 20 | 486 | 52 | 7 |
| C DBD NN global threshold | 43 | 33 | 0.767 (0.623–0.868) | 0.072 | 0.055 | 2 | 551 | 27 | 4 |
| D DBD NN family thresholds | 162 | 137 | 0.846 (0.782–0.893) | **0.272** | **0.230** | 20 | 414 | 78 | 11 |
| E DBD specificity-group | 44 | 38 | 0.864 (0.733–0.936) | 0.074 | 0.064 | 42 | 510 | 43 | 6 |

Figure:
[strategy-precision-coverage.svg](figures/phase2_6/strategy-precision-coverage.svg).

Homology-cluster bootstrap intervals (500 replicates) are wider than Wilson
intervals:

| method | precision (cluster CI) | coverage (cluster CI) |
|---|---:|---:|
| A | 0.879 (0.811–0.938) | 0.233 (0.179–0.289) |
| B | 0.935 (0.846–1.000) | 0.152 (0.084–0.229) |
| D | 0.844 (0.737–0.934) | 0.274 (0.176–0.374) |
| E | 0.859 (0.667–0.981) | 0.075 (0.033–0.128) |

C is not cluster-bootstrapped because a global threshold was selected in only
one outer fold; the estimate would be too sparse to interpret.

### What this means

1. **Reciprocity buys precision, not coverage.** B vs D: precision 0.933 vs
   0.846, coverage 0.151 vs 0.272, correct calls 84 vs 137. Under a high-
   confidence criterion B is better; under a candidate-recall criterion D
   recovers more correct calls but with more false transfers.
2. **A single global DBD threshold does not work across families.** C selected a
   threshold in only one outer fold; its 0.072 coverage is an abstention result,
   not evidence that global DBD NN is intrinsically worse at the same coverage.
3. **Specificity-group awareness is conservative.** E has 7.0% ambiguity and
   0.074 coverage; the unique-group rule rejects too many queries to be useful
   as a primary mapping strategy.
4. **Full-length strict RBH remains competitive.** A's precision 0.878 is higher
   than D's 0.846, with slightly lower coverage. DBD similarity is not a
   universal improvement over full-length similarity.

## Validation setting is not interchangeable

The same methods evaluated under donor-, species- and homology-cluster-held-out
folds:

| method | donor precision / coverage | species precision / coverage | homology-cluster precision / coverage |
|---|---:|---:|---:|
| A | 0.937 / 0.455 | 0.878 / 0.233 | 0.870 / 0.077 |
| B | 0.945 / 0.369 | 0.933 / 0.151 | 0.825 / 0.096 |
| C | 0.911 / 0.681 | 0.767 / 0.072 | no calls |
| D | 0.922 / 0.622 | 0.846 / 0.272 | 0.727 / 0.258 |
| E | 0.978 / 0.225 | 0.864 / 0.074 | 0.750 / 0.054 |

Figure:
[validation-settings.svg](figures/phase2_6/validation-settings.svg).

Donor-held-out folds retain same-species and homologous donors and are strongly
optimistic. Homology-cluster-held-out is the strictest setting. Species-held-out
is the primary Phase 2.6 setting, but it is not equivalent to an external
dataset. No external target benchmark is available in this phase.

## Family-level behavior

Key families (species-held-out, precision / coverage):

| family | A | B | D | E |
|---|---:|---:|---:|---:|
| HD-ZIP | 1.00 / 0.12 | 1.00 / 0.26 | 1.00 / 0.44 | 0.92 / 0.38 |
| TCP | 1.00 / 0.66 | 1.00 / 0.63 | 0.97 / 0.91 | 1.00 / 0.46 |
| bHLH | 1.00 / 0.17 | 1.00 / 0.44 | 0.89 / 0.78 | no calls |
| WRKY | 1.00 / 0.33 | 1.00 / 0.11 | 0.91 / 0.20 | 1.00 / 0.04 |
| Dof | 1.00 / 0.11 | 0.88 / 0.49 | 0.83 / 0.69 | no calls |
| ERF | 0.83 / 0.21 | 0.83 / 0.06 | 0.82 / 0.16 | no calls |
| bZIP | 0.73 / 0.23 | 0.79 / 0.29 | 0.64 / 0.69 | 0.67 / 0.25 |
| MYB | 0.40 / 0.08 | no calls | no calls | no calls |

"No calls" means the method could not select a family threshold meeting the
precision target in the outer training folds; it is an abstention, not a
demonstration that the family is unmappable.

## CIS-BP comparison

The full review is in
[phase2.6-cisbp-method-review.md](phase2.6-cisbp-method-review.md).

CIS-BP 2014 used Clustal Omega DBD alignments, per-DBD-class %ID thresholds,
8-mer E-score similarity, a 0.75 precision target, and simple nearest-neighbor
transfer. It did not require reciprocal best hits. Its 70% fallback was a
consequence of insufficient per-class data; it was not a universal biological
threshold.

ReguTrace differs: PlantTFHunter family labels, BLASTP DBD alignments, PWM
Pearson correlation, a 0.90 precision target, and explicit family abstention.
Method B adds a reciprocal filter that CIS-BP does not have. The ablation shows
this filter is a useful confidence control, but it is not a reproduction of
CIS-BP.

## Calibration methodology audit

### Leakage and pseudoreplication

- 596 unique donor protein sequences; no repeated protein sequences.
- 484 calibration motifs occur in both JASPAR 2024 and 2026; they are one
  record, not independent evidence.
- 2 exact duplicate motif pairs remain after merging.
- 158 homology clusters at ≥70% DBD identity / ≥70% coverage; the largest
  cluster has 102 donors.
- 6,250 / 11,784 evaluable donor pairs (53%) are within a homology cluster;
  6,569 (56%) are same-species pairs.

Pairwise bootstrap comparison:

| predictor | pair bootstrap CI width | donor bootstrap | homology-cluster bootstrap |
|---|---:|---:|---:|
| DBD identity AUC 0.696 | 0.018 | 0.062 | 0.099 |
| full-length protein AUC 0.608 | 0.021 | 0.064 | 0.132 |

Figure:
[pairwise-bootstrap.svg](figures/phase2_6/pairwise-bootstrap.svg).

The Phase 2.5 pair-level CI was pseudoreplicated: treating donor pairs as
independent understates uncertainty by roughly 5–6× relative to a
homology-cluster bootstrap. Phase 2.6 reports cluster bootstrap intervals for
method metrics as well.

### Family support floor

| min donors | min species | min exact groups | supported families | supported donors |
|---:|---:|---:|---:|---:|
| 10 | 2 | 1 | 13 | 515 |
| 15 | 3 | 2 | 11 | 490 |
| 20 | 3 | 2 | 11 | 490 |
| 30 | 3 | 2 | 9 | 446 |

Across the grid, changing the exact-group minimum from 1 to 3 never changed the
supported set; the species minimum changed it only at the 10-donor floor. The
donor count is the binding constraint. Phase 2.6 retains the Phase 2.5 floor for
comparability and does not lower it to increase coverage.

All 596 calibration donors are `full_architecture`; single-domain fallback
donors were excluded before the ablation. Domain-architecture compatibility is
therefore a fixed eligibility gate in this experiment, not a varied parameter.

### Motif similarity, overlap and information content

At the primary 0.80 similarity label and 0.50 overlap, species-held-out
precision was A 0.878, B 0.933, D 0.846, E 0.864. Tightening the label to 0.90
lowered precision to A 0.713, B 0.856, D 0.747, E 0.773. Changing minimum
overlap from 0.30 to 0.70 lowered precision (for example B 0.944 → 0.900).

Information-content stratification was not monotonic: B quartile 2 had precision
0.739 while quartiles 3–4 were 1.00; A quartile 4 was 0.838. The per-quartile
sample sizes are 9–43 calls, so this is inconclusive rather than evidence for an
IC threshold.

Using all donor motifs versus the maximum-IC representative changed the
correctness label for 0/143 A, 0/90 B, 0/162 D and 1/44 E calls. Representative
selection is not driving the ablation, but this does not justify collapsing
multi-motif donors.

## HF1 implications

Phase 2.5 HF1 results are unchanged:

```text
strict 2024                       260 TFs
strict 2026                       271 TFs
uncalibrated DBD                1,125 TFs / 32,765 rows
calibrated DBD                     19 TFs / 57 rows / 11 ambiguous
strict 2026 + calibrated DBD      280 TFs (+9 HD-ZIP)
```

Phase 2.6 explains why the calibrated DBD set is small: only HD-ZIP and TCP met
the Phase 2.5 nested validation, and the reciprocal filter further restricts
calls. Non-reciprocal DBD transfer would increase candidate coverage, but under
species- and homology-held-out evaluation its precision is below the Phase 2.5
0.90 target. It is therefore not promoted to a release mapping mode.

No new HF1 map was generated from an unvalidated non-reciprocal threshold. The
9 additional HD-ZIP TFs remain the only calibrated contribution over strict
2026.

## Runtime

Ablation wall time on the calibration set (threshold selection included for
B–E):

| method | seconds |
|---|---:|
| A | <0.1 |
| B | 16.9 |
| C | 2.6 |
| D | 17.0 |
| E | 18.2 |

Peak memory was not instrumented for the ablation. The Phase 2.5 calibrated HF1
mapping measured 78.8 MiB peak sampled working set; Phase 2.6 does not add a
production run.

## Tests

Rust suite unchanged from Phase 2.5:

```text
40 CLI unit tests passed
7 CLI integration tests passed
69 core unit tests passed
0 failed
```

Python tests:

```text
Phase 2.5: 7 passed
Phase 2.6: 5 passed
```

Phase 2.6 tests cover disjoint species folds, reciprocal vs non-reciprocal top
donor selection, group-aware ambiguity, homology clustering and motif-overlap
matching.

## Artifacts

Machine-readable ablation and audit outputs:

```text
ablation-methods.tsv
ablation-family.tsv
ablation-oof.tsv
ablation-selected-thresholds.tsv
ablation-sweep.tsv
ablation-similarity.tsv
ablation-runtime.json
folds-donor/
folds-homology/
audit/validation-setting-comparison.tsv
audit/method-cluster-bootstrap.tsv
audit/motif-similarity-sensitivity.tsv
audit/information-content-stratification.tsv
audit/representative-motif-sensitivity.tsv
audit/family-support-grid.tsv
audit/leakage-audit.tsv
audit/pairwise-bootstrap-comparison.tsv
```

The bundle is copied to the sync workspace at:

```text
analysis/regutrace-yeast-validation/phase2.6-ablation-2026-11-24/
```

## Limitations

1. The calibration set is 596 donors / 619 experimental motifs; most families
   have 1–2 species and cannot support a robust species-held-out threshold.
2. The 0.80 PWM-similarity label is operational; 0.90 sensitivity is reported.
   Neither is biological equivalence.
3. No external DAP-seq/ChIP-seq target benchmark was available.
4. Method A uses fixed legacy thresholds; B–E recalibrate. This is intentional:
   changing strict thresholds would no longer be the legacy strict method.
5. The homology-cluster grouping threshold (70% identity / 70% coverage) is a
   leakage-control rule, not a biological threshold.
6. "No calls" in a family means no threshold met the precision target, not that
   the family cannot be mapped with more data.
7. No commit, push or installation package rebuild was performed.

## Next steps

1. Phase 3: do not scan the exploratory 32,765 HF1 rows; implement mapping /
   site-hit / aggregation separation and exact scan reuse.
2. Phase 5: obtain an external plant DAP-seq/ChIP-seq benchmark before promoting
   any non-reciprocal DBD mode or replacing strict 2026 as the default.
3. If a non-reciprocal mode is revisited, require a pre-registered precision
   target under homology-cluster- and species-held-out validation, not HF1
   coverage.
