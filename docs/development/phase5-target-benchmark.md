# Phase 5 — Independent target-level benchmark

**Date:** 2026-11-25

## 1. Data

| Item | Value |
|---|---|
| Assay | O'Malley et al. 2016, *Cell*, Arabidopsis DAP-seq (GEO GSE60141) |
| Peaks | per-TF GEM `narrowPeak.gz` files from the GEO supplementary directory |
| TFs used | MYB44, WRKY21, ERF15, GBF3, TGA1, TCP24, RAV1, HSFA6B (PDF2 explored but no direct JASPAR motif) |
| Genome | *Arabidopsis thaliana* TAIR10 (Ensembl plants release 60 GFF3) |
| Promoter | 1,000 bp upstream, 0 downstream, both strands |
| Positive gene | promoter interval overlaps at least one assay peak |
| Unknown | genes without a peak — **not** established negatives |
| Testing universe | 33,309 promoters with a TAIR10 gene model |
| Leakage control | setting 2 removes all benchmark TF symbols from the donor reference map before DBD transfer |

Peak files are 201 bp GEM events; the benchmark uses the peak interval, not the summit alone. This is promoter-overlap evaluation; distal peaks are not associated to genes.

## 2. Settings

### Setting 1 — direct experimental motif

Use the TF's JASPAR 2024 motif from the plant reference map. This evaluates promoter scanning and target ranking when the motif is known; it is not a test of motif inference.

### Setting 2 — end-to-end transferred motif

Remove all benchmark TFs from the reference map, extract each query DBD with `plant-dbd-models-v1`, run `map-dbd-motifs` against the remaining donors, choose the top-ranked full-architecture motif by coverage/identity/bitscore, and scan. This is a held-out motif-transfer test; the query TF's own motif is not available to the mapper.

## 3. Direct-motif results (uniform background)

| TF | Peaks | Positive genes | Predicted targets | TP | PPV_assay | Recall | Enrichment | AUPRC | Top-100 recovery |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TCP24 | 3,352 | 1,054 | 1,366 | 569 | 0.417 | 0.540 | 13.2× | 0.369 | 0.93 |
| WRKY21 | 1,537 | 723 | 12,281 | 621 | 0.051 | 0.859 | 2.33× | 0.330 | 0.81 |
| ERF15 | 9,908 | 2,695 | 1,761 | 1,245 | 0.707 | 0.462 | 8.74× | 0.434 | 0.74 |
| GBF3 | 6,267 | 3,302 | 3,973 | 1,926 | 0.485 | 0.583 | 4.89× | 0.466 | 0.71 |
| HSFA6B | 25,949 | 6,824 | 8,194 | 2,454 | 0.299 | 0.360 | 1.46× | 0.272 | 0.34 |
| TGA1 | 4,089 | 1,807 | 2,362 | 883 | 0.374 | 0.489 | 6.89× | 0.255 | 0.22 |
| MYB44 | 15,418 | 5,910 | 4,902 | 1,843 | 0.376 | 0.312 | 2.12× | 0.297 | 0.68 |

`PPV_assay` is TP / predicted; it is not true precision because unpeaked genes are unknown, not negatives. `Enrichment` is PPV_assay / assay prevalence. AUPRC uses all 33,309 genes with unpeaked genes scored as non-positive; it is also optimistic/ambiguous.

## 4. Transferred-motif results

| TF | Motif donor | Predicted targets | TP | PPV_assay | Recall | Enrichment | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|
| RAV1 | Q9C6M5 / MA1800.1 | 7,197 | 1,285 | 0.179 | 0.326 | 1.51× | 0.151 |
| TCP24 | Q9LEZ9 / MA1290.1 | 1,692 | 413 | 0.244 | 0.392 | 7.71× | 0.250 |
| WRKY21 | O04609 / MA1303.2 | 10,129 | 617 | 0.061 | 0.853 | 2.81× | 0.335 |
| ERF15 | Q9LND1 / MA1049.2 | 1,071 | 745 | 0.696 | 0.276 | 8.60× | 0.298 |
| HSFA6B | O81821 / MA2019.2 | 7,518 | 2,285 | 0.304 | 0.335 | 1.48× | 0.261 |
| TGA1 | Q39162 / MA1335.1 | 4,807 | 1,226 | 0.255 | 0.678 | 4.70× | 0.449 |
| MYB44 | Q9SN12 / MA0575.2 | 2,400 | 869 | 0.362 | 0.147 | 2.04× | 0.231 |

GBF3 had no supported DBD transfer (`BELOW_DBD_THRESHOLD`); RAV1's direct motif MA0582.2 was unattainable at p=1e-4, so only the transferred setting evaluates RAV1. PDF2 lacked a direct JASPAR motif in the benchmark map and is a data limitation, not a negative result.

Common-TF comparison (6 TFs):

| Metric | Direct | Transfer |
|---|---:|---:|
| Median AUPRC | 0.314 | 0.279 |
| Median enrichment | 4.61× | 3.75× |
| Median recall | 0.475 | 0.363 |
| Median PPV_assay | 0.375 | 0.279 |

Direct motif lookup is better on average, but transfer retains useful signal and is better for TGA1. The difference is a realistic estimate of the cost of not knowing the query motif.

## 5. Uniform vs empirical background

Direct-motif scans were repeated with `--background-mode empirical` (zero-order promoter A/C/G/T):

| Metric | Uniform | Empirical |
|---|---:|---:|
| Median AUPRC | 0.330 | 0.346 |
| Median enrichment | 4.89× | 3.20× |
| Median recall | 0.489 | 0.685 |
| Median PPV_assay | 0.376 | 0.236 |

Empirical background increases sensitivity/recall and predicted target count, lowers PPV_assay and enrichment, and gives a small AUPRC increase. It is not a universal improvement; the default remains uniform.

## 6. Interpretation and limitations

- DAP-seq binding is not proof of functional transcriptional regulation.
- Unpeaked genes are unknown; no true negatives exist. AUPRC/PPV values are comparative, not biological precision.
- Promoter-overlap association ignores distal enhancers and peak-gene looping.
- The benchmark is 7 TFs from one assay series; it cannot validate all 1,143 exploratory TFs or all families.
- A separate Arabidopsis strict-RBH vs calibrated-DBD arm was not run: the plant pack is donor-centric and the benchmark TFs are not all reference donors. The direct-motif arm is the motif-given comparator; the transfer arm is the held-out exploratory comparator. This is a data limitation, not a negative result.
- No threshold was optimized on these results. The scanner's p=1e-4 and default background are unchanged.
- The direct setting uses the TF's own JASPAR motif; the transfer setting excludes all benchmark TFs from donors but still uses the same reference pack and DBD model.
- The benchmark confirms that motif quality and target ranking contain real signal, while also showing that raw target counts are not a precision estimate.

## 7. Reproduce

Inputs and outputs live under `%LOCALAPPDATA%\rsWGCNA-build\phase5\dap\`:

- `peaks/` — nine downloaded GSE60141 narrowPeak files;
- `run1/tf_map.tsv`, `tf_genes.tsv` — direct motif map;
- `run1/evidence.tsv`, `sites.tsv`, `promoter-summary.tsv` — uniform scan;
- `run1/evidence-empirical.tsv` — empirical scan;
- `transfer/` — held-out DBD mapping, best map, transfer evidence;
- `run1/dap-benchmark*.tsv/json` — results.

`scripts/phase5/benchmark_dap.py` computes all metrics from the promoter summary and evidence TSV.
