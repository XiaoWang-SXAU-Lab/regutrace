# Phase 5 — Multiple testing and promoter-level statistics

**Date:** 2026-11-25

## 1. FIMO audit

FIMO documentation (https://meme-suite.org/meme/doc/fimo.html) and the 2011 paper (Grant, Bailey & Noble, Bioinformatics) establish:

- every motif position is scored with a log-likelihood ratio;
- exact p-values are computed from the motif score distribution (Staden dynamic programming);
- p-values are converted to q-values with Benjamini–Hochberg;
- both strands are scored by default for complementable DNA motifs (`--norc` disables);
- `--thresh` is an output filter; `--qv-thresh` switches the threshold to q-values;
- q-values are approximate if `--max-stored-scores` is reached; the default cap is 100,000 stored matches; the least significant 50% are dropped when the cap is hit.

Important boundary: FIMO q-values are sequence-level occurrence statistics. They are not a biological FDR for regulatory interactions, and they do not model promoter-level dependence or multiple motifs assigned to one TF.

The local scanner computes exact p-values for its discrete score and zero-order background. Its default `p=1e-4` is a prefilter, not a genome-wide FDR.

## 2. Complete tested-window universe

For a motif of length `m`, on 37,415 promoters of length 1,000, both strands:

```
windows_per_promoter = 2 * (1000 - m + 1)
total_tests = 37,415 * windows_per_promoter
```

Examples:

| m | Windows per promoter | Total tests |
|---:|---:|---:|
| 6 | 1,990 | 74,455,850 |
| 10 | 1,982 | 74,156,530 |
| 17 | 1,968 | 73,632,720 |

At p=1e-4, a 10 bp motif has ~7,416 expected chance hits per motif before any correction.

## 3. Site-level BH

Implementation:

- collect exact site p-values per motif from the Phase 3 site table;
- sort hit p-values; rank `i` among hits is the exact rank among all tests for `p≤cutoff`, because every window with `p≤cutoff` is in the site table;
- BH `q_i = min_{j≥i} total_tests * p_j / j`;
- cap at 1;
- mark targets/sites with q≤0.05.

Non-hit windows have p>cutoff and are omitted from the rank list; their terms can only make a hit's q smaller, so the reported hit q is a conservative upper bound. This is **not** ordinary BH on the prefiltered table: the complete tested-window count is in the denominator.

## 4. Promoter-level statistic

For each motif and promoter:

```
p_promoter = min(1, best_site_p * windows_per_promoter)
```

This is a Bonferroni bound over the promoter's tested windows. Then BH across all 37,415 eligible promoters (ranking among promoters with a hit; non-hit promoters have p=1 and can only lower q, so the reported q is conservative).

This ignores dependence between overlapping windows and between motifs. It is a conservative promoter-level screen, not a regulatory FDR.

## 5. HF1 comparison

| Policy | Original p≤1e-4 | Site BH q≤0.05 | Promoter BH q≤0.05 |
|---|---:|---:|---:|
| All candidate motifs | 1,124 TFs / 58.94% | 429 / 22.50% | 417 / 21.87% |
| Best-supported | 991 / 51.97% | 86 / 4.51% | 81 / 4.25% |
| Representative groups | 1,104 / 57.89% | 326 / 17.09% | 287 / 15.05% |

Pair totals:

| Policy | Original | Site FDR | Promoter FDR |
|---|---:|---:|---:|
| All | 25,481,360 | 2,842,315 | 2,174,252 |
| Best | 7,394,435 | 393,946 | 264,936 |
| Reps | 19,848,826 | 1,293,478 | 945,875 |

The collapse is expected: 74M tests at p=1e-4 allow thousands of chance hits per motif. Correction at q=0.05 selects only the strongest long motifs. This is evidence that raw exploratory target counts are candidate background, not validated interactions.

## 6. Dinucleotide-preserving promoter null

Subset: first 200 full-length promoters (200 kb), three benchmark motifs (`MA1079.3`, `MA1664.3`, `MA2027.2`), exact dinucleotide counts preserved via Eulerian shuffle, same scanner, uniform background, p=1e-4.

| Metric | Observed | Shuffled |
|---|---:|---:|
| Promoters with ≥1 hit | 117 / 200 (58.5%) | 123 / 200 (61.5%) |
| Median best p | 3.05e-5 | 4.51e-5 |

The shuffled sequences produce at least as many promoter hits as the real promoters. At this subset/scale, raw promoter hit presence is **not** enriched over composition-matched sequence. This is a warning against interpreting hit counts as biological evidence without a promoter-level null.

## 7. Multiple motif assignments per TF

The site/promoter correction is per motif. A TF with 17 motifs has 17 chances; even if each motif is individually corrected, the TF-level union is not. The best-supported and representative policies address this by reducing the number of hypotheses; the all-candidate layer remains available for downstream methods that model dependence.

We did not multiply across motifs by a Bonferroni factor because the motifs are correlated and the appropriate dependence model is unresolved. Reporting the per-motif corrected layer plus the policy comparison is the defensible choice.

## 8. Reproduce

`scripts/phase5/aggregate_audit.py` computes both corrections from the existing site table. `scripts/phase5/null_shuffle.py` generates the shuffled promoters; the same `regutrace scan` binary with the three-motif map produces the null sites.

## 9. Limitations

- Conservative q-values, not exact FDR under dependence.
- Non-hit p-values are bounded only by the cutoff, not materialized; this is conservative for hit q-values but does not estimate the full null distribution.
- Promoter-level Bonferroni is dependence-naive.
- The null subset is 200 promoters/3 motifs, not a genome-wide null.
- FIMO was audited from documentation/paper; FIMO itself was not installed or executed.
