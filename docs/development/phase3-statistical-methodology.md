# Phase 3 statistical methodology

Date: 2026-11-25

## What the scanner p-value means

For a PWM with probability matrix `p[b,k]` and background `q[b]`, the scanner
uses

```text
score(k) = sum_b [ log2(p[b,k] / q[b]) * 1{sequence[k]=b} ]
```

and a discretized score `round(10 * score)`. It builds the exact null distribution
of the discrete score under independent positions with probabilities `q`, then
reports

```text
p = P(null discrete score >= observed discrete score)
```

This is a valid p-value for **one pre-specified window** under the model that
positions are independent and identically distributed with `q`. It is not a
p-value for a promoter, a gene, a motif family or a TF.

The implementation skips windows containing non-ACGT bases. Skipping is not a
biological null; it is an explicit missing-data rule.

## Site, promoter and gene-set layers

### Site level

The site p-value is exact for the discretized score and the supplied background.
Because both strands and overlapping windows are scanned, the number of tests is
larger than the number of independent bases. A site table with 3.37 million rows
is not 3.37 million independent tests.

### Promoter level

`promoter-summary.tsv` reports the number of sites, the best site p-value, and
the best PWM. The minimum p-value is **not** converted to a promoter p-value.
Under the null, the minimum of many dependent, discrete p-values is not uniform;
a Bonferroni correction would be conservative, and an independence correction
would be anti-conservative. A valid promoter-level p-value needs a sequence-based
null that preserves length, GC and dinucleotide composition (for example,
dinucleotide-preserving shuffles or a Markov background) plus a dependence-aware
multiple-testing procedure.

**Decision:** no promoter-level q-values in Phase 3. The promoter table is
descriptive.

### Gene-set level

`regutrace enrich` implements a one-sided Fisher exact test on a binary
motif-hit matrix:

- universe = genes in `promoter-summary.tsv`;
- hit = at least one unique site for the motif in the gene's promoter;
- gene set = GMT membership intersected with the universe;
- BH q-values across motifs within one gene set.

This test is valid for the binary table it is given, but the binary table ignores
promoter length, GC, site count, motif similarity and gene-set overlap. It is
therefore labeled **exploratory**. It is implemented because it is useful for
ranking hypotheses; it is not a calibrated regulatory enrichment claim.

### Multiple testing

No genome-wide site or promoter FDR is reported. Naive BH over millions of
overlapping windows would not control the FDR for the biological question, and
BH over only the significant rows would be invalid. The correct future
approaches are:

1. a sequence-based null that matches promoter length and composition;
2. block resampling of promoters or genes, not sites;
3. hierarchical testing (promoter → gene set) with dependence explicitly modeled;
4. separate nulls for motif families and shared PWMs.

These are deferred, not silently replaced with a nominal q-value.

## Empirical background

`--background-mode empirical` estimates a global strand-symmetric zero-order
model from all supplied promoters:

```text
A = T = (count_A + count_T) / 2
C = G = (count_C + count_G) / 2
p = (symmetrized_count + 1) / (2*(A+C) + 4)
```

Properties:

- **Strand symmetry** avoids a strand-dependent background.
- **Pseudocount 1 per base** keeps every log-odds finite and avoids zero
  probabilities.
- **Ambiguous bases** are not counted; windows containing them are skipped.
- **No parameter-uncertainty correction**: the same promoters estimate the
  background and are scanned. The p-values are conditional on the estimate.
- **Zero-order only**: dinucleotide and GC-window structure are not modeled.

On HF1 the estimate was A/T 0.32555, C/G 0.17445 (65.1% AT), based on
24,361,105 A/T and 13,053,895 C/G bases in 37,415 promoters. Relative to uniform,
the empirical background reduced unique sites by 16.9% (3,369,011 → 2,798,756),
TF–target pairs by 8.5% (1,960,433 → 1,793,597), and unattainable motifs from 29
to 19. The mapping was identical; only the scan changed.

This is a real model choice, not a cosmetic option. A uniform background
over-calls GC-rich motifs in GC-rich promoters; an empirical background reduces
that bias but is still not a dinucleotide null.

## Attainable p-value floor

For a PWM of length `L`, the best possible discrete score is the sum of the
per-position maximum quantized log-odds. Under the null, even that score has a
positive probability. `Scanner::min_pvalue()` returns that probability.

If `min_pvalue > requested_p`, no promoter can produce a hit. Phase 3:

- reports `min_pvalue` and `attainable` for every motif;
- marks these motifs `threshold_unattainable`;
- skips their scan jobs;
- keeps the requested threshold unchanged;
- never reports a q-value for them.

On HF1 JASPAR 2026 uniform, 29 mapped motifs were unattainable at `p=1e-4`; with
the empirical background, 19 were. The complete list and floor are in
`motif-testability.tsv`.

## Independent synthetic validation

The validation script reimplements the scanner in Python from the contract
(probability smoothing, log2 odds, scale 10, exact convolution, upper tail) and
compares it with the CLI on fixed windows of random promoters. It is not a wrapper
around Rust.

### Exact reproduction

| background | motif | IC (bits) | min p | max abs p difference | pass |
|---|---|---:|---:|---:|---|
| uniform | AT6 | 3.83 | 2.44e-4 | 0 | yes |
| uniform | GC6 | 3.83 | 2.44e-4 | 0 | yes |
| uniform | LOW4 | 0.31 | 3.91e-3 | 0 | yes |
| uniform | HIGH10 | 9.53 | 9.54e-7 | 0 | yes |
| AT-rich | AT6 | 3.83 | 1.45e-4 | 1.1e-16 | yes |
| AT-rich | GC6 | 3.83 | 1.14e-5 | 0 | yes |
| AT-rich | LOW4 | 0.31 | 2.03e-3 | 1.1e-16 | yes |
| AT-rich | HIGH10 | 9.53 | 3.99e-7 | 3.3e-16 | yes |

The Python and Rust p-values are identical to floating-point precision, and the
minimum p-values match exactly.

### Null calibration on fixed windows

For 1,000 iid uniform promoters and a fixed window:

| motif | fraction p≤0.05 | binomial tail | fraction p≤0.01 | fraction p≤1e-4 |
|---|---:|---:|---:|---:|
| AT6 | 0.035 | 0.99 | 0.004 | 0 |
| GC6 | 0.036 | 0.99 | 0.006 | 0 |

The fractions are at or below nominal. The test is conservative because the
p-value is discrete (many tied scores). A Kolmogorov–Smirnov test is not
appropriate for a discrete null; the table reports the binomial tail for the
event `p ≤ α`, which is the relevant calibration statement. The KS D values are
retained only as descriptive diagnostics.

### GC-rich promoters and background choice

GC-rich motif GCGCGC, 1,000 iid promoters per GC level, fixed window:

| promoter GC | background | fraction p≤0.05 | binomial tail | fraction p≤0.01 |
|---:|---|---:|---:|---:|
| 0.2 | uniform | 0.002 | 1.00 | 0.000 |
| 0.2 | empirical | 0.027 | 1.00 | 0.003 |
| 0.5 | uniform | 0.035 | 0.99 | 0.003 |
| 0.5 | empirical | 0.035 | 0.99 | 0.003 |
| 0.8 | uniform | **0.169** | **3.7e-43** | 0.046 |
| 0.8 | empirical | **0.047** | 0.69 | 0.004 |

At 80% GC, the uniform model produces 16.9% false positives at the 5% threshold;
the empirical model produces 4.7%. This is the strongest justification for
offering `--background-mode empirical`. It does not make the empirical model a
complete null; it only removes the first-order base-composition bias.

### Overlap dependence

A 5,000 bp random sequence, 6 bp motif, 4,995 overlapping windows:

| motif | lag-1 score correlation | lag-5 | lag-10 |
|---|---:|---:|---:|
| ACGTAC | -0.278 | -0.054 | -0.003 |
| AAAAAA | 0.840 | 0.198 | 0.018 |

The sign depends on motif periodicity, but the magnitude shows that adjacent
windows are not independent. This is why site-level BH over all windows is not a
valid genome-wide FDR. The negative correlation for ACGTAC is not an error: a
one-base shift changes which position-specific base is preferred.

### Discretization

For ACGTAC (6 bp), ACGT (4 bp) and ACGTACGT (8 bp), exhaustive enumeration of all
4^L sequences gave a maximum absolute difference of 0 between the discretized
p-value and the continuous-score tail. This does not prove zero error for every
matrix, but it shows that the scale-10 discretization is not a material source of
error for the tested motifs. The scanner's p-value remains exact for the
discretized score by construction.

### Attainable floor

| motif length | IC (bits) | minimum p | attainable at 1e-4 |
|---:|---:|---:|---|
| 4 | 3.14 | 3.91e-3 | no |
| 6 | 4.72 | 2.44e-4 | no |
| 8 | 6.29 | 1.53e-5 | yes |
| 10 | 7.86 | 9.54e-7 | yes |
| 12 | 9.43 | 5.96e-8 | yes |
| 16 | 12.58 | 2.33e-10 | yes |
| 20 | 15.72 | 9.09e-13 | yes |

Length and information content jointly determine the floor; length alone is not
sufficient.

## FIMO comparison

FIMO is not installed in the local environment. The HF1 workspace contains a
micromamba package cache, not a runnable `fimo` binary. No FIMO comparison was
run; no claim of agreement is made. A future comparison must match the
background, pseudocount, strand mode and threshold; otherwise differences are not
interpretable.

## Limitations

1. Zero-order background; no dinucleotide, GC-window or chromatin model.
2. Empirical background is estimated from the scanned promoters; no uncertainty
   propagation.
3. No site/promoter q-values; only exact single-window p-values and descriptive
   aggregation.
4. Discrete p-values are conservative; continuous approximations would be
   anti-conservative at the tail.
5. Overlapping windows and shared PWMs induce dependence; no independence
   assumption is made in the claims.
6. Enrichment is Fisher/BH on a binary hit matrix, not a calibrated regulatory
   model.
7. No external DAP-seq/ChIP-seq calibration in Phase 3.
