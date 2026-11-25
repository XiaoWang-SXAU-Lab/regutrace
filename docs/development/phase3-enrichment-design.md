# Phase 3 motif enrichment interface

Date: 2026-11-25  
Status: implemented as an explicitly exploratory interface

## Motivation

A site table answers "where is this motif?" It does not answer "is this motif
over-represented in a gene set?" Users will ask the second question. Providing a
documented, tested interface is better than letting ad-hoc spreadsheet tests
become the de facto answer.

## Command

```text
regutrace enrich \
  --sites <motif-site-hits.tsv> \
  --promoter-summary <promoter-summary.tsv> \
  --gene-sets <sets.gmt> \
  --out <enrichment.tsv> \
  [--summary-out <summary.tsv>] \
  [--min-set-size 5] [--max-set-size N]
```

## Test definition

- **Universe**: genes in `promoter-summary.tsv`, i.e. promoters that were actually
  scanned. Genes absent from the scan are not silently treated as negative.
- **Hit**: a gene has at least one unique site for the motif in its promoter.
  Site count is collapsed to a binary; two sites are not two independent
  observations.
- **Motif unit**: each motif ID in the site table. Exact-duplicate motif IDs are
  not collapsed by `pwm_hash` in this interface; that would hide accession-level
  ambiguity. A future option can aggregate by PWM hash explicitly.
- **Gene set**: GMT membership intersected with the universe.
- **Test**: one-sided Fisher exact test for over-representation:

```text
a = genes in set with motif hit
b = genes in set without hit
c = genes outside set with motif hit
d = genes outside set without hit
p = P(X >= a | X ~ Hypergeometric(a+b+c+d, a+c, a+b))
```

- **Effect size**: fold enrichment = observed / expected, and an odds ratio with
  a 0.5 Haldane correction for zero cells.
- **Multiple testing**: Benjamini–Hochberg across motifs within one gene set. The
  output records the number of tested motifs in the set.

## Output columns

```text
gene_set, description, set_size_in_universe, universe_size,
motif_id, motif_hit_genes, overlap, expected_overlap,
fold_enrichment, odds_ratio, fisher_p, bh_q, tested_motifs_in_set
```

The optional summary table has one row per gene set:

```text
gene_set, set_size_in_universe, tested_motifs, min_bh_q, motifs_q_le_0.05
```

## Why this is exploratory

The Fisher test is exact for the 2×2 table, but the table is a simplification:

1. **Promoter length and GC**: longer or GC-rich promoters have more sites for
   many motifs. The binary hit reduces but does not remove this bias.
2. **Site dependence**: overlapping windows and multiple sites per promoter are
   collapsed to one bit; the bit is not a calibrated binding probability.
3. **Motif similarity**: correlated motifs (for example MYB family matrices) are
   tested separately. BH does not account for this dependence.
4. **Gene-set overlap**: overlapping gene sets are tested separately. BH is
   within-set, not across sets.
5. **Shared PWMs**: two TFs with the same matrix produce the same site set; their
   motif IDs are not independent evidence.
6. **No expression or chromatin data**: enrichment is sequence-only.

A q-value from this interface is a ranking aid, not a claim that the gene set is
regulatorily enriched in the cell.

## What would make it calibrated

- a length/GC/dinucleotide-matched null distribution for each motif;
- promoter-level or gene-level block resampling;
- motif clustering by PWM similarity with a family-level test;
- gene-set overlap-aware multiple testing;
- independent expression or DAP-seq/ChIP-seq validation.

These are deferred; the interface does not pretend to have them.

## Tests

- `hypergeometric_matches_known_values`: N=10, K=4, n=3 against the exact
  combinatorial values.
- `bh_qvalues_are_monotone`: BH step-up behavior.
- `enrichment_interface_uses_scanned_universe_and_reports_bh`: end-to-end CLI test
  with a 10-gene universe, a 5-gene set, one enriched motif and one non-enriched
  motif.

## HF1 use

No curated HF1 gene-set GMT was supplied with Phase 3. The interface was tested
synthetically and is available for Phase 4/5 once gene sets are defined. Running
it on arbitrary random splits would produce uninterpretable q-values and was not
done.
