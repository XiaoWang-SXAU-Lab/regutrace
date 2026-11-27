# Phase 6 design note — target admission (why one MYB "regulates" half the genome)

Status: investigation complete, first step implemented as a read-path filter.
Not a policy change: `phase5-v1` semantics, the schema-v3 export and `pairs` are untouched.

## Observation

HF1 job `Fta_1G0019230` (MYB) was reported with ~36,130 targets. Measured on
`hf1-gui-background-v3-A.db`:

| policy | targets | share of 37,415 promoters |
|---|---:|---:|
| `all` (union over its motifs) | 36,085 | 96.4 % |
| `reps` (16 motifs) | 32,132 | 85.9 % |
| `best` (1 motif: MA1039.1) | 4,099 | 11.0 % |

The TF has 47 motifs that all map to **two** query DBD windows (14–110 and
14–111 — one residue apart), i.e. they do overlap in one region. Coordinate
grouping already collapses them for `best`. De-duplication is therefore not the
missing step: it moves 96 % to 11 %, and 11 % is still implausible.

## Why it happens

1. **Per-site p = 1e-4 against a uniform i.i.d. null admits chance-level
   coverage.** For a 1 kb promoter, both strands and a length-L motif there are
   `2*(1000-L+1)` tests, so the expected share of promoters with at least one
   hit is `1-(1-1e-4)^(2*(1000-L+1)) ≈ 18.0 %`. Observed per-motif rates for this
   TF span 0–39.7 % and sit around 18 %; the `best` motif (MA1039.1) is at
   11.0 %, *below* that expectation.
2. **Each TF inherits the whole family's motif catalogue.** Every one of this
   TF's 47 motifs is shared by 104–131 TFs (mean 117.7). Two sibling MYBs share
   47/48 and 46/52 motifs (Jaccard 0.98 / 0.87). 249 of the 1,379 `best` rows
   (18 %) use a motif shared by ≥100 TFs. Unioning a family catalogue over 1 kb
   promoters saturates the genome; 33,896 of 36,085 targets are hit by ≥2 of the
   47 motifs (max 43).
3. **Short motifs can only be called by imperfect matches at this threshold.**
   Min attainable p for an 8-mer is `4^-8 = 1.53e-5`. All 4,099 targets of
   MA1039.1 have `best_p` in 1.53e-5…9.16e-5 — no exact match anywhere. For the
   same reason a promoter-level BH correction is **structurally unattainable**
   here: rank-1 threshold at q=0.05 over 37,415 promoters is 1.34e-6, below the
   smallest p an 8-mer can produce.

## What already exists (and is not sufficient)

Exact-PWM de-duplication, coordinate-overlap grouping at 0.8 (`best`), and
specificity-group representatives (`reps`) are all *within-TF* motif
de-duplication. None of them addresses acceptance of a single marginal site, and
`reps` barely helps (16 groups → 85.9 %).

## Measured options (best policy, whole background)

| rule | pairs | TFs with ≥1 target | this MYB | median targets/TF |
|---|---:|---:|---:|---:|
| R0 baseline | 7,510,563 | 996 | 4,099 | 5,300 |
| **R1 same motif ≥2 sites** | **2,455,867** | **996** | **360** | **947** |
| R2 score ≥ 0.85·max | 2,357,020 | 980 | 4,099 | 1,375 |
| R3 score ≥ 0.90·max | 1,786,653 | 977 | 2,900 | 841 |
| R4 score ≥ 0.95·max | 1,380,183 | 963 | 1,603 | 556 |
| R5 = R1+R2 | 727,393 | 975 | 360 | 295 |

Notes:

* R1 removes 67 % of pairs and 82 % of the median target count **without losing
  a single TF** — the signature of removing redundancy rather than signal.
* A relative score threshold does not help for the motifs that cause the
  problem: MA1039.1's hits already score 11.11–12.88 against a 12.88 maximum, so
  0.85·max = 10.95 passes everything.
* Tightening p is not a usable substitute: the MYB drops to 0 targets at
  p ≤ 1e-5 while TF coverage falls from 996 to 591.

## Implemented now (read path only)

`TFTargetOptions.MinSites` / API `min_sites` / TF-profile checkbox: a target is
admitted only if **one and the same motif** contributes that many distinct sites
to its promoter. Sites from different motifs deliberately do not add up — that
is what made the union worthless. Default in the TF profile is 2; `0`/`1`
disables it.

Untouched: `tf_motif_policy`, `pairs`, the schema-v3 export, `phase5-v1`.
The exported `pairs` remains the best-policy set (4,099 for this TF), so a
downstream WGCNA prior is unchanged until the decision below is taken.

## Open decision (needs a policy version, i.e. Phase 6)

To make the **export and the WGCNA prior** benefit, one of these must be
versioned as `phase6-v1`:

* adopt R1 as the admission rule for `pairs` (largest effect, smallest loss);
* replace the uniform i.i.d. null with a composition-matched / empirical
  background (`--background-mode empirical` already exists; Phase 5 measured
  16.9 % → 4.7 % false positives at GC 80 %);
* decide whether homology-transferred motif catalogues may be presented as
  TF-specific evidence at all, or only as family-level candidates.

Any of these changes every target count and therefore requires re-running the
DAP-seq benchmark before acceptance. `scripts/phase5/benchmark_dap.py` is the
harness.

## Reproduction

```sql
-- per-motif coverage vs the i.i.d. null (L from motif_registry.matrix_length)
SELECT m.motif_id, (SELECT COUNT(DISTINCT target_id) FROM motif_target x
                    WHERE x.motif_id=m.motif_id) n
FROM (SELECT DISTINCT motif_id FROM tf_motif WHERE tf_id='Fta_1G0019230') m;

-- R1 on the exported DB
SELECT COUNT(*) FROM (SELECT s.tf_id, mt.target_id
  FROM tf_motif_selected s JOIN motif_target mt ON mt.motif_id=s.motif_id
  WHERE mt.n_sites>=2 GROUP BY 1,2);
```

## Resolution (2026-11-26) — measured; R1 is not the fix

Measured on the **frozen** Phase 5 background (`7a33d3fc…cbbee4`, schema v3), not on the GUI job
database used above. Full write-up: `phase6-target-admission-closure.md`; artifacts under
`analysis/regutrace-hf1-phase6-2026-11-26/`.

**Baselines reproduce exactly**, so the measurement is comparable to Phase 5: 25,481,360 /
7,394,435 / 19,848,826 pairs and 1,124 / 991 / 1,104 TFs; medians 21,715 / 4,091 / 14,657 on the
1,143-TF denominator. Note the median denominator: medians over TFs *with* targets are higher
(21,785 / 4,643 / 14,955), which is why the 4,099-vs-4,091 difference above is not a discrepancy.

**R1 at `min_sites=2`** loses no TF in any policy: `all` 49.84 % of pairs and 32.7 % of the median,
`best` 32.92 % / 17.26 %, `reps` 41.58 % / 24.44 %. Threshold 3 already drops TFs (best: −4), so 2
is the only defensible value — and it is what already ships in the UI.

**But R1 does not reach the cases it was written for.** The design-note TF does move
(`Fta_1G0019230`: 96.7 % → 54.6 % of the genome at `all`, 11.0 % → 1.0 % at `best`), yet the worst
offenders keep most of their targets: `Fta_4G0017950`/`Fta_8G0042250` (Dof, MA1274.1) 33,257 →
30,492 (91.7 %), `Fta_8G0002060` 85.7 %, `Fta_1G0030520` 85.4 %. Seven of the ten largest
`best`-policy target sets are Dof and the top ten retain 72–92 %. 57 motifs hit ≥50 % of promoters;
they occupy 84 of the 1,273 best-policy selections and 79 of the 991 TFs.

**The cause is the null model, not admission.** MA1274.1 is a 21 bp, 15.59-bit motif with minimum
attainable p = 2.27e-13, so shortness cannot explain 88.9 % promoter coverage. Recomputing the
scanner's exact score distribution under an AT-rich background explains it: expected coverage
0.174 (uniform) → 0.598 (60 % AT) → 0.832 (65 % AT) → 0.963 (70 % AT) against 0.889 observed. The
8-mer MYB behaves in the opposite direction: 0.166 uniform → 0.100 at 65 % AT against 0.110
observed. Across all 614 motifs with targets, observed coverage correlates with the
composition-corrected expectation at r = 0.937 (uniform: 0.238) and 544/614 fall within ±10 pp of
the 65 % AT expectation.

**Decision.** Do not freeze R1 as `phase6-v1` for `pairs`, the export or the prior: it re-labels
the same saturated predictions under a version number without touching the tail. Keep `min_sites`
as the shipped read-path filter (default 2). `phase5-v1` and schema v3 stay the contract. A
corrected null changes every p-value in every scan and therefore needs a new policy version plus a
fresh DAP-seq run — that is Phase 6 work, not a closure decision.

**Still open.** *(item 1 resolved 2026-11-27 — see below)*
1. ~~`n_sites` is `COUNT(*)` over site rows with no overlap merge (`background.rs::
   aggregate_motif_target`), so adjacent windows of one degenerate match count twice. R1's
   "two independent sites" rationale is therefore unproven until checked against the site table.~~
   **Resolved.** Measured on six motifs: overlap inflation is 1.01–1.36× for five of them and
   **5.00×** for `MA1065.3`. `regutrace background` now merges overlapping windows by default
   (`--site-count independent`; `--site-count overlapping` reproduces the old counting so frozen
   artifacts stay rebuildable), and `n_sites` means independent sites. See
   `phase6-target-admission-closure.md` §8.
2. No DAP-seq re-run: every number here describes the HF1 candidate set, not assay performance.
3. HF1's measured promoter base composition is not stored, so the AT references are hypothetical
   thresholds rather than a fitted value. **Partly resolved:** the empirical-null rescan reports
   `background.values_acgt`, measured at **65.1 % AT**
   (`analysis/regutrace-hf1-phase6-rescan-2026-11-26/`).

**Independent verification.** The composition analysis required a Python mirror of the scanner
null; it reproduces `motif_testability.min_pvalue` for **614/614** motifs exactly. That required
reproducing `pwm.rs::normalize_counts`, which applies a **0.5 pseudocount per base** (denominator
`sum + 2`); without it 4 motifs disagree by one probability quantum. Any future change to that
prior changes the stored `min_pvalue` values.

