# Phase 6 — target admission closure

**Date:** 2026-11-26
**Status:** measurement complete. Recommendation: do **not** version R1 as the `pairs` contract.
**Scope:** analysis of the frozen Phase 5 background only. No algorithm, threshold, schema, or
`phase5-v1` change. No DAP-seq re-run (see §7).

Companion: `regutrace/docs/development/phase6-target-admission.md` (the original design note).

## 1. Question

The design note observed that the HF1 MYB `Fta_1G0019230` is reported with ~36,000 targets
(96 % of the genome) and proposed **R1**: admit a TF→target pair only when one and the same
motif contributes at least *N* distinct sites to that promoter. The open decision was whether
R1 should be frozen as `phase6-v1` and applied to `pairs`, the background export and the WGCNA
prior, rather than remaining a UI read-path filter.

This document answers that with measurements on the **frozen** background
(`hf1-background.db`, SHA-256 `7a33d3fc…cbbee4`, schema v3), not on the newer GUI job database
that the design note used.

## 2. The measurement reproduces the Phase 5 baseline

Required before any R1 number is trusted:

| Policy | Phase 5 frozen | Recomputed here | TFs (frozen / here) |
|---|---:|---:|---|
| `all` | 25,481,360 | 25,481,360 | 1,124 / 1,124 |
| `best` | 7,394,435 | 7,394,435 | 991 / 991 |
| `reps` | 19,848,826 | 19,848,826 | 1,104 / 1,104 |

Median targets per TF on the Phase 5 denominator (1,143 mapped TFs) also reproduce exactly:
21,715 / 4,091 / 14,657. The denominator matters: medians over TFs *with* at least one target
are higher (21,785 / 4,643 / 14,955), so R1 numbers below are always quoted on the 1,143
denominator to stay comparable with `phase5-policy-summary.json`.

## 3. R1 does what the design note said — for the median

`n_sites` here is `motif_target.n_sites`, i.e. the number of distinct site rows of that motif in
that promoter; sites from different motifs deliberately do not add up.

| Policy | min sites | Pairs | % of baseline | TFs with ≥1 target | Median targets | % of baseline median | TFs lost |
|---|---:|---:|---:|---:|---:|---:|---:|
| all | 1 | 25,481,360 | 100.0 | 1,124 | 21,715 | 100.0 | 0 |
| all | 2 | 12,700,976 | 49.84 | 1,124 | 7,100 | 32.7 | 0 |
| all | 3 | 5,884,018 | 23.09 | 1,124 | 2,900 | 13.35 | 0 |
| all | 4 | 3,582,072 | 14.06 | 1,122 | 1,020 | 4.70 | 2 |
| best | 1 | 7,394,435 | 100.0 | 991 | 4,091 | 100.0 | 0 |
| **best** | **2** | **2,434,280** | **32.92** | **991** | **706** | **17.26** | **0** |
| best | 3 | 976,655 | 13.21 | 987 | 100 | 2.44 | 4 |
| best | 4 | 517,926 | 7.00 | 951 | 21 | 0.51 | 40 |
| reps | 1 | 19,848,826 | 100.0 | 1,104 | 14,657 | 100.0 | 0 |
| reps | 2 | 8,253,905 | 41.58 | 1,104 | 3,582 | 24.44 | 0 |
| reps | 3 | 3,309,993 | 16.68 | 1,103 | 1,179 | 8.04 | 1 |
| reps | 4 | 1,730,150 | 8.72 | 1,099 | 435 | 2.97 | 5 |

At `min_sites=2` no TF is lost in any policy, which is what makes R1 attractive. At 3 and above
TFs start disappearing (best: −4 at 3, −40 at 4), so 2 is the only defensible threshold — and it
is also the threshold that already ships in the UI.

960 of the 991 `best`-policy TFs change at R1, so this is not a mild filter; it is a global
rescaling of the candidate set.

## 4. R1 does not fix the cases it was invented for

The design note was written for `Fta_1G0019230`. On the frozen background R1 does transform that
TF (96.7 % → 54.6 % of the genome at `all`; 11.0 % → 1.0 % at `best`). But the ranking of the
worst offenders is unchanged, because the motifs behind them are dense:

| TF (family) | Motif(s) | Len | Baseline targets (`best`) | After R1 | Retained |
|---|---|---:|---:|---:|---:|
| `Fta_4G0017950` (Dof) | MA1274.1 | 21 | 33,257 | 30,492 | 91.7 % |
| `Fta_8G0042250` (Dof) | MA1274.1 | 21 | 33,257 | 30,492 | 91.7 % |
| `Fta_8G0002060` (Dof) | MA1267.2 | 19 | 32,189 | 27,586 | 85.7 % |
| `Fta_1G0030520` (Dof) | MA1281.2 | 19 | 32,127 | 27,422 | 85.4 % |
| `Fta_1G0002020` (Dof) | MA1279.2; MA2708.1 | 16; 11 | 31,105 | 22,450 | 72.2 % |
| `Fta_7G0014660` (CPP) | MA1161.1 | 15 | 30,694 | 25,964 | 84.6 % |

Seven of the ten largest `best`-policy target sets are Dof TFs, and the top ten keep 72–92 % of
their targets. `MA1274.1` (DOF3.6) still "regulates" 81.5 % of the genome after R1.

57 of 614 scanned motifs hit at least half of all promoters. They occupy **84 of the 1,273
`best`-policy motif selections, covering 79 of the 991 TFs** — 8 % of TFs. So the problem is
concentrated, but it is exactly the part of the distribution R1 cannot reach.

## 5. Root cause: the null model, not the admission rule

`MA1274.1` is not a short or information-poor motif: 21 bp, information content 15.59 bits,
minimum attainable p = 2.27e-13. It cannot be "too short to be specific". It hits 88.9 % of
promoters anyway.

The frozen scan scores against a **uniform i.i.d. null** (0.25 per base). Plant promoters are
not 0.25/0.25/0.25/0.25. Recomputing the scanner's exact score distribution under an AT-rich
background gives the expected hit fraction under a composition-corrected null:

| Motif | Len | Observed | Expected @ uniform | Expected @ 60 % AT | @ 65 % AT | @ 70 % AT |
|---|---:|---:|---:|---:|---:|---:|
| MA1274.1 (Dof3.6) | 21 | 0.889 | 0.174 | 0.598 | 0.832 | 0.963 |
| MA1267.2 (Dof) | 19 | 0.860 | 0.176 | 0.533 | 0.737 | 0.886 |
| MA1281.2 (Dof) | 19 | 0.859 | 0.171 | 0.509 | 0.710 | 0.868 |
| MA1279.2 (Dof) | 16 | 0.750 | 0.175 | 0.448 | 0.608 | 0.749 |
| MA1039.1 (MYB, the design-note motif) | 8 | 0.110 | 0.166 | 0.128 | 0.100 | 0.071 |

Two directions, one cause. The Dof motifs hit about five times more promoters than a uniform null
predicts, and at ~65 % AT the expectation rises to 0.61–0.83 against 0.75–0.89 observed — most of
the excess, with 3–14 pp left over. The 8-mer MYB hits *fewer* than a uniform null predicts and is
already at chance under the same corrected null. Neither is a length effect, which is why
tightening `p` or gating on attainable p cannot fix this.

Systematically, for all 614 motifs with targets in the frozen background:

| Background | Median expected hit fraction | Pearson r vs observed | Motifs within ±10 pp | Median observed |
|---|---:|---:|---:|---:|
| uniform 0.25 | 0.172 | 0.238 | 376 / 614 | 0.130 |
| 60 % AT | 0.157 | 0.924 | 536 / 614 | 0.130 |
| **65 % AT** | **0.140** | **0.937** | **544 / 614** | **0.130** |
| 70 % AT | 0.115 | 0.938 | 504 / 614 | 0.130 |

The observed per-motif hit rate is reproduced to within 10 percentage points for 89 % of motifs
by composition alone. The remaining spread is where real signal, dinucleotide structure, repeats,
and overlapping windows live.

**Conclusion:** the uniform i.i.d. null is the defect. It inflates Dof-family motifs and deflates
A/T-poor ones, and R1 is a downstream patch that inherits the inflation. Even under a correct
null, ~2,000 tested windows per promoter at p ≤ 1e-4 gives a ≈ 17.8 % chance-coverage floor per
motif (`1 − (1 − p)^(2·(1000 − L + 1))`, ≈ 17.8 % for any realistic `L`) — so a corrected null is
necessary but not sufficient; the target count still has to be read against that floor.

## 6. Independent verification of the arithmetic

The composition analysis required recomputing the scanner's null in Python. The mirror is
bit-exact: `min_pvalue` recomputed for all **614/614** motifs equals `motif_testability.min_pvalue`
stored by the Rust scanner.

Getting there required reproducing `pwm.rs::normalize_counts`, which adds a **0.5 pseudocount per
base** (denominator `sum + 2`). A plain column normalization disagreed for 4 of 614 motifs
(MA0097.1, MA0584.2, MA0940.1, MA1169.1) by exactly one probability quantum; with the pseudocount
the disagreement is zero. This also documents, for the first time, that the stored `min_pvalue`
values depend on the pseudocount prior and will change if that prior changes.

## 7. Decision

1. **Do not freeze R1 as `phase6-v1` for `pairs`, the export, or the WGCNA prior.** It halves the
   candidate set and cuts the median 83 % without losing a TF, but it leaves the 79 TFs that
   motivated the rule essentially untouched (80–92 % retention). Promoting it would re-label the
   same saturated Dof predictions with a version number and no measurable precision gain.
2. **Keep `min_sites` as the shipped read-path filter with default 2.** It is a legitimate,
   transparent way to suppress single-site calls in the UI. It must not be described as an
   admission standard for the background.
3. **`phase5-v1` and schema v3 remain the frozen contract.** `pairs` stays the best-policy set.
4. The real fix — a composition-aware or Markov-1 null — changes every p-value in every scan and
   therefore needs a new policy version plus a fresh DAP-seq benchmark. That is Phase 6 work, not
   a closure decision.

**Update, same day.** The composition-aware default has now been implemented: `--background-mode`
defaults to `empirical` in the CLI, the pipeline and the dashboard, and the scan manifest records
which null produced a result. See [phase6-null-model.md](phase6-null-model.md). Rust (131) and Go
(74) tests pass. A controlled six-motif probe on the frozen HF1 promoter set is in that document's
§7: a uniform control reproduces the frozen MA1274.1 coverage exactly (0.8889), and under the
empirical null it falls to **0.546**, while five other motifs converge to 0.14–0.24 around the
17.8 % chance level. The full 927-motif rescan and the DAP re-run are still the acceptance gate.

Two things follow for this document. First, the rejection of R1 in §7 was made under the *wrong*
null: MA1274.1 loses 83 % of its sites (362,144 → 60,249) even though its promoter coverage only
halves, so `min_sites` on **non-overlapping** sites deserves a fresh test now that the p-values
mean what they say. Second, the ≈ 17.8 % floor is not a modelling artefact to be tuned away — it
is why no admission rule at `p = 1e-4` can be the whole answer.


## 8. Open items before any admission rule is frozen

- **`n_sites` counts overlapping windows.** `background.rs::aggregate_motif_target` sets
  `n_sites = COUNT(*)` over site rows with no overlap merge, so two adjacent windows of one
  degenerate match count as two sites. R1's stated rationale ("two independent sites") is
  therefore not guaranteed. This must be checked against the site table (local, ~2 GB) before any
  rule treats ≥2 as independent evidence. **Measured 2026-11-26** on six motifs under the corrected
  null: the inflation is 1.01–1.36× for five of them but **5.00×** for `MA1065.3` (5,296 targets
  passing "≥2 sites" against 1,059 with two genuinely independent sites). The bug is real but
  motif-specific. See `analysis/regutrace-hf1-phase6-rescan-2026-11-26/README.md`. Correct query:

  ```sql
  -- per motif: how many promoters have >=2 non-overlapping sites
  SELECT motif_id, gene_id, COUNT(*) AS overlapping_windows
  FROM sites WHERE motif_ids = :motif GROUP BY 1,2 ORDER BY 3 DESC LIMIT 20;
  ```

- **`n_sites` overlap counting — fixed 2026-11-27.** `regutrace background` now merges overlapping
  hit windows before counting. `--site-count independent` (the default) counts mutually
  non-overlapping windows by greedy interval scheduling; `--site-count overlapping` reproduces the
  old `COUNT(*)` exactly, so a frozen artifact can still be rebuilt byte-for-byte. `best_p` and
  `best_score` are still taken over the raw windows — merging changes how many sites there are, not
  which one scored best. The mode is recorded in `meta.site_count_mode`, an optional key so that
  pre-2026-11 databases keep validating. Consequences to carry forward: `motif_target.n_sites` and
  `pairs.n_sites` change meaning, so a rebuilt background needs a new version label and SHA-256 —
  `7a33d3fc…cbbee4` stays a `COUNT(*)` artifact; and `min_sites` in the dashboard now counts
  independent sites, which is what its label always claimed. Four unit tests cover degenerate
  sliding, opposite strands at one position, adjacent-but-disjoint windows, shared-PWM rows and
  coordinate-less rows.

- **No DAP-seq re-run.** Every R1 and null-model number here is a property of the HF1 candidate
  set. Whether R1 (or a corrected null) helps or hurts assays cannot be answered from the frozen
  background; `scripts/phase5/benchmark_dap.py` is the harness and still has to be run.
- **`Fta_1G0019230` is not the worst case.** It was the design note's example, but the Dof TFs in
  §4 are. Future admission work should use a saturated Dof motif as the test case.
- The composition figures are computed at a fixed 65 % AT reference, not at HF1's measured
  promoter composition. The DB does not store per-promoter base frequencies, so the exact AT
  content of the HF1 promoter set is still unmeasured.

## 9. Reproduce

Artifacts and scripts: `analysis/regutrace-hf1-phase6-2026-11-26/` (README, SHA256SUMS).

```bash
python3 -I scripts/phase6_r1_part.py <db> r1.jsonl best 1 2 3 4     # per policy/threshold
python3 -I scripts/phase6_r1_part.py <db> r1.jsonl reps 1 2 3 4
python3 -I scripts/phase6_r1_part.py <db> r1.jsonl all 1 2 3 4
python3 -I scripts/phase6_summary.py <db> r1.jsonl r1-summary.json
python3 -I scripts/phase6_composition_all.py <db> <JASPAR2026_plants> composition-all.json composition-all.tsv
python3 -I scripts/phase6_artifacts.py                              # tables + SHA256SUMS
```

`phase6_composition_all.py` exits with `min_pvalue_mismatches` in its summary; a non-zero value
means the scanner's normalization changed and the mirror must be updated before the composition
numbers are used again.
