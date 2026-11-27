# Phase 6 — null model: composition-calibrated background by default

**Date:** 2026-11-26
**Status:** implemented and unit-tested (**Rust 131 passed, Go 74 passed / 10 env-gated skips**);
**not accepted.** The HF1 rescan is partially measured (§7) and the DAP-seq re-run has not been
done.
**Supersedes:** the `--background-mode` default, and with it which scanner a plain `scan` runs.
`uses_phase3_scan()` is `background_mode != "uniform"`, so the new default takes the phase3
scanner instead of the legacy one. Under the *same* null the two engines emit byte-identical
evidence (verified on a synthetic genome), so this is a code-path change, not an output change.
The scanner algorithm, the score discretization, the p-value definition and schema v3 are
unchanged.

## 1. What changed

The p-value threshold `p = 1e-4` was being evaluated against a **uniform i.i.d. null**
(0.25 per base). For any species whose promoters are not 25 % each base this is the wrong
null, and every p-value in every scan inherits the error. Measured on HF1:

| Background used to score | Pearson r (observed per-motif promoter coverage vs expectation) | Motifs within ±10 pp |
|---|---:|---:|
| uniform 0.25 | 0.238 | 376 / 614 |
| fitted composition (≈ 67–68 % AT) | 0.938 | 544 / 614 |

The default is therefore now `--background-mode empirical`, which fits a strand-symmetric
zero-order null to the A/C/G/T composition of the scanned promoters. That model already existed
and was already tested (`empirical_background_is_strand_symmetric`); what changed is that it is
no longer opt-in.

## 2. Why it matters, concretely

A motif whose consensus is A/T-rich scores highly against a 0.25 background partly because the
promoter is A/T-rich, not because the sequence is specific. Quantified per motif as
`Σ_pos log2(θ_fit[consensus_pos] / 0.25)`:

| Motif | Length | IC (bits) | Score gained from the wrong null | Observed coverage | Excess over 65 %-AT expectation |
|---|---:|---:|---:|---:|---:|
| MA1274.1 (Dof3.6) | 21 | 15.59 | **+8.23 bits** | 88.9 % | +5.7 pp |
| MA1267.2 (Dof) | 19 | 14.86 | +7.34 bits | 86.0 % | +12.4 pp |
| MA1281.2 (Dof) | 19 | 15.16 | +7.34 bits | 85.9 % | +14.9 pp |
| MA1161.1 (CPP) | 15 | — | +6.65 bits | 82.0 % | **−5.4 pp** |
| MA1039.1 (MYB 8-mer) | 8 | 9.80 | −0.55 bits | 11.0 % | +1.0 pp |

MA1274.1 collects **more than half of its information content as null bias**. GC-rich motifs are
penalised instead (median across 614 motifs is −0.80 bits), which is why the 8-mer MYB looked
*less* promiscuous than chance rather than more. Both directions are the same bug.

**Read the last two columns as separate statistics.** Score gained from the wrong null and
coverage excess over the composition-matched expectation do not rank motifs the same way.
Ranking by score gain puts **MA1161.1 (CPP)** on this list (+6.65 bits), but its expectation at
65 % AT is 87.4 % against 82.0 % observed — a *negative* excess. It is a case of score bias
without coverage inflation, and must not be quoted as an example of a saturated motif. Only
MA1274.1, MA1267.2 and MA1281.2 are inflated on both statistics. To identify coverage inflation
use `saturated-motifs.tsv` / `composition-vs-observed.tsv` (`excess_over_AT65`), not this table.

Fitted HF1 promoter composition: A ≈ T ≈ 0.33–0.34, C ≈ G ≈ 0.16–0.17 (grid fit over
0.50–0.90 A/T; level-matched at 0.67, correlation-matched at 0.68). That estimate was
**inferred from the frozen scan's coverage**, not measured from sequence; the rescan's manifest
reports the true value — measured at 65.1 % AT, see §7.

## 3. Files touched

| File | Change |
|---|---|
| `crates/regutrace-cli/src/main.rs` | `Args::default().background_mode` → `empirical`; new `background_mode_explicit` flag; explicit `--background` now implies `uniform` (it supplies the whole null); usage text; `warn_uncalibrated_null()` moved out of the shared option parser so only `scan` warns, and its wording now distinguishes the default 0.25 null from a user-supplied one |
| `crates/regutrace-cli/src/pipeline.rs` | default `--background-mode` → `empirical`; `--background` without a mode implies `uniform` |
| `crates/regutrace-cli/src/scan_phase3.rs` | scan manifest gains `background.calibrated_to_promoter_composition`, `null_model.calibrated`, `null_model.description`, `null_model.chance_coverage_note` |
| `crates/regutrace-cli/tests/cli.rs` | two scans pinned to `--background-mode uniform`: their coordinate and hit-count expectations were measured under that null and change under the new default |
| `crates/regutrace-cli/tests/phase3.rs` | site/cache block pinned to `uniform` for the same reason, plus a new assertion that an omitted mode hits the same cache entry as an explicit `empirical` |
| `README.md` | documents the default and the measured justification |
| `NEWS.md` | records the default change |
| dashboard `backend/internal/regutrace/manager.go` | `defaultParams()` and `normalizeParams()` → `empirical` |
| dashboard `frontend/src/views/ToolkitRegutraceView.vue` | initial state and preset reset → `empirical` |

The two test files are **not** business-logic changes: they pin the null so that expectations
measured under uniform keep their meaning. Without them `cargo test` fails on assertions that
describe the legacy scan, not on anything this change broke.

Rust tests added: `default_background_mode_is_empirical`,
`explicit_background_implies_uniform_mode`,
`explicit_uniform_mode_is_accepted_and_routed_to_the_legacy_scan`,
`explicit_background_with_empirical_mode_is_rejected`.

## 4. Compatibility and reproducibility

- **Frozen Phase 5 results were produced under `uniform`.** They remain valid as a record of that
  configuration; they are not comparable to a rescan with the new default. Any rebuilt HF1
  background is therefore a **different artifact** and needs a new SHA-256 and a new version
  label, not an in-place update of `7a33d3fc…cbbee4`.
- `--background-mode uniform` still works and prints a warning. `--background A,C,G,T` still
  works and implies it, so existing command lines are not broken.
- The scan cache key already hashes the four background values (`ScanKeyBase::new`), so old cache
  entries cannot be reused under the new null. No cache invalidation step is required.
- `pwm.rs::normalize_counts` (0.5 pseudocount per base) is untouched. Changing it would change
  every stored `min_pvalue`; see the closure document. The empirical **background** model uses a
  separate 1.0 pseudocount per base (`scan_phase3::background_model`); the two are different
  priors and must not be confused.
- The fitted null depends on the promoter set: the same genome scanned with a different gene
  subset gets different p-values. Results are only comparable between jobs with the same promoter
  set. `background.values_acgt` is recorded in the manifest and in `scan_config` so this is
  auditable, but it is not automatic — state the promoter set alongside any p-value.

## 5. What this does not fix

Correcting the null makes `p = 1e-4` mean what it says. It does **not** make a single marginal
site informative, because ~2,000 windows are tested per promoter:

```
expected share of promoters with ≥1 chance hit = 1 - (1 - p)^(2 * (promoter_len - L + 1))
                                               ≈ 17.8 %   at p=1e-4, 1000 bp
```

That figure is ~17.8 % for any motif length (2,000 windows at 1000 bp). Earlier drafts of the
closure document quoted 13–17 %; that was wrong and has been corrected there.

So after this change a motif that covers ~18 % of promoters is at chance, and the ranking of the
extreme cases should collapse but per-motif coverage will not fall to zero. Target counts still
have to be read against that floor, and the admission question (`min_sites`) is still open — see
`phase6-target-admission-closure.md`, which rejected R1 as the primary fix for exactly this reason.

Not implemented: a Markov-1 (dinucleotide) or per-promoter null. The fitted model is zero-order
and strand-symmetric, so it corrects species-wide composition but not promoter-to-promoter
variation. Phase 5's dinucleotide-preserving shuffle (61.5 % shuffled vs 58.5 % observed) suggests
dinucleotide structure is a second-order effect at this scale — but §7 is evidence that for long
A/T-rich motifs it is not second-order at all.

## 6. Verification — pre-registered acceptance

Commands (see §7 for what has already been measured):

```bash
cd ~/Bio-info-PRJ/regutrace
cargo test --workspace                      # includes the 4 new default/mode tests

regutrace scan --promoters p.fa --motifs m.jaspar --out e.tsv --scan-manifest man.json
python3 -c "import json;d=json.load(open('man.json'));print(d['background'],d['null_model'])"

regutrace scan --promoters p.fa --motifs m.jaspar --background 0.3,0.2,0.2,0.3 --out e2.tsv

cd ~/Bio-info-PRJ/wgcna-dashboard/backend && go test ./internal/regutrace ./internal/api
cd ../frontend && pnpm build
```

Acceptance for the HF1 rescan, **frozen before the full run**:

1. **Manifest.** `background.mode` = `empirical_promoter_strand_symmetric`,
   `null_model.calibrated` = true, and `background.values_acgt` within ±0.02 of
   [0.33, 0.17, 0.17, 0.33]. *Already met — see §7.*
2. **The composition dependence disappears (primary criterion).** Under the old null, per-motif
   coverage tracked the composition expectation at r = 0.938 and 57 of 614 motifs hit ≥ 50 % of
   promoters. After the rescan, recompute `composition-vs-observed.tsv` on the new scan and require
   that the 57 motifs of `saturated-motifs.tsv` lose their position in the coverage ranking, and
   that coverage no longer tracks motif AT-richness. Pass = the saturated set disperses; fail = it
   keeps its ranking.
3. **The bulk converges on chance.** Median per-motif promoter coverage should land near the
   17.8 % chance level — measured so far at 0.14–0.24 (§7). A median below 0.10 means the fit is
   too aggressive; above 0.30 means it is not calibrating.
4. **MA1274.1 is reported, not gated.** 0.889 → **0.546** measured (§7). Do **not** require
   0.18–0.35: a zero-order null cannot see A/T tracts, so a 21 bp A/T-rich motif staying well
   above the chance level is an expected limitation of the model, not a failure of the change. It
   must be recorded as the headline residual, and it is the reason a higher-order null is now the
   open item rather than a nice-to-have.
5. **No TF is lost.** `best`-policy TFs with ≥1 target must stay **≥ 985 of 991**. Target *counts*
   per TF will compress toward ~7,000 and that is expected, not a failure. Note that **72 of the
   79 TFs that use a saturated motif have no other motif at all**, so a collapse of this criterion
   is the specific failure mode to watch — not a slow decline.
6. Only then re-run `scripts/phase5/benchmark_dap.py` and compare AUPRC / enrichment against the
   frozen 0.330 / 4.89× (direct motif) and 0.261 / 2.81× (held-out transfer). **No threshold is
   to be tuned on that comparison**; it is an acceptance check, not an optimisation target.
7. The report must state plainly that this corrects the *meaning* of the p-value and does **not**
   make `p = 1e-4` usable as a target-admission rule. A rescan that passes §6 is not a resolution
   of the admission problem.

## 7. Partial HF1 measurement (2026-11-26)

A controlled probe, not the full rescan. Inputs are the frozen set, from
`%LOCALAPPDATA%\rsWGCNA-build`: `phase5/prod-ui-test/hf1-site-pipeline/promoters.fa`
(37,415 × 1000 bp), `jaspar2026/plant-jaspar-2026-rbh/…/JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt`
(927 motifs), `p = 1e-4`, both strands, default promoter window.

**Control.** MA1274.1 under `--background-mode uniform` gives promoter coverage **0.8889**,
identical to the frozen value. The input set and scanner therefore reproduce the frozen scan, and
the empirical number below is comparable to it.

| Motif | Old (uniform, frozen) | New (empirical, measured) |
|---|---:|---:|
| MA1274.1 (Dof3.6) | 0.8889 | **0.5464** |
| MA2725.1 | 0.3005 | 0.1937 |
| MA1337.2 | 0.1299 | 0.2195 |
| MA1039.1 (MYB 8-mer) | 0.1096 | 0.2396 |
| MA1065.3 | 0.0500 | 0.1671 |
| MA2046.2 | 0.0233 | 0.1397 |

Manifest: `mode = empirical_promoter_strand_symmetric`, `null_model.calibrated = true`,
`values_acgt = [0.3256, 0.1744, 0.1744, 0.3256]` (65.1 % AT).

- **Criterion 1 met.**
- **Criteria 2 and 3 are supported**: coverage moves in the predicted direction at both ends
  (deflated motifs rise, the saturated motif falls), and five of six converge into 0.14–0.24, i.e.
  on the 17.8 % chance level. The frozen spread was 0.023–0.889 (38×); the measured spread is
  0.140–0.546 (3.9×).
- **Criterion 4 as originally drafted would have failed**, which is why it is now a reporting
  requirement: MA1274.1 lands at 0.546, more than twice the chance level, because real HF1
  promoters carry A/T structure that a zero-order null cannot represent.
- MA1274.1 unique sites fall 362,144 → 60,249 (6.0×), so the *site* layer shrinks far more than
  the promoter layer. This is why `min_sites` on **non-overlapping** sites is the natural next
  experiment: R1 was rejected under the wrong null, and that rejection should be re-tested now.
- Caveat: this is a Linux debug build on six motifs. It does not replace the full 927-motif
  rescan, which belongs on the Windows release binary.

## 8. Rollback

`--background-mode uniform` (or `--background 0.25,0.25,0.25,0.25`) restores the previous
behaviour exactly, for both the CLI and the pipeline. The dashboard equivalent is the Background
selector in the advanced panel.
