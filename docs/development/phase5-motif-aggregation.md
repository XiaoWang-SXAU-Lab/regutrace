# Phase 5 — Candidate motif aggregation and coordinate-aware selection

**Date:** 2026-11-25

## 1. Question

Exploratory HF1 has 1,143 mapped TFs but a median of 17 candidate motifs and 21,715 targets per TF. Are these motifs independent biological DBD hypotheses, or redundant donor/PWM hypotheses for the same DBD region?

## 2. Coordinate audit

A coordinate-aware rerun of `map-dbd-motifs` produced 28,465 exploratory rows with:

- query protein ID;
- query DBD start/end in protein coordinates;
- query alignment start/end;
- donor DBD coordinates;
- DBD identity, bitscore, E-value, query/donor coverage, alignment length;
- mapping status, extraction family/models, domain mode/count;
- motif ID, reference TF, specificity group.

Audit result:

| Threshold | TFs with >1 query DBD group | Max groups/TF |
|---|---:|---:|
| 0.5 | 0 | 1 |
| 0.7 | 0 | 1 |
| 0.8 | 0 | 1 |
| 0.9 | 0 | 1 |

For every coordinate-covered TF, all candidate alignments overlap **one** query DBD region. The overlap threshold is therefore not the active lever on HF1; it is retained as a safeguard for multi-domain TFs in other species.

Composite architectures: RAV is kept as one functional architecture group even if its AP2 and B3 intervals do not overlap; not splitting composite DBDs is deliberate.

## 3. Regression case `Fta_1G0005130`

| Property | Value |
|---|---:|
| Family | bHLH |
| Candidate mapping rows | 19 |
| Query DBD | 152–198 |
| Unique motifs / PWMs | 19 / 19 |
| Specificity groups (0.80) | 7 |
| Query DBD groups | 1 |
| All-candidate targets | 8,709 |
| Best-supported targets | 2,290 |
| Representative-group targets | 4,112 |
| Median pairwise target Jaccard | 0.0018 |
| Max pairwise Jaccard | 1.0 (exact PWM pairs) |

The 19 motifs are donor/PWM alternatives for one bHLH region. Their target sets are mostly disjoint, so taking the union inflates target count almost linearly.

## 4. Policies

Implemented in `scripts/phase5/aggregate_audit.py`:

### all

All original candidate mappings, unchanged. This is the reusable broad background.

### best@0.8

Within each compatible DBD architecture group:

1. remove exact PWM duplicates by PWM hash;
2. rank by architecture compatibility (full before fallback), min(query coverage, donor coverage), DBD identity, then bitscore;
3. choose the top motif; deterministic motif-ID tie-break; ties are recorded in the audit, not silently merged.

Ranking never uses target count or promoter p-value.

### reps@0.8

Within each DBD group, choose the best motif per specificity group (`group_080`). This keeps alternative specificity hypotheses but removes near-duplicate donor motifs.

## 5. HF1 comparison

Denominator 1,907 identified TFs; all policies start from the same 1,143 mapped TFs.

| Metric | All | Best@0.8 | Reps@0.8 |
|---|---:|---:|---:|
| TFs with candidate motifs | 1,143 | 1,143 | 1,143 |
| TFs with ≥1 target | 1,124 | 991 | 1,104 |
| Actual coverage | 58.94% | 51.97% | 57.89% |
| Median motifs/TF | 17 | 1 | 6 |
| Median effective PWMs/TF | 17 | 1 | 6 |
| Median targets/TF | 21,715 | 4,091 | 14,657 |
| Target Q1 / Q3 | 13,891 / 34,348 | 2,290 / 8,385 | 7,644 / 28,711 |
| Total TF–target pairs | 25,481,360 | 7,394,435 | 19,848,826 |
| Untestable selected motifs (sum over TFs) | 2,487 | 162 | 736 |

Best-supported reduces the pair total by 71.0%; representative groups reduce it by 22.1% while keeping 6 motifs. The best median TF retains 25.2% of its all-candidate targets.

Threshold results for 0.5/0.7/0.8/0.9 are identical on HF1 because there is one DBD region per TF. The audit retains all four thresholds for reproducibility and future multi-domain genomes.

## 6. Interpretation

- The large exploratory sets are dominated by the **union of alternative donor-derived motif hypotheses**, not by multiple DBD regions.
- Target overlap between motifs is low, so similar-looking motifs still produce largely different target sets.
- Removing exact PWM duplicates alone is insufficient; specificity-group representatives are the meaningful intermediate.
- Representative selection is an analysis/query/export policy. The original mappings and all-candidate union remain in the background database.
- Ambiguous mappings are visible in the all view but excluded from best/reps target sets by default.
- Do not interpret best-supported as the biologically correct donor; it is a deterministic ranking with preserved tie/ambiguity metadata.

## 7. Reproduce

```bash
python scripts/phase5/aggregate_audit.py \
  --coord-map <phase5>/coord/hf1-exploratory-coord-map.tsv \
  --base-map <phase3>/hf1-exploratory-map.tsv \
  --sites <phase4>/scan-mem-run/sites.tsv \
  --testability <phase4>/scan-mem-run/testability.tsv \
  --pwm-hash <phase2_5>/calibration/motif-specificity-groups.tsv \
  --specificity-groups <phase2_5>/calibration/motif-specificity-groups.tsv \
  --out-dir <phase5>/analysis2
```

Outputs: `phase5-coord-audit.tsv`, `phase5-policy-targets.tsv`, `phase5-policy-summary.json`, `phase5-motif-aggregation.json`.

## 8. Limitations

- Coordinate rerun used 2024 JASPAR and the 2025/2026 donor DBD preparation; 30 base TFs lacked coordinate-row matches (strict/calibrated or donor-set differences) and were handled as opaque single groups in best/reps.
- Overlap is computed on query protein coordinates; donor coordinates are retained but not used for grouping.
- No DNA-contact residue conservation was reliably available; it was not fabricated.
- Specificity groups are PWM-similarity groups, not experimental specificity classes.
- Composite multidomain families beyond RAV may need family-specific rules.
