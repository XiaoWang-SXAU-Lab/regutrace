# Phase 5 scripts

Reproducible analysis and background export. These do not alter the scanner defaults or legacy strict RBH behavior.

- `aggregate_audit.py` — coordinate-aware grouping, best/representative/all policies, site-level and promoter-level BH, coverage and target counts.
- `build_background.py` — standalone SQLite background (`tf_catalogue`, `tf_motif`, `motif_target`, `motif_registry`, canonical `hits`/`pairs`).
- `add_background_layers.py` — adds normalized `motif_testability`, `promoter`, and `scan_config` tables to an existing background DB (schema v2).
- `benchmark_dap.py` — O'Malley 2016 DAP-seq direct/transferred/empirical benchmark.
- `null_shuffle.py` — dinucleotide-preserving promoter null.
- `test_phase5.py` — grouping/overlap/dedup/policy unit checks.

Typical run (paths abbreviated):

```bash
python scripts/phase5/aggregate_audit.py \
  --coord-map "$LOCALAPPDATA/rsWGCNA-build/phase5/coord/hf1-exploratory-coord-map.tsv" \
  --base-map "$LOCALAPPDATA/rsWGCNA-build/phase3/hf1-exploratory-map.tsv" \
  --sites "$LOCALAPPDATA/rsWGCNA-build/phase4/scan-mem-run/sites.tsv" \
  --testability "$LOCALAPPDATA/rsWGCNA-build/phase4/scan-mem-run/testability.tsv" \
  --pwm-hash "$LOCALAPPDATA/rsWGCNA-build/phase2_5/calibration/motif-specificity-groups.tsv" \
  --specificity-groups "$LOCALAPPDATA/rsWGCNA-build/phase2_5/calibration/motif-specificity-groups.tsv" \
  --out-dir "$LOCALAPPDATA/rsWGCNA-build/phase5/analysis2"
```

`build_background.py --policy best@0.8|reps@0.8|all` selects the importable `pairs` policy. `all` can materialize the large exploratory union and should be used deliberately.

The scripts reuse the Phase 3 site table; no promoter rescanning is required for the policy/FDR comparison.
