# Phase 2.6 ablation scripts

Run from the repository root. Inputs are the Phase 2.5 calibration directory and
its BLAST tables.

```bash
CAL=.../phase2_5/calibration
OUT=.../phase2_6

# Species-held-out ablation (primary)
python scripts/phase2_6/ablation.py \
  --dataset "$CAL" --dbd-blast "$CAL/dbd-all-vs-all.tsv" \
  --protein-blast "$CAL/protein-all-vs-all.tsv" \
  --out "$OUT" --fold-scheme species --folds 5

# Alternative held-out units
python scripts/phase2_6/ablation.py ... --out "$OUT/folds-donor" --fold-scheme donor
python scripts/phase2_6/ablation.py ... --out "$OUT/folds-homology" --fold-scheme homology

# Audit
python scripts/phase2_6/audit.py \
  --dataset "$CAL" \
  --ablation-species "$OUT" \
  --ablation-donor "$OUT/folds-donor" \
  --ablation-homology "$OUT/folds-homology" \
  --dbd-blast "$CAL/dbd-all-vs-all.tsv" \
  --out "$OUT/audit"

# Figures
python scripts/phase2_6/make_figures.py --ablation "$OUT" --audit "$OUT/audit" --out "$OUT/figures"

# Synthetic checks
python scripts/phase2_6/test_phase2_6.py
```

The ablation script selects B–E thresholds only on outer training folds. Method
A is fixed at the legacy strict thresholds. Method C uses one global threshold;
D/E use family-specific thresholds. No threshold is copied across methods.

The audit's pairwise and method CIs use homology-cluster bootstrap units. Do not
replace them with pair-level bootstrap intervals.

The exact inputs, counts and outputs are recorded in
`docs/development/phase2.6-report.md`.
