# Phase 2.5 calibration scripts

All scripts are standard-library Python 3. Use the repository root as the
working directory. Inputs are the two plant reference packs, their JASPAR
metadata, the donor audits and the donor DBD extractions.

## Inputs used for the 2026-11-24 run

```text
plant-jaspar-2024-rbh/
  reference-map.tsv
  reference-proteins.faa
  motifs/JASPAR2024_CORE_plants.txt

plant-jaspar-2026-rbh/
  reference-map.tsv
  reference-proteins.faa
  motifs/JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt

JASPAR 2024 metadata: local jaspar.sqlite from the 2024 release
JASPAR 2026 metadata: jaspar2026-plant-metadata.tsv (fetched by script)
2024 donor audit: donor-audit.tsv
2026 donor audit: donor-audit.tsv
2024 donor DBDs: donor-dbd.tsv / donor-dbd.faa
2026 donor DBDs: donor2026-dbd.tsv / donor2026-dbd.faa
```

The calibration bundle in the sync workspace contains the exact outputs and
hashes; `docs/development/phase2.5-report.md` records the actual commands and
counts.

## Command order

```bash
# 1. Donor DBD extraction (once per pack)
regutrace identify-plant-tfs ... --out donor2026-tf_genes.tsv
regutrace dbd-extract ... --out-tsv donor2026-dbd.tsv --out-faa donor2026-dbd.faa

# 2. Build merged calibration dataset
python scripts/phase2_5/build_calibration_dataset.py \
  --pack2024 <pack2024> --metadata2024db <jaspar.sqlite> \
  --audit2024 <donor-audit-2024.tsv> \
  --dbd2024tsv <donor-dbd.tsv> --dbd2024faa <donor-dbd.faa> \
  --pack2026 <pack2026> --metadata2026 <jaspar2026-metadata.tsv> \
  --audit2026 <donor-audit-2026.tsv> \
  --dbd2026tsv <donor2026-dbd.tsv> --dbd2026faa <donor2026-dbd.faa> \
  --out <calibration>

# 3. All-vs-all BLAST on DBD and full-length proteins
blastp -query <calibration>/calibration-dbd.fasta \
  -subject <calibration>/calibration-dbd.fasta \
  -outfmt "6 qseqid sseqid bitscore pident length qlen slen evalue qstart qend sstart send" \
  -max_target_seqs 1000 -evalue 10 -out <calibration>/dbd-all-vs-all.tsv
blastp -query <calibration>/calibration-proteins.faa \
  -subject <calibration>/calibration-proteins.faa \
  -outfmt "6 qseqid sseqid bitscore pident length qlen slen evalue qstart qend sstart send" \
  -max_target_seqs 1000 -evalue 10 -out <calibration>/protein-all-vs-all.tsv

# 4. Pairwise, specificity groups, predictor comparison
python scripts/phase2_5/build_pairwise.py --dataset <calibration> \
  --dbd-blast <calibration>/dbd-all-vs-all.tsv \
  --protein-blast <calibration>/protein-all-vs-all.tsv \
  --pfm2024 <pack2024>/motifs/JASPAR2024_CORE_plants.txt \
  --pfm2026 <pack2026>/motifs/JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt \
  --out <calibration>/dbd-motif-pairwise.tsv
python scripts/phase2_5/build_specificity_groups.py --dataset <calibration> \
  --pfm2024 <pack2024>/motifs/JASPAR2024_CORE_plants.txt \
  --pfm2026 <pack2026>/motifs/JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt \
  --out <calibration>
python scripts/phase2_5/analyze_pairwise.py \
  --pairwise <calibration>/dbd-motif-pairwise.tsv --out <calibration>

# 5. Calibration and figures
python scripts/phase2_5/calibrate.py --dataset <calibration> \
  --dbd-blast <calibration>/dbd-all-vs-all.tsv \
  --out <calibration> \
  --config config/dbd-models/plant-calibrated-v1.json \
  --threshold-config config/dbd-thresholds/plant-calibrated-v1.json
python scripts/phase2_5/make_figures.py --dataset <calibration> --out <calibration>/figures

# 6. HF1 calibrated mapping
regutrace map-dbd-motifs ... \
  --threshold-config config/dbd-thresholds/plant-calibrated-v1.json \
  --specificity-groups <calibration>/motif-specificity-groups.tsv \
  --require-reciprocal-dbd true \
  --out <hf1>/hf1-dbd-calibrated.tsv \
  --diagnostics-out <hf1>/hf1-dbd-calibrated-diagnostics.tsv
python scripts/phase2_5/hf1_reassess.py ... --out <hf1>

# 7. Synthetic checks
python scripts/phase2_5/test_phase2_5.py
```

## Notes

- The calibration scripts do not download JASPAR or UniProt data. Use
  `scripts/fetch-jaspar-metadata.py` and
  `scripts/build-jaspar-2026-audit.py` for the JASPAR 2026 metadata/donor
  preparation, with local caches.
- Do not add JASPAR matrices or UniProt sequences to the source repository.
- `plant-calibrated-v1.json` intentionally contains only HD-ZIP and TCP. All
  other families must return `NO_CALIBRATED_THRESHOLD`.
