#!/usr/bin/env bash
# Fetch a minimal, versioned Arabidopsis TF–Target data pack.
#
# The script only downloads public reference data; it does not run the tool.
# Usage:
#   bash scripts/fetch-data.sh [destination] [ensembl_release]
#
# Default destination: %USERPROFILE%/_Bio-info-PRJ-local/regutrace-data
# (use a local path, not the Synology Drive folder).

set -euo pipefail

DEST="${1:-$HOME/Bio-info-PRJ/_rsWGCNA-local/regutrace-data}"
RELEASE="${2:-60}"
BASE="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${RELEASE}"
JASPAR="https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt"

mkdir -p "$DEST"
cd "$DEST"

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "exists: $out"
    return
  fi
  echo "download: $url"
  curl -L --fail --retry 3 --connect-timeout 20 -o "$out.part" "$url"
  mv "$out.part" "$out"
}

fetch "$BASE/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz" \
  "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
fetch "$BASE/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${RELEASE}.gff3.gz" \
  "Arabidopsis_thaliana.TAIR10.${RELEASE}.gff3.gz"
fetch "$BASE/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz" \
  "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
fetch "$JASPAR" "JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt"

sha256sum *.gz *.txt > SHA256SUMS 2>/dev/null || true

cat <<EOF

Data pack: $DEST
Genome/GFF: Ensembl Plants release ${RELEASE}
Motifs: JASPAR 2024 CORE plants (CC BY 4.0; cite JASPAR)
SHA256: $DEST/SHA256SUMS

Example (promoter 1 kb, HSF motifs):
  regutrace --genome "$DEST/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz" \\
    --gff "$DEST/Arabidopsis_thaliana.TAIR10.${RELEASE}.gff3.gz" \\
    --motifs "$DEST/JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt" \\
    --tf-map tf_map.tsv --promoter-len 1000 --pvalue 1e-4 --out evidence.tsv
EOF
