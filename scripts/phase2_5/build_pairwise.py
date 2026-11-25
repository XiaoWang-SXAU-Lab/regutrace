#!/usr/bin/env python3
"""Build within-family donor-pair DBD, protein and motif similarity table."""
import argparse
import collections
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    motif_similarity,
    read_pfm,
    read_tsv,
    same_genus,
    write_tsv,
)


def best_alignments(path):
    best = {}
    for line in open(path, encoding="utf-8"):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 12:
            continue
        query, subject = fields[0], fields[1]
        if query == subject:
            continue
        bitscore = float(fields[2])
        identity = float(fields[3])
        length = float(fields[4])
        qlen = float(fields[5])
        slen = float(fields[6])
        evalue = float(fields[7])
        qstart, qend, sstart, send = map(float, fields[8:12])
        coverage = min((qend - qstart + 1) / qlen, (send - sstart + 1) / slen)
        entry = {
            "query": query,
            "subject": subject,
            "bitscore": bitscore,
            "identity": identity,
            "coverage": coverage,
            "length": length,
            "evalue": evalue,
        }
        key = (query, subject)
        if key not in best or (entry["bitscore"], entry["identity"]) > (
            best[key]["bitscore"],
            best[key]["identity"],
        ):
            best[key] = entry
    return best


def symmetric_pair(best, a, b):
    left = best.get((a, b))
    right = best.get((b, a))
    if left is None and right is None:
        return None
    if right is None or (left is not None and left["bitscore"] >= right["bitscore"]):
        return left
    return right


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--dbd-blast", type=Path, required=True)
    parser.add_argument("--protein-blast", type=Path, required=True)
    parser.add_argument("--pfm2024", type=Path, required=True)
    parser.add_argument("--pfm2026", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    donors = {
        row["donor_id"]: row
        for row in read_tsv(args.dataset / "calibration-donors.tsv")
    }
    motif_rows = {
        row["motif_id"]: row
        for row in read_tsv(args.dataset / "calibration-motifs.tsv")
    }
    matrices = read_pfm(args.pfm2024)
    matrices.update(read_pfm(args.pfm2026))
    dbd_best = best_alignments(args.dbd_blast)
    protein_best = best_alignments(args.protein_blast)
    by_family = collections.defaultdict(list)
    for donor_id, row in donors.items():
        by_family[row["family"]].append(donor_id)
    cache = {}
    rows = []
    for family, members in sorted(by_family.items()):
        members.sort()
        for a, b in itertools.combinations(members, 2):
            dbd = symmetric_pair(dbd_best, a, b)
            protein = symmetric_pair(protein_best, a, b)
            if dbd is None:
                continue
            motifs_a = [value for value in donors[a]["motif_ids"].split(";") if value]
            motifs_b = [value for value in donors[b]["motif_ids"].split(";") if value]
            best_motif = None
            exact_duplicate = False
            for motif_a in motifs_a:
                for motif_b in motifs_b:
                    if motif_rows.get(motif_a, {}).get("exact_group") == motif_rows.get(
                        motif_b, {}
                    ).get("exact_group"):
                        exact_duplicate = True
                    key = tuple(sorted((motif_a, motif_b)))
                    if key not in cache:
                        if key[0] == key[1]:
                            cache[key] = {
                                "correlation": 1.0,
                                "offset": 0,
                                "overlap": len(matrices[motif_a]),
                                "orientation": "forward",
                            }
                        else:
                            cache[key] = motif_similarity(
                                matrices[motif_a], matrices[motif_b]
                            )
                    similarity = cache[key]
                    if similarity is None:
                        continue
                    if best_motif is None or similarity["correlation"] > best_motif[
                        "correlation"
                    ]:
                        best_motif = {
                            **similarity,
                            "motif_a": motif_a,
                            "motif_b": motif_b,
                        }
            same_species = donors[a]["species"] == donors[b]["species"]
            genus = same_genus(donors[a]["species"], donors[b]["species"])
            rows.append(
                {
                    "family": family,
                    "donor_a": a,
                    "donor_b": b,
                    "species_a": donors[a]["species"],
                    "species_b": donors[b]["species"],
                    "species_relation": (
                        "same_species"
                        if same_species
                        else "same_genus"
                        if genus
                        else "different_genus"
                    ),
                    "dbd_identity": f"{dbd['identity']:.4f}",
                    "dbd_coverage": f"{dbd['coverage']:.4f}",
                    "dbd_bitscore": f"{dbd['bitscore']:.1f}",
                    "dbd_length_a": donors[a]["dbd_length"],
                    "dbd_length_b": donors[b]["dbd_length"],
                    "protein_identity": (
                        f"{protein['identity']:.4f}" if protein else ""
                    ),
                    "protein_coverage": (
                        f"{protein['coverage']:.4f}" if protein else ""
                    ),
                    "protein_bitscore": (
                        f"{protein['bitscore']:.1f}" if protein else ""
                    ),
                    "motif_a": best_motif["motif_a"] if best_motif else "",
                    "motif_b": best_motif["motif_b"] if best_motif else "",
                    "motif_similarity": (
                        f"{best_motif['correlation']:.4f}" if best_motif else ""
                    ),
                    "motif_orientation": best_motif["orientation"] if best_motif else "",
                    "motif_offset": best_motif["offset"] if best_motif else "",
                    "motif_overlap": best_motif["overlap"] if best_motif else "",
                    "motif_overlap_fraction": (
                        f"{best_motif['overlap'] / min(int(motif_rows[best_motif['motif_a']]['length']), int(motif_rows[best_motif['motif_b']]['length'])):.4f}"
                        if best_motif
                        else ""
                    ),
                    "motif_pair_count": len(motifs_a) * len(motifs_b),
                    "exact_motif_duplicate": "yes" if exact_duplicate else "no",
                    "methods_a": donors[a]["methods"],
                    "methods_b": donors[b]["methods"],
                    "motifs_a": ";".join(motifs_a),
                    "motifs_b": ";".join(motifs_b),
                }
            )
    write_tsv(args.out, rows)
    print(f"within-family donor pairs: {len(rows)}")
    print(f"families: {len(by_family)}")


if __name__ == "__main__":
    main()
