#!/usr/bin/env python3
"""Cluster calibration motifs into operational specificity groups.

Primary grouping uses TOMTOM-style Pearson correlation >= 0.90 with at least
50% overlap of the shorter motif. Exact PWM hashes are always merged. The
threshold is an operational grouping rule, not a claim of biological
equivalence; sensitivity groups at 0.80/0.85/0.95 are also written.
"""
import argparse
import collections
import itertools
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import motif_similarity, read_pfm, read_tsv, write_tsv  # noqa: E402

THRESHOLDS = (0.80, 0.85, 0.90, 0.95)
MIN_OVERLAP_FRACTION = 0.50


def cluster_motifs(motifs, similarities, threshold):
    clusters = [[motif] for motif in motifs]
    while True:
        best = None
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                cross = [
                    similarities.get(tuple(sorted((a, b))), -1.0)
                    for a in clusters[i]
                    for b in clusters[j]
                ]
                if not cross:
                    continue
                score = min(cross)
                if score >= threshold and (best is None or score > best[0]):
                    best = (score, i, j)
        if best is None:
            break
        _, i, j = best
        clusters[i].extend(clusters[j])
        clusters.pop(j)
    clusters.sort(key=lambda group: sorted(group))
    return clusters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--pfm2024", type=Path, required=True)
    parser.add_argument("--pfm2026", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    donors = read_tsv(args.dataset / "calibration-donors.tsv")
    motifs = {row["motif_id"]: row for row in read_tsv(args.dataset / "calibration-motifs.tsv")}
    matrices = read_pfm(args.pfm2024)
    matrices.update(read_pfm(args.pfm2026))
    motif_family = collections.defaultdict(collections.Counter)
    for donor in donors:
        for motif in donor["motif_ids"].split(";"):
            if motif:
                motif_family[motif][donor["family"]] += 1
    families = collections.defaultdict(list)
    for motif, counts in motif_family.items():
        families[counts.most_common(1)[0][0]].append(motif)

    similarities = {}
    pair_rows = []
    for family, members in sorted(families.items()):
        members.sort()
        for a, b in itertools.combinations(members, 2):
            if a == b:
                continue
            key = tuple(sorted((a, b)))
            similarity = motif_similarity(matrices[a], matrices[b])
            if similarity is None:
                continue
            overlap_fraction = similarity["overlap"] / min(
                int(motifs[a]["length"]), int(motifs[b]["length"])
            )
            score = (
                similarity["correlation"]
                if overlap_fraction >= MIN_OVERLAP_FRACTION
                else -1.0
            )
            similarities[key] = score
            pair_rows.append(
                {
                    "family": family,
                    "motif_a": a,
                    "motif_b": b,
                    "correlation": f"{similarity['correlation']:.4f}",
                    "orientation": similarity["orientation"],
                    "offset": similarity["offset"],
                    "overlap": similarity["overlap"],
                    "overlap_fraction": f"{overlap_fraction:.4f}",
                    "exact_duplicate": (
                        "yes"
                        if motifs[a]["exact_group"] == motifs[b]["exact_group"]
                        else "no"
                    ),
                }
            )
    write_tsv(args.out / "motif-pairwise-similarity.tsv", pair_rows)

    group_rows = []
    summary_rows = []
    for family, members in sorted(families.items()):
        members.sort()
        by_threshold = {}
        for threshold in THRESHOLDS:
            clusters = cluster_motifs(members, similarities, threshold)
            by_threshold[threshold] = clusters
            sizes = [len(group) for group in clusters]
            summary_rows.append(
                {
                    "family": family,
                    "threshold": f"{threshold:.2f}",
                    "motifs": len(members),
                    "groups": len(clusters),
                    "singleton_groups": sum(1 for size in sizes if size == 1),
                    "largest_group": max(sizes) if sizes else 0,
                    "median_group_size": f"{statistics.median(sizes):.1f}" if sizes else "0",
                }
            )
        for motif in members:
            row = {
                "motif_id": motif,
                "family": family,
                "pwm_hash": motifs[motif]["pwm_hash"],
                "exact_group": motifs[motif]["exact_group"],
                "length": motifs[motif]["length"],
                "information_content": motifs[motif]["information_content"],
            }
            for threshold in THRESHOLDS:
                for index, group in enumerate(by_threshold[threshold], 1):
                    if motif in group:
                        row[f"group_{int(threshold * 100):03d}"] = (
                            f"{family}_T{int(threshold * 100):03d}_{index:03d}"
                        )
                        row[f"group_{int(threshold * 100):03d}_size"] = len(group)
                        break
            group_rows.append(row)
    write_tsv(args.out / "motif-specificity-groups.tsv", group_rows)
    write_tsv(args.out / "specificity-group-summary.tsv", summary_rows)
    print(f"motifs grouped: {len(group_rows)}")
    print(f"families: {len(families)}")
    print(f"motif pairs: {len(pair_rows)}")


if __name__ == "__main__":
    main()
