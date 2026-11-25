#!/usr/bin/env python3
"""Compare DBD and full-length protein similarity as motif-specificity predictors."""
import argparse
import collections
import json
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_tsv, write_tsv  # noqa: E402


def auc(rows, score_field, label_field="label"):
    positives = sorted(row[score_field] for row in rows if row[label_field])
    negatives = sorted(row[score_field] for row in rows if not row[label_field])
    if not positives or not negatives:
        return None
    # Rank-sum AUC with average ranks for ties.
    ranked = sorted(
        [(score, 1) for score in positives] + [(score, 0) for score in negatives],
        key=lambda item: item[0],
    )
    rank_sum = 0.0
    index = 0
    while index < len(ranked):
        end = index
        while end + 1 < len(ranked) and ranked[end + 1][0] == ranked[index][0]:
            end += 1
        average_rank = (index + end + 2) / 2.0  # 1-based average rank
        rank_sum += sum(average_rank for _, label in ranked[index : end + 1] if label == 1)
        index = end + 1
    n_pos = len(positives)
    n_neg = len(negatives)
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def bootstrap_auc(rows, score_field, replicates=500):
    values = []
    for _ in range(replicates):
        sample = [random.choice(rows) for _ in rows]
        value = auc(sample, score_field)
        if value is not None:
            values.append(value)
    if not values:
        return (None, None)
    values.sort()
    return (values[int(0.025 * len(values))], values[int(0.975 * len(values)) - 1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairwise", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--similarity", type=float, default=0.80)
    args = parser.parse_args()
    random.seed(20261124)
    rows = []
    for row in read_tsv(args.pairwise):
        try:
            dbd = float(row["dbd_identity"])
            protein = float(row["protein_identity"])
            similarity = float(row["motif_similarity"])
            overlap = float(row["motif_overlap_fraction"])
        except (ValueError, KeyError):
            continue
        if overlap < 0.50:
            continue
        rows.append(
            {
                "family": row["family"],
                "dbd_identity": dbd,
                "protein_identity": protein,
                "motif_similarity": similarity,
                "label": similarity >= args.similarity,
                "species_relation": row["species_relation"],
            }
        )
    overall = {
        "n_pairs": len(rows),
        "positives": sum(row["label"] for row in rows),
        "auc_dbd": auc(rows, "dbd_identity"),
        "auc_protein": auc(rows, "protein_identity"),
    }
    overall["auc_dbd_ci"] = bootstrap_auc(rows, "dbd_identity")
    overall["auc_protein_ci"] = bootstrap_auc(rows, "protein_identity")
    by_family = collections.defaultdict(list)
    for row in rows:
        by_family[row["family"]].append(row)
    family_rows = []
    for family, group in sorted(by_family.items()):
        if len(group) < 20:
            continue
        positives = sum(row["label"] for row in group)
        if positives == 0 or positives == len(group):
            continue
        family_rows.append(
            {
                "family": family,
                "n_pairs": len(group),
                "positives": positives,
                "auc_dbd": f"{auc(group, 'dbd_identity'):.4f}",
                "auc_protein": f"{auc(group, 'protein_identity'):.4f}",
                "delta_dbd_minus_protein": f"{auc(group, 'dbd_identity') - auc(group, 'protein_identity'):.4f}",
            }
        )
    write_tsv(args.out / "pairwise-predictor-auc.tsv", family_rows)
    bins = collections.defaultdict(lambda: {"n": 0, "positive": 0})
    for row in rows:
        key = (int(row["dbd_identity"] // 5) * 5, row["species_relation"])
        bins[key]["n"] += 1
        bins[key]["positive"] += int(row["label"])
    bin_rows = []
    for (identity, relation), values in sorted(bins.items()):
        if values["n"] < 10:
            continue
        bin_rows.append(
            {
                "dbd_identity_bin": identity,
                "species_relation": relation,
                "n_pairs": values["n"],
                "fraction_motif_similar_080": f"{values['positive'] / values['n']:.4f}",
            }
        )
    write_tsv(args.out / "pairwise-dbd-bins.tsv", bin_rows)
    (args.out / "pairwise-predictor-summary.json").write_text(
        json.dumps(overall, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(overall, indent=2))
    print("family rows", len(family_rows))


if __name__ == "__main__":
    main()
