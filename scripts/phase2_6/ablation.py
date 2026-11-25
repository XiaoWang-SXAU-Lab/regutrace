#!/usr/bin/env python3
"""Phase 2.6 mapping-strategy ablation on the Phase 2.5 calibration set.

Methods:
A  full-length strict RBH (legacy thresholds, reciprocal)
B  DBD reciprocal best hit, family-specific thresholds
C  DBD nearest neighbor, no reciprocity, one global threshold
D  DBD nearest neighbor, no reciprocity, family-specific thresholds
E  DBD specificity-group-aware inference, no reciprocity

All methods use the same donor set and the same global species-grouped outer
folds. B-E are recalibrated independently inside each outer training fold.
"""
import argparse
import collections
import csv
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "phase2_5"))
from common import read_tsv, write_tsv  # noqa: E402
from calibrate import (  # noqa: E402
    load_alignments,
    load_motif_similarity,
    motif_match,
    wilson,
)

IDENTITY_GRID = [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95]
COVERAGE_GRID = [0.5, 0.6, 0.7, 0.8, 0.9]
MAX_EVALUE = 1e-5
STRICT_IDENTITY = 35.0
STRICT_COVERAGE = 0.60
PRECISION_TARGET = 0.90
MIN_CALLS = 5
MIN_TRAIN_FRACTION = 0.20
SIMILARITY_TARGET = 0.80
HIGH_SIMILARITY_TARGET = 0.90


def global_species_folds(donors, species_of, k=5):
    by_species = collections.defaultdict(list)
    for donor in donors:
        by_species[species_of[donor]].append(donor)
    species = sorted(by_species, key=lambda name: (-len(by_species[name]), name))
    fold_count = min(k, len(species))
    if fold_count < 2:
        return []
    folds = [
        {"name": f"species_fold_{index + 1}", "species": set(), "donors": []}
        for index in range(fold_count)
    ]
    for name in species:
        target = min(folds, key=lambda fold: (len(fold["donors"]), fold["name"]))
        target["species"].add(name)
        target["donors"].extend(by_species[name])
    return folds


def donor_random_folds(donors, k=5, seed=20261124):
    shuffled = sorted(donors)
    random.Random(seed).shuffle(shuffled)
    folds = [{"name": f"donor_fold_{index + 1}", "species": set(), "donors": []} for index in range(k)]
    for index, donor in enumerate(shuffled):
        folds[index % k]["donors"].append(donor)
    return [fold for fold in folds if fold["donors"]]


def global_homology_folds(donors, alignments, k=5, identity=70.0, coverage=0.70):
    parent = {donor: donor for donor in donors}

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a, b):
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    member_set = set(donors)
    for donor in donors:
        for row in alignments.get(donor, []):
            if (
                row["subject"] in member_set
                and row["identity"] >= identity
                and row["coverage"] >= coverage
            ):
                union(donor, row["subject"])
    clusters = collections.defaultdict(list)
    for donor in donors:
        clusters[find(donor)].append(donor)
    clusters = sorted(clusters.values(), key=lambda group: (-len(group), sorted(group)[0]))
    fold_count = min(k, len(clusters))
    if fold_count < 2:
        return []
    folds = [{"name": f"homology_fold_{index + 1}", "species": set(), "donors": []} for index in range(fold_count)]
    for cluster in clusters:
        target = min(folds, key=lambda fold: (len(fold["donors"]), fold["name"]))
        target["donors"].extend(cluster)
    return folds


def reverse_best(donor, query, panel, reverse_index):
    best = max(
        (
            (candidate, score)
            for candidate, score in reverse_index.get(donor, [])
            if candidate == query or candidate in panel
        ),
        key=lambda item: item[1],
        default=(None, -1.0),
    )
    return best[0]


def top_donors(query, references, alignments, family_of, thresholds, reciprocal, reverse_index):
    identity_min, coverage_min = thresholds
    rows = [
        row
        for row in alignments.get(query, [])
        if row["subject"] in references
        and family_of.get(row["subject"]) == family_of.get(query)
        and row["identity"] >= identity_min
        and row["coverage"] >= coverage_min
        and row["evalue"] <= MAX_EVALUE
    ]
    if not rows:
        return [], [], "no_candidate"
    top_score = rows[0]["bitscore"]
    if not reciprocal:
        top = sorted(
            {row["subject"] for row in rows if abs(row["bitscore"] - top_score) < 1e-9}
        )
        return top, rows, "ok" if len(top) == 1 else "multiple_top_donors"
    # Reciprocal methods: rank by forward bitscore first, then require the top
    # donor to return the query as its best hit. A lower reciprocal donor is not
    # substituted (same behaviour as the legacy strict map).
    top = sorted({row["subject"] for row in rows if abs(row["bitscore"] - top_score) < 1e-9})
    if len(top) > 1:
        return top, rows, "multiple_top_donors"
    donor = top[0]
    if reverse_best(donor, query, references, reverse_index) != query:
        return [], rows, "not_reciprocal"
    return [donor], rows, "ok"


def top_donors_legacy(query, references, alignments, family_of, thresholds, reverse_index):
    """Legacy strict ranking: best raw bitscore first, then eligibility.

    A lower-scoring eligible row is not substituted, matching the Phase 0
    audit finding and the original strict implementation.
    """
    identity_min, coverage_min = thresholds
    rows = [
        row
        for row in alignments.get(query, [])
        if row["subject"] in references
        and family_of.get(row["subject"]) == family_of.get(query)
    ]
    if not rows:
        return [], "no_candidate"
    top_score = rows[0]["bitscore"]
    top = [row for row in rows if abs(row["bitscore"] - top_score) < 1e-9]
    donors = sorted({row["subject"] for row in top})
    if len(donors) > 1:
        return donors, "multiple_top_donors"
    row = top[0]
    if (
        row["identity"] < identity_min
        or row["coverage"] < coverage_min
        or row["evalue"] > MAX_EVALUE
    ):
        return [], "below_identity_or_coverage"
    donor = donors[0]
    if reverse_best(donor, query, references, reverse_index) != query:
        return [], "not_reciprocal"
    return [donor], "ok"


def predict_one(
    method,
    query,
    references,
    alignments,
    reverse_index,
    family_of,
    motifs_by_donor,
    group_of,
    pair_similarity,
    thresholds,
    reciprocal,
):
    query_motifs = motifs_by_donor[query]
    if method == "E":
        # All eligible candidates, not just the top donor. A call requires all
        # candidate motifs to collapse to one specificity group.
        identity_min, coverage_min = thresholds
        candidates = [
            row
            for row in alignments.get(query, [])
            if row["subject"] in references
            and family_of.get(row["subject"]) == family_of.get(query)
            and row["identity"] >= identity_min
            and row["coverage"] >= coverage_min
            and row["evalue"] <= MAX_EVALUE
        ]
        if not candidates:
            return result_row(method, query, family_of, query_motifs, "no_call", "", [], pair_similarity, "no_candidate")
        groups = collections.defaultdict(set)
        missing = False
        for row in candidates:
            for motif in motifs_by_donor[row["subject"]]:
                group = group_of.get(motif)
                if group is None:
                    missing = True
                    continue
                groups[group].add(motif)
        if missing or len(groups) != 1:
            return result_row(method, query, family_of, query_motifs, "ambiguous", ";".join(sorted({row["subject"] for row in candidates})), [], pair_similarity, "multiple_or_missing_specificity_groups")
        group, motifs = next(iter(groups.items()))
        return result_row(method, query, family_of, query_motifs, "called", ";".join(sorted({row["subject"] for row in candidates})), sorted(motifs), pair_similarity, f"specificity_group={group}")
    if method == "A":
        donors, reason = top_donors_legacy(
            query, references, alignments, family_of, thresholds, reverse_index
        )
        rows = []
    else:
        donors, rows, reason = top_donors(
            query, references, alignments, family_of, thresholds, reciprocal, reverse_index
        )
    if len(donors) > 1:
        return result_row(method, query, family_of, query_motifs, "ambiguous", ";".join(donors), [], pair_similarity, reason)
    if not donors:
        return result_row(method, query, family_of, query_motifs, "no_call", "", [], pair_similarity, reason)
    donor = donors[0]
    return result_row(method, query, family_of, query_motifs, "called", donor, motifs_by_donor[donor], pair_similarity, reason)


def result_row(method, query, family_of, query_motifs, status, donor, predicted, pair_similarity, reason):
    correct = False
    high = False
    exact = False
    similarity = ""
    if status == "called":
        correct = motif_match(predicted, query_motifs, pair_similarity, SIMILARITY_TARGET)
        high = motif_match(predicted, query_motifs, pair_similarity, HIGH_SIMILARITY_TARGET)
        exact = bool(set(predicted) & set(query_motifs))
        values = []
        for motif in predicted:
            for query_motif in query_motifs:
                if motif == query_motif:
                    values.append(1.0)
                elif (motif, query_motif) in pair_similarity:
                    values.append(pair_similarity[(motif, query_motif)][0])
                elif (query_motif, motif) in pair_similarity:
                    values.append(pair_similarity[(query_motif, motif)][0])
        similarity = max(values) if values else ""
    return {
        "method": method,
        "tf_id": query,
        "family": family_of.get(query, ""),
        "species": "",
        "status": status,
        "reason": reason,
        "top_donor": donor,
        "predicted_motifs": ";".join(predicted),
        "query_motifs": ";".join(query_motifs),
        "correct": correct,
        "high_similarity": high,
        "exact_motif": exact,
        "max_motif_similarity": similarity,
    }


def summarize(predictions):
    total = len(predictions)
    called = [row for row in predictions if row["status"] == "called"]
    correct = [row for row in called if row["correct"]]
    high = [row for row in called if row["high_similarity"]]
    ambiguous = [row for row in predictions if row["status"] == "ambiguous"]
    no_call = [row for row in predictions if row["status"] == "no_call"]
    precision, p_low, p_high = wilson(len(correct), len(called))
    coverage, c_low, c_high = wilson(len(called), total)
    recall, r_low, r_high = wilson(len(correct), total)
    donors = set()
    donor_counts = collections.Counter()
    for row in called:
        for donor in row["top_donor"].split(";"):
            if donor:
                donors.add(donor)
                donor_counts[donor] += 1
    similarities = [row["max_motif_similarity"] for row in called if row["max_motif_similarity"] != ""]
    return {
        "total": total,
        "accepted_calls": len(called),
        "correct_calls": len(correct),
        "high_similarity_calls": len(high),
        "ambiguous": len(ambiguous),
        "abstained": len(no_call),
        "precision": precision,
        "precision_ci_low": p_low,
        "precision_ci_high": p_high,
        "coverage": coverage,
        "coverage_ci_low": c_low,
        "coverage_ci_high": c_high,
        "recall": recall,
        "recall_ci_low": r_low,
        "recall_ci_high": r_high,
        "unique_donors_used": len(donors),
        "max_queries_per_donor": max(donor_counts.values()) if donor_counts else 0,
        "ambiguous_rate": len(ambiguous) / total if total else 0.0,
        "abstention_rate": len(no_call) / total if total else 0.0,
        "median_motif_similarity": statistics.median(similarities) if similarities else "",
        "mean_motif_similarity": statistics.mean(similarities) if similarities else "",
    }


def select_family_thresholds(
    method,
    train,
    inner_folds,
    family_of,
    alignments,
    reverse_index,
    motifs_by_donor,
    group_of,
    pair_similarity,
    family,
):
    family_train = [donor for donor in train if family_of.get(donor) == family]
    if len(family_train) < MIN_CALLS:
        return None
    best = None
    for identity in IDENTITY_GRID:
        for coverage in COVERAGE_GRID:
            predictions = []
            for fold in inner_folds:
                inner_train = [donor for donor in train if donor not in set(fold["donors"])]
                for query in fold["donors"]:
                    if family_of.get(query) != family:
                        continue
                    predictions.append(
                        predict_one(
                            method,
                            query,
                            set(inner_train),
                            alignments,
                            reverse_index,
                            family_of,
                            motifs_by_donor,
                            group_of,
                            pair_similarity,
                            (identity, coverage),
                            False,
                        )
                    )
            summary = summarize(predictions)
            if summary["accepted_calls"] < max(MIN_CALLS, math.ceil(MIN_TRAIN_FRACTION * len(family_train))):
                continue
            if summary["precision"] + 1e-12 < PRECISION_TARGET:
                continue
            rank = (summary["coverage"], identity, coverage)
            if best is None or rank > best[0]:
                best = (rank, identity, coverage, summary)
    return best


def select_global_threshold(
    method,
    train,
    inner_folds,
    family_of,
    alignments,
    reverse_index,
    motifs_by_donor,
    group_of,
    pair_similarity,
):
    if len(train) < MIN_CALLS:
        return None
    best = None
    for identity in IDENTITY_GRID:
        for coverage in COVERAGE_GRID:
            predictions = []
            for fold in inner_folds:
                inner_train = [donor for donor in train if donor not in set(fold["donors"])]
                for query in fold["donors"]:
                    predictions.append(
                        predict_one(
                            method,
                            query,
                            set(inner_train),
                            alignments,
                            reverse_index,
                            family_of,
                            motifs_by_donor,
                            group_of,
                            pair_similarity,
                            (identity, coverage),
                            False,
                        )
                    )
            summary = summarize(predictions)
            if summary["accepted_calls"] < max(MIN_CALLS, math.ceil(MIN_TRAIN_FRACTION * len(train))):
                continue
            if summary["precision"] + 1e-12 < PRECISION_TARGET:
                continue
            rank = (summary["coverage"], identity, coverage)
            if best is None or rank > best[0]:
                best = (rank, identity, coverage, summary)
    return best


def fixed_loso_sweep(method, donors, species_of, family_of, alignments, reverse_index, motifs_by_donor, group_of, pair_similarity):
    """Fixed-threshold leave-one-species-out curve for methods B-E."""
    rows = []
    families = sorted(set(family_of[donor] for donor in donors))
    for family in families:
        members = [donor for donor in donors if family_of[donor] == family]
        folds = global_species_folds(members, species_of, k=5)
        if len(folds) < 2:
            continue
        for identity in IDENTITY_GRID:
            for coverage in COVERAGE_GRID:
                predictions = []
                for fold in folds:
                    train = [donor for donor in members if donor not in set(fold["donors"])]
                    for query in fold["donors"]:
                        predictions.append(
                            predict_one(
                                method,
                                query,
                                set(train),
                                alignments,
                                reverse_index,
                                family_of,
                                motifs_by_donor,
                                group_of,
                                pair_similarity,
                                (identity, coverage),
                                method == "B",
                            )
                        )
                summary = summarize(predictions)
                rows.append(
                    {
                        "method": method,
                        "family": family,
                        "identity": identity,
                        "coverage": coverage,
                        "accepted_calls": summary["accepted_calls"],
                        "correct_calls": summary["correct_calls"],
                        "precision": f"{summary['precision']:.4f}",
                        "coverage_rate": f"{summary['coverage']:.4f}",
                        "ambiguous": summary["ambiguous"],
                    }
                )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--dbd-blast", type=Path, required=True)
    parser.add_argument("--protein-blast", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--fold-scheme", choices=["species", "donor", "homology"], default="species")
    parser.add_argument("--seed", type=int, default=20261124)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    donors = read_tsv(args.dataset / "calibration-donors.tsv")
    motifs_by_donor = {
        row["donor_id"]: [motif for motif in row["motif_ids"].split(";") if motif]
        for row in donors
    }
    family_of = {row["donor_id"]: row["family"] for row in donors}
    species_of = {row["donor_id"]: row["species"] for row in donors}
    group_of = {
        row["motif_id"]: row["group_090"]
        for row in read_tsv(args.dataset / "motif-specificity-groups.tsv")
    }
    pair_similarity = load_motif_similarity(args.dataset / "motif-pairwise-similarity.tsv")
    dbd_alignments, dbd_reverse = load_alignments(args.dbd_blast)
    protein_alignments, protein_reverse = load_alignments(args.protein_blast)
    all_donors = [row["donor_id"] for row in donors]
    def make_folds(donor_list):
        if args.fold_scheme == "species":
            return global_species_folds(donor_list, species_of, k=args.folds)
        if args.fold_scheme == "donor":
            return donor_random_folds(donor_list, k=args.folds, seed=args.seed)
        return global_homology_folds(donor_list, dbd_alignments, k=args.folds)

    outer = make_folds(all_donors)
    if len(outer) < 2:
        raise SystemExit(f"not enough {args.fold_scheme} folds")

    all_predictions = []
    threshold_rows = []
    runtimes = {}
    methods = ["A", "B", "C", "D", "E"]
    for method in methods:
        start = time.perf_counter()
        for outer_fold in outer:
            test = outer_fold["donors"]
            train = [donor for donor in all_donors if donor not in set(test)]
            inner = make_folds(train)
            if method == "A":
                thresholds = (STRICT_IDENTITY, STRICT_COVERAGE)
                alignments = protein_alignments
                reverse = protein_reverse
                reciprocal = True
            elif method == "C":
                selected = select_global_threshold(
                    method,
                    train,
                    inner,
                    family_of,
                    dbd_alignments,
                    dbd_reverse,
                    motifs_by_donor,
                    group_of,
                    pair_similarity,
                )
                if selected is None:
                    thresholds = None
                else:
                    _, identity, coverage, _ = selected
                    thresholds = (identity, coverage)
                    threshold_rows.append(
                        {
                            "method": method,
                            "outer_fold": outer_fold["name"],
                            "family": "GLOBAL",
                            "identity": identity,
                            "coverage": coverage,
                        }
                    )
                alignments = dbd_alignments
                reverse = dbd_reverse
                reciprocal = False
            else:
                selected_by_family = {}
                for family in sorted(set(family_of[donor] for donor in train)):
                    selected = select_family_thresholds(
                        method,
                        train,
                        inner,
                        family_of,
                        dbd_alignments,
                        dbd_reverse,
                        motifs_by_donor,
                        group_of,
                        pair_similarity,
                        family,
                    )
                    if selected is not None:
                        _, identity, coverage, _ = selected
                        selected_by_family[family] = (identity, coverage)
                        threshold_rows.append(
                            {
                                "method": method,
                                "outer_fold": outer_fold["name"],
                                "family": family,
                                "identity": identity,
                                "coverage": coverage,
                            }
                        )
                alignments = dbd_alignments
                reverse = dbd_reverse
                reciprocal = method == "B"
            for query in test:
                family = family_of[query]
                if method in {"B", "D", "E"}:
                    threshold = selected_by_family.get(family)
                    if threshold is None:
                        row = result_row(method, query, family_of, motifs_by_donor[query], "no_call", "", [], pair_similarity, "no_family_threshold")
                    else:
                        row = predict_one(
                            method,
                            query,
                            set(train),
                            alignments,
                            reverse,
                            family_of,
                            motifs_by_donor,
                            group_of,
                            pair_similarity,
                            threshold,
                            reciprocal,
                        )
                elif method == "C" and thresholds is None:
                    row = result_row(method, query, family_of, motifs_by_donor[query], "no_call", "", [], pair_similarity, "no_global_threshold")
                else:
                    row = predict_one(
                        method,
                        query,
                        set(train),
                        alignments,
                        reverse,
                        family_of,
                        motifs_by_donor,
                        group_of,
                        pair_similarity,
                        thresholds,
                        reciprocal,
                    )
                row["species"] = species_of[query]
                row["outer_fold"] = outer_fold["name"]
                all_predictions.append(row)
        runtimes[method] = time.perf_counter() - start

    write_tsv(args.out / "ablation-oof.tsv", all_predictions)
    method_rows = []
    for method in methods:
        predictions = [row for row in all_predictions if row["method"] == method]
        summary = summarize(predictions)
        method_rows.append(
            {
                "method": method,
                **summary,
                "runtime_seconds": f"{runtimes.get(method, 0.0):.2f}",
            }
        )
    write_tsv(args.out / "ablation-methods.tsv", method_rows)

    family_rows = []
    for method in methods:
        for family in sorted(set(family_of.values())):
            predictions = [
                row
                for row in all_predictions
                if row["method"] == method and row["family"] == family
            ]
            if not predictions:
                continue
            summary = summarize(predictions)
            family_rows.append({"method": method, "family": family, **summary})
    write_tsv(args.out / "ablation-family.tsv", family_rows)
    write_tsv(args.out / "ablation-selected-thresholds.tsv", threshold_rows)

    similarity_rows = []
    for method in methods:
        bins = collections.Counter()
        for row in all_predictions:
            if row["method"] != method or row["status"] != "called" or row["max_motif_similarity"] == "":
                continue
            value = float(row["max_motif_similarity"])
            bins[min(9, max(0, int(value * 10)))] += 1
        for index in range(10):
            similarity_rows.append(
                {
                    "method": method,
                    "bin": f"{index / 10:.1f}-{(index + 1) / 10:.1f}",
                    "calls": bins.get(index, 0),
                }
            )
    write_tsv(args.out / "ablation-similarity.tsv", similarity_rows)

    sweep_rows = []
    for method in ("B", "C", "D", "E"):
        if method == "C":
            # Global threshold curve: one threshold, all families.
            for identity in IDENTITY_GRID:
                for coverage in COVERAGE_GRID:
                    predictions = []
                    folds = global_species_folds(all_donors, species_of, k=args.folds)
                    for fold in folds:
                        train = [donor for donor in all_donors if donor not in set(fold["donors"])]
                        for query in fold["donors"]:
                            predictions.append(
                                predict_one(
                                    method,
                                    query,
                                    set(train),
                                    dbd_alignments,
                                    dbd_reverse,
                                    family_of,
                                    motifs_by_donor,
                                    group_of,
                                    pair_similarity,
                                    (identity, coverage),
                                    False,
                                )
                            )
                    summary = summarize(predictions)
                    sweep_rows.append(
                        {
                            "method": method,
                            "family": "GLOBAL",
                            "identity": identity,
                            "coverage": coverage,
                            "accepted_calls": summary["accepted_calls"],
                            "correct_calls": summary["correct_calls"],
                            "precision": f"{summary['precision']:.4f}",
                            "coverage_rate": f"{summary['coverage']:.4f}",
                            "ambiguous": summary["ambiguous"],
                        }
                    )
        else:
            sweep_rows.extend(
                fixed_loso_sweep(
                    method,
                    all_donors,
                    species_of,
                    family_of,
                    dbd_alignments,
                    dbd_reverse,
                    motifs_by_donor,
                    group_of,
                    pair_similarity,
                )
            )
    write_tsv(args.out / "ablation-sweep.tsv", sweep_rows)
    (args.out / "ablation-runtime.json").write_text(json.dumps(runtimes, indent=2) + "\n", encoding="utf-8")
    print("methods:")
    for row in method_rows:
        print(row)
    print("family rows", len(family_rows), "sweep rows", len(sweep_rows))


if __name__ == "__main__":
    main()
