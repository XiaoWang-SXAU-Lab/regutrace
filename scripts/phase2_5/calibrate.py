#!/usr/bin/env python3
"""Calibrate family-specific DBD identity/coverage thresholds.

The calibration task is donor-held-out motif transfer: for each query donor,
rank same-family reference donors by DBD BLAST bitscore and transfer the motif
set of the unique top donor. A prediction is correct when the predicted and
query motif sets share an operational specificity group (TOMTOM-style
correlation >= 0.90, >= 50% overlap) or an exact motif ID.

Threshold selection uses training folds only. Species-grouped cross-validation
is the primary leakage control; homology-cluster-grouped cross-validation is
reported as a secondary sensitivity analysis. Final deployment thresholds are
selected on all calibration donors with the same pre-registered criterion.
"""
import argparse
import collections
import csv
import json
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_tsv, write_tsv  # noqa: E402

IDENTITY_GRID = [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95]
COVERAGE_GRID = [0.5, 0.6, 0.7, 0.8, 0.9]
MAX_EVALUE = 1e-5
PRECISION_TARGET = 0.90
SIMILARITY_TARGET = 0.80
HIGH_SIMILARITY_TARGET = 0.90
MIN_TRAINING_COVERAGE = 0.20
MIN_COVERED = 5
HOMOLOGY_IDENTITY = 0.70
HOMOLOGY_COVERAGE = 0.70


def wilson(correct, total, z=1.96):
    if total == 0:
        return (0.0, 0.0, 0.0)
    p = correct / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return (p, max(0.0, center - half), min(1.0, center + half))


def load_alignments(path):
    by_query = collections.defaultdict(list)
    reverse = collections.defaultdict(list)
    for line in open(path, encoding="utf-8"):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 12:
            continue
        query, subject = fields[0], fields[1]
        if query == subject:
            continue
        bitscore = float(fields[2])
        identity = float(fields[3])
        qlen = float(fields[5])
        slen = float(fields[6])
        evalue = float(fields[7])
        qstart, qend, sstart, send = map(float, fields[8:12])
        coverage = min((qend - qstart + 1) / qlen, (send - sstart + 1) / slen)
        by_query[query].append(
            {
                "subject": subject,
                "bitscore": bitscore,
                "identity": identity,
                "coverage": coverage,
                "evalue": evalue,
            }
        )
        reverse[subject].append((query, bitscore))
    for query in by_query:
        by_query[query].sort(
            key=lambda row: (-row["bitscore"], -row["identity"], row["subject"])
        )
    return by_query, reverse


def load_motif_similarity(path):
    similar = {}
    for row in read_tsv(path):
        key = tuple(sorted((row["motif_a"], row["motif_b"])))
        similar[key] = (
            float(row["correlation"]),
            float(row["overlap_fraction"]),
        )
    return similar


def motif_match(motifs_a, motifs_b, pair_similarity, threshold=0.90):
    for a in motifs_a:
        for b in motifs_b:
            if a == b:
                return True
            correlation, overlap = pair_similarity.get(tuple(sorted((a, b))), (-1.0, 0.0))
            if correlation >= threshold and overlap >= 0.50:
                return True
    return False


def predict(query, references, alignments, reverse_index, thresholds, motifs_by_donor, pair_similarity):
    identity_min, coverage_min = thresholds
    candidates = [
        row
        for row in alignments.get(query, [])
        if row["subject"] in references
        and row["identity"] >= identity_min
        and row["coverage"] >= coverage_min
        and row["evalue"] <= MAX_EVALUE
    ]
    if not candidates:
        return {
            "query": query,
            "covered": False,
            "calibrated": False,
            "reciprocal": False,
            "ambiguous": False,
            "reason": "no_candidate",
            "top_donor": "",
            "predicted_motifs": [],
            "query_motifs": motifs_by_donor[query],
            "correct": False,
            "exact": False,
            "dbd_identity": "",
            "dbd_coverage": "",
        }
    top_score = candidates[0]["bitscore"]
    top = [row for row in candidates if abs(row["bitscore"] - top_score) < 1e-9]
    top_donors = sorted(row["subject"] for row in top)
    query_motifs = motifs_by_donor[query]
    if len(top_donors) > 1:
        predicted = sorted(
            {motif for donor in top_donors for motif in motifs_by_donor[donor]}
        )
        return {
            "query": query,
            "covered": True,
            "calibrated": False,
            "reciprocal": False,
            "ambiguous": True,
            "reason": "multiple_top_donors",
            "top_donor": ";".join(top_donors),
            "predicted_motifs": predicted,
            "query_motifs": query_motifs,
            "correct": False,
            "exact": False,
            "high_similarity": False,
            "dbd_identity": candidates[0]["identity"],
            "dbd_coverage": candidates[0]["coverage"],
        }
    donor = top_donors[0]
    predicted = motifs_by_donor[donor]
    reverse_best = max(
        (
            (candidate_query, score)
            for candidate_query, score in reverse_index.get(donor, [])
            if candidate_query == query or candidate_query in references
        ),
        key=lambda item: item[1],
        default=(None, -1.0),
    )[0]
    reciprocal = reverse_best == query
    exact = bool(set(predicted) & set(query_motifs))
    correct = reciprocal and motif_match(
        predicted, query_motifs, pair_similarity, SIMILARITY_TARGET
    )
    high = reciprocal and motif_match(
        predicted, query_motifs, pair_similarity, HIGH_SIMILARITY_TARGET
    )
    return {
        "query": query,
        "covered": True,
        "calibrated": reciprocal,
        "reciprocal": reciprocal,
        "ambiguous": False,
        "reason": "calibrated" if reciprocal else "not_reciprocal_dbd_best",
        "top_donor": donor,
        "predicted_motifs": predicted,
        "query_motifs": query_motifs,
        "correct": correct,
        "exact": exact,
        "dbd_identity": candidates[0]["identity"],
        "dbd_coverage": candidates[0]["coverage"],
    }


def summarize(predictions):
    covered = [row for row in predictions if row["covered"]]
    calibrated = [row for row in covered if row.get("calibrated", not row["ambiguous"])]
    correct = [row for row in calibrated if row["correct"]]
    exact = [row for row in calibrated if row["exact"]]
    high_rows = [row for row in calibrated if row.get("high_similarity")]
    ambiguous = [row for row in covered if row["ambiguous"]]
    non_reciprocal = [
        row for row in covered if row.get("reciprocal") is False and not row["ambiguous"]
    ]
    precision, low, high = wilson(len(correct), len(calibrated))
    coverage, cov_low, cov_high = wilson(len(calibrated), len(predictions))
    return {
        "n": len(predictions),
        "covered": len(covered),
        "calibrated": len(calibrated),
        "correct": len(correct),
        "exact": len(exact),
        "high_similarity_correct": len(high_rows),
        "high_similarity_rate": (len(high_rows) / len(calibrated)) if calibrated else 0.0,
        "ambiguous": len(ambiguous),
        "non_reciprocal": len(non_reciprocal),
        "precision": precision,
        "precision_ci_low": low,
        "precision_ci_high": high,
        "coverage": coverage,
        "coverage_ci_low": cov_low,
        "coverage_ci_high": cov_high,
    }


def loso_predictions(donors, donor_species, alignments, reverse_index, motifs_by_donor, pair_similarity, thresholds):
    folds = species_folds(donors, donor_species)
    if len(folds) < 2:
        return None
    predictions = []
    for fold in folds:
        train = set(donors) - set(fold["donors"])
        for query in fold["donors"]:
            predictions.append(
                predict(
                    query,
                    train - {query},
                    alignments,
                    reverse_index,
                    thresholds,
                    motifs_by_donor,
                    pair_similarity,
                )
            )
    return predictions


def select_threshold_loso(donors, donor_species, alignments, reverse_index, motifs_by_donor, pair_similarity):
    best = None
    for identity in IDENTITY_GRID:
        for coverage in COVERAGE_GRID:
            predictions = loso_predictions(
                donors,
                donor_species,
                alignments,
                reverse_index,
                motifs_by_donor,
                pair_similarity,
                (identity, coverage),
            )
            if predictions is None:
                continue
            summary = summarize(predictions)
            if summary["calibrated"] < max(
                MIN_COVERED, math.ceil(MIN_TRAINING_COVERAGE * len(donors))
            ):
                continue
            if summary["precision"] + 1e-12 < PRECISION_TARGET:
                continue
            rank = (summary["coverage"], identity, coverage)
            if best is None or rank > best[0]:
                best = (rank, identity, coverage, summary)
    return best


def species_folds(members, donor_species, max_folds=5):
    by_species = collections.defaultdict(list)
    for donor in members:
        by_species[donor_species[donor]].append(donor)
    species = sorted(by_species, key=lambda name: (-len(by_species[name]), name))
    fold_count = min(max_folds, len(species))
    if fold_count < 2:
        return []
    folds = [{"name": f"species_fold_{index + 1}", "species": set(), "donors": []} for index in range(fold_count)]
    for name in species:
        target = min(folds, key=lambda fold: (len(fold["donors"]), fold["name"]))
        target["species"].add(name)
        target["donors"].extend(by_species[name])
    return folds


def homology_clusters(members, alignments, max_folds=5):
    parent = {donor: donor for donor in members}

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a, b):
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    member_set = set(members)
    for donor in members:
        for row in alignments.get(donor, []):
            if (
                row["subject"] in member_set
                and row["identity"] >= HOMOLOGY_IDENTITY * 100
                and row["coverage"] >= HOMOLOGY_COVERAGE
            ):
                union(donor, row["subject"])
    clusters = collections.defaultdict(list)
    for donor in members:
        clusters[find(donor)].append(donor)
    clusters = sorted(clusters.values(), key=lambda group: (-len(group), sorted(group)[0]))
    fold_count = min(max_folds, len(clusters))
    if fold_count < 2:
        return []
    folds = [{"name": f"homology_fold_{index + 1}", "clusters": [], "donors": []} for index in range(fold_count)]
    for cluster in clusters:
        target = min(folds, key=lambda fold: (len(fold["donors"]), fold["name"]))
        target["clusters"].append(cluster)
        target["donors"].extend(cluster)
    return folds


def nested_species_cv(family, members, donor_species, alignments, reverse_index, motifs_by_donor, pair_similarity):
    outer = species_folds(members, donor_species)
    if len(outer) < 2:
        return [], []
    performance = []
    oof = []
    for fold in outer:
        train = [donor for donor in members if donor not in set(fold["donors"])]
        test = fold["donors"]
        selected = select_threshold_loso(
            train, donor_species, alignments, reverse_index, motifs_by_donor, pair_similarity
        )
        if selected is None:
            performance.append(
                {
                    "family": family,
                    "outer_fold": fold["name"],
                    "selected": "no",
                    "identity": "",
                    "coverage": "",
                    "train_n": len(train),
                    "test_n": len(test),
                    **{key: "" for key in summarize([])},
                }
            )
            continue
        _, identity, coverage, inner_summary = selected
        predictions = [
            predict(
                query,
                set(train) - {query},
                alignments,
                reverse_index,
                (identity, coverage),
                motifs_by_donor,
                pair_similarity,
            )
            for query in test
        ]
        test_summary = summarize(predictions)
        for query, prediction in zip(test, predictions):
            oof.append(
                {
                    "family": family,
                    "outer_fold": fold["name"],
                    "donor": query,
                    "species": donor_species[query],
                    **prediction,
                }
            )
        performance.append(
            {
                "family": family,
                "outer_fold": fold["name"],
                "selected": "yes",
                "identity": identity,
                "coverage": coverage,
                "train_n": len(train),
                "test_n": len(test),
                "inner_precision": f"{inner_summary['precision']:.4f}",
                "inner_coverage": f"{inner_summary['coverage']:.4f}",
                **test_summary,
            }
        )
    return performance, oof


def fixed_fold_eval(members, folds, donor_species, alignments, reverse_index, motifs_by_donor, pair_similarity, thresholds):
    if not folds:
        return summarize([])
    predictions = []
    for fold in folds:
        train = set(members) - set(fold["donors"])
        for query in fold["donors"]:
            predictions.append(
                predict(
                    query,
                    train - {query},
                    alignments,
                    reverse_index,
                    thresholds,
                    motifs_by_donor,
                    pair_similarity,
                )
            )
    return summarize(predictions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--dbd-blast", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--threshold-config", type=Path, required=True)
    parser.add_argument("--model-config", type=Path, default=Path("config/dbd-models/plant-v1.json"))
    args = parser.parse_args()
    donors = read_tsv(args.dataset / "calibration-donors.tsv")
    family_rows = {
        row["family"]: row
        for row in read_tsv(args.dataset / "calibration-family-summary.tsv")
    }
    pair_similarity = load_motif_similarity(args.dataset / "motif-pairwise-similarity.tsv")
    motifs_by_donor = {
        row["donor_id"]: [motif for motif in row["motif_ids"].split(";") if motif]
        for row in donors
    }
    donor_species = {row["donor_id"]: row["species"] for row in donors}
    donor_family = {row["donor_id"]: row["family"] for row in donors}
    alignments, reverse_index = load_alignments(args.dbd_blast)
    by_family = collections.defaultdict(list)
    for row in donors:
        by_family[row["family"]].append(row["donor_id"])

    performance_rows = []
    threshold_rows = []
    abstention_rows = []
    oof_rows = []
    sweep_rows = []
    for family, members in sorted(by_family.items()):
        members.sort()
        support = family_rows.get(family, {})
        if support.get("calibration_status") != "supported":
            abstention_rows.append(
                {
                    "family": family,
                    "donors": len(members),
                    "species": support.get("species", ""),
                    "reason": support.get("abstain_reason", "family_not_supported"),
                }
            )
            continue
        for identity in IDENTITY_GRID:
            for coverage in COVERAGE_GRID:
                sweep_predictions = loso_predictions(
                    members,
                    donor_species,
                    alignments,
                    reverse_index,
                    motifs_by_donor,
                    pair_similarity,
                    (identity, coverage),
                )
                if sweep_predictions is None:
                    continue
                sweep_summary = summarize(sweep_predictions)
                sweep_rows.append(
                    {
                        "family": family,
                        "min_identity": identity,
                        "min_coverage": coverage,
                        "calibrated": sweep_summary["calibrated"],
                        "correct": sweep_summary["correct"],
                        "precision": f"{sweep_summary['precision']:.4f}",
                        "precision_ci_low": f"{sweep_summary['precision_ci_low']:.4f}",
                        "coverage": f"{sweep_summary['coverage']:.4f}",
                        "coverage_ci_low": f"{sweep_summary['coverage_ci_low']:.4f}",
                        "ambiguous": sweep_summary["ambiguous"],
                    }
                )
        nested_perf, nested_oof = nested_species_cv(
            family, members, donor_species, alignments, reverse_index, motifs_by_donor, pair_similarity
        )
        performance_rows.extend(nested_perf)
        if nested_oof:
            oof_rows.extend(nested_oof)
        nested_summary = (
            summarize(
                [
                    {
                        "covered": row["covered"],
                        "calibrated": row.get("calibrated", not row["ambiguous"]),
                        "correct": row["correct"],
                        "exact": row["exact"],
                        "high_similarity": row.get("high_similarity", False),
                        "ambiguous": row["ambiguous"],
                    }
                    for row in nested_oof
                ]
            )
            if nested_oof
            else summarize([])
        )
        if not nested_oof:
            abstention_rows.append(
                {
                    "family": family,
                    "donors": len(members),
                    "species": support.get("species", ""),
                    "reason": "no_nested_species_calls",
                }
            )
            continue
        minimum_calls = max(5, math.ceil(0.10 * len(members)))
        if nested_summary["calibrated"] < minimum_calls:
            abstention_rows.append(
                {
                    "family": family,
                    "donors": len(members),
                    "species": support.get("species", ""),
                    "reason": f"insufficient_nested_calls(<{minimum_calls})",
                }
            )
            continue
        if (
            nested_summary["precision"] + 1e-12 < PRECISION_TARGET
            or nested_summary["precision_ci_low"] + 1e-12 < 0.70
        ):
            abstention_rows.append(
                {
                    "family": family,
                    "donors": len(members),
                    "species": support.get("species", ""),
                    "reason": "nested_validation_precision_below_target",
                }
            )
            continue
        selected = select_threshold_loso(
            members, donor_species, alignments, reverse_index, motifs_by_donor, pair_similarity
        )
        if selected is None:
            abstention_rows.append(
                {
                    "family": family,
                    "donors": len(members),
                    "species": support.get("species", ""),
                    "reason": "no_final_threshold_meets_criterion",
                }
            )
            continue
        _, identity, coverage, full_summary = selected
        homology = homology_clusters(members, alignments)
        homology_summary = fixed_fold_eval(
            members,
            homology,
            donor_species,
            alignments,
            reverse_index,
            motifs_by_donor,
            pair_similarity,
            (identity, coverage),
        )
        threshold_rows.append(
            {
                "family": family,
                "min_identity": identity,
                "min_coverage": coverage,
                "max_evalue": MAX_EVALUE,
                "calibration_status": "calibrated",
                "selection_criterion": (
                    f"nested_species_cv_precision>={PRECISION_TARGET},"
                    f"-ci_low>=0.70,nested_calibrated>={minimum_calls};"
                    f"final_threshold_selected_by_loso_coverage_max"
                ),
                "donors": len(members),
                "species": support.get("species", ""),
                "final_loso_covered": full_summary["covered"],
                "final_loso_calibrated": full_summary["calibrated"],
                "final_loso_correct": full_summary["correct"],
                "final_loso_precision": f"{full_summary['precision']:.4f}",
                "final_loso_precision_ci_low": f"{full_summary['precision_ci_low']:.4f}",
                "final_loso_coverage": f"{full_summary['coverage']:.4f}",
                "nested_oof_n": nested_summary["n"],
                "nested_oof_calibrated": nested_summary["calibrated"],
                "nested_oof_correct": nested_summary["correct"],
                "nested_oof_precision": f"{nested_summary['precision']:.4f}",
                "nested_oof_precision_ci_low": f"{nested_summary['precision_ci_low']:.4f}",
                "nested_oof_precision_ci_high": f"{nested_summary['precision_ci_high']:.4f}",
                "nested_oof_coverage": f"{nested_summary['coverage']:.4f}",
                "nested_oof_coverage_ci_low": f"{nested_summary['coverage_ci_low']:.4f}",
                "nested_oof_ambiguous": nested_summary["ambiguous"],
                "nested_oof_exact": nested_summary["exact"],
                "nested_oof_high_similarity_rate": f"{nested_summary['high_similarity_rate']:.4f}",
                "homology_precision": f"{homology_summary['precision']:.4f}",
                "homology_precision_ci_low": f"{homology_summary['precision_ci_low']:.4f}",
                "homology_coverage": f"{homology_summary['coverage']:.4f}",
                "homology_ambiguous": homology_summary["ambiguous"],
                "notes": "nested species-grouped CV for evaluation; final threshold selected on all donors",
            }
        )
    write_tsv(args.out / "calibration-performance.tsv", performance_rows)
    write_tsv(args.out / "calibration-thresholds.tsv", threshold_rows)
    write_tsv(args.out / "family-abstention-reasons.tsv", abstention_rows)
    write_tsv(args.out / "oof-predictions.tsv", oof_rows)
    write_tsv(args.out / "calibration-threshold-sweep.tsv", sweep_rows)

    calibrated = {row["family"]: row for row in threshold_rows}
    config = {
        "version": "plant-calibrated-v1",
        "calibration_status": "calibrated",
        "source": "Phase 2.5 donor-held-out DBD calibration; see docs/development/dbd-calibration-methodology.md",
        "motif_similarity_target": SIMILARITY_TARGET,
        "high_similarity_target": HIGH_SIMILARITY_TARGET,
        "families": {
            family: {
                "min_identity": row["min_identity"],
                "min_coverage": row["min_coverage"],
                "max_evalue": MAX_EVALUE,
                "source": "phase2.5 nested species-grouped cross-validation",
                "donors": row["donors"],
                "nested_oof_precision": row["nested_oof_precision"],
                "nested_oof_precision_ci_low": row["nested_oof_precision_ci_low"],
                "nested_oof_coverage": row["nested_oof_coverage"],
                "nested_oof_ambiguous": row["nested_oof_ambiguous"],
            }
            for family, row in calibrated.items()
        },
    }
    args.threshold_config.parent.mkdir(parents=True, exist_ok=True)
    args.threshold_config.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    model_config = json.loads(args.model_config.read_text(encoding="utf-8"))
    combined = {
        "version": "plant-calibrated-v1",
        "source": "Phase 2.5 calibration; DBD models from plant-dbd-models-v1",
        "calibration_status": "calibrated",
        "motif_similarity_target": SIMILARITY_TARGET,
        "high_similarity_target": HIGH_SIMILARITY_TARGET,
        "families": {
            family: {
                "models": model_config.get("families", {}).get(family, {}).get("models", []),
                "status": "calibrated",
                "min_identity": row["min_identity"],
                "min_coverage": row["min_coverage"],
                "max_evalue": MAX_EVALUE,
                "source": "phase2.5 nested species-grouped cross-validation",
                "nested_oof_precision": row["nested_oof_precision"],
                "nested_oof_coverage": row["nested_oof_coverage"],
            }
            for family, row in calibrated.items()
        },
    }
    args.config.parent.mkdir(parents=True, exist_ok=True)
    args.config.write_text(json.dumps(combined, indent=2) + "\n", encoding="utf-8")
    print(f"thresholds: {len(threshold_rows)} families")
    print(f"abstentions: {len(abstention_rows)} families")
    print(f"performance rows: {len(performance_rows)}")
    for row in threshold_rows:
        print(
            row["family"],
            row["min_identity"],
            row["min_coverage"],
            "nested OOF precision",
            row["nested_oof_precision"],
            "coverage",
            row["nested_oof_coverage"],
        )


if __name__ == "__main__":
    main()
