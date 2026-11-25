#!/usr/bin/env python3
"""Phase 2.6 calibration-methodology audit.

Produces validation-setting comparisons, parameter sensitivity tables, leakage
counts, and pair/donor/homology-cluster bootstrap intervals.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "phase2_5"))
from common import read_fasta, read_tsv, write_tsv  # noqa: E402
from calibrate import load_alignments, load_motif_similarity, wilson  # noqa: E402
from analyze_pairwise import auc  # noqa: E402

SIMILARITY_THRESHOLDS = [0.70, 0.75, 0.80, 0.85, 0.90]
OVERLAP_THRESHOLDS = [0.30, 0.50, 0.70]


def load_oof(path):
    rows = []
    for row in read_tsv(path):
        row["correct"] = row["correct"] == "True"
        row["high_similarity"] = row["high_similarity"] == "True"
        rows.append(row)
    return rows


def summarize_oof(rows):
    called = [row for row in rows if row["status"] == "called"]
    correct = [row for row in called if row["correct"]]
    ambiguous = [row for row in rows if row["status"] == "ambiguous"]
    precision, p_low, p_high = wilson(len(correct), len(called))
    coverage, c_low, c_high = wilson(len(called), len(rows))
    recall, r_low, r_high = wilson(len(correct), len(rows))
    return {
        "total": len(rows),
        "accepted_calls": len(called),
        "correct_calls": len(correct),
        "ambiguous": len(ambiguous),
        "abstained": len(rows) - len(called) - len(ambiguous),
        "precision": precision,
        "precision_ci_low": p_low,
        "coverage": coverage,
        "coverage_ci_low": c_low,
        "recall": recall,
        "recall_ci_low": r_low,
        "ambiguous_rate": len(ambiguous) / len(rows) if rows else 0.0,
    }


def motif_match_at(predicted, query, pair_similarity, correlation_threshold, overlap_threshold):
    for a in predicted:
        for b in query:
            if a == b:
                return True
            correlation, overlap = pair_similarity.get(tuple(sorted((a, b))), (-1.0, 0.0))
            if correlation >= correlation_threshold and overlap >= overlap_threshold:
                return True
    return False


def homology_cluster_map(donors, alignments, identity=70.0, coverage=0.70):
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
    mapping = {}
    for index, (_, members) in enumerate(sorted(clusters.items(), key=lambda item: (-len(item[1]), sorted(item[1])[0])), 1):
        for donor in members:
            mapping[donor] = f"H{index:04d}"
    return mapping


def bootstrap_auc(rows, score_field, mode, donor_cluster, replicates=500, seed=20261124):
    rng = random.Random(seed)
    pairs_by_donor = collections.defaultdict(list)
    for row in rows:
        pairs_by_donor[row["donor_a"]].append(row)
        pairs_by_donor[row["donor_b"]].append(row)
    donors = sorted(pairs_by_donor)
    clusters = collections.defaultdict(list)
    for donor in donors:
        clusters[donor_cluster.get(donor, donor)].append(donor)
    cluster_names = sorted(clusters)
    values = []
    for _ in range(replicates):
        if mode == "pair":
            sample = [rng.choice(rows) for _ in rows]
        elif mode == "donor":
            sampled = {rng.choice(donors) for _ in range(len(donors))}
            sample = [
                row
                for row in rows
                if row["donor_a"] in sampled and row["donor_b"] in sampled
            ]
        else:
            sampled_clusters = [rng.choice(cluster_names) for _ in range(len(cluster_names))]
            sampled_donors = set()
            for cluster in sampled_clusters:
                sampled_donors.update(clusters[cluster])
            sample = [
                row
                for row in rows
                if row["donor_a"] in sampled_donors and row["donor_b"] in sampled_donors
            ]
        normalized = [
            {"score": float(row[score_field]), "label": row["label"]}
            for row in sample
        ]
        value = auc(normalized, "score")
        if value is not None:
            values.append(value)
    values.sort()
    if not values:
        return (None, None, None)
    return (
        values[int(0.025 * len(values))],
        values[int(0.975 * len(values)) - 1],
        statistics.mean(values),
    )



def cluster_bootstrap_metrics(rows, cluster_of, replicates=500, seed=20261124):
    rng = random.Random(seed)
    by_cluster = collections.defaultdict(list)
    for row in rows:
        by_cluster[cluster_of.get(row["tf_id"], row["tf_id"])].append(row)
    names = sorted(by_cluster)
    precision_values = []
    coverage_values = []
    for _ in range(replicates):
        sampled = [rng.choice(names) for _ in names]
        sample = []
        for cluster in sampled:
            sample.extend(by_cluster[cluster])
        called = [row for row in sample if row["status"] == "called"]
        correct = [row for row in called if row["correct"]]
        if called:
            precision_values.append(len(correct) / len(called))
        if sample:
            coverage_values.append(len(called) / len(sample))
    def interval(values):
        if not values:
            return ("", "", "")
        values.sort()
        return (
            f"{statistics.mean(values):.4f}",
            f"{values[int(0.025 * len(values))]:.4f}",
            f"{values[int(0.975 * len(values)) - 1]:.4f}",
        )
    return interval(precision_values), interval(coverage_values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--ablation-species", type=Path, required=True)
    parser.add_argument("--ablation-donor", type=Path, required=True)
    parser.add_argument("--ablation-homology", type=Path, required=True)
    parser.add_argument("--dbd-blast", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    donors = read_tsv(args.dataset / "calibration-donors.tsv")
    donor_ids = [row["donor_id"] for row in donors]
    species_of = {row["donor_id"]: row["species"] for row in donors}
    family_of = {row["donor_id"]: row["family"] for row in donors}
    motif_rows = {row["motif_id"]: row for row in read_tsv(args.dataset / "calibration-motifs.tsv")}
    pair_similarity = load_motif_similarity(args.dataset / "motif-pairwise-similarity.tsv")
    alignments, _ = load_alignments(args.dbd_blast)
    cluster_of = homology_cluster_map(donor_ids, alignments)

    # 1. Validation-setting comparison.
    validation_rows = []
    for scheme, path in (
        ("donor_held_out", args.ablation_donor),
        ("homology_cluster_held_out", args.ablation_homology),
        ("species_held_out", args.ablation_species),
    ):
        oof = load_oof(path / "ablation-oof.tsv")
        for method in ("A", "B", "C", "D", "E"):
            rows = [row for row in oof if row["method"] == method]
            summary = summarize_oof(rows)
            validation_rows.append({"fold_scheme": scheme, "method": method, **summary})
    write_tsv(args.out / "validation-setting-comparison.tsv", validation_rows)
    species_oof = load_oof(args.ablation_species / "ablation-oof.tsv")
    cluster_rows = []
    for method in ("A", "B", "D", "E"):
        rows = [row for row in species_oof if row["method"] == method]
        precision, coverage = cluster_bootstrap_metrics(rows, cluster_of)
        cluster_rows.append(
            {
                "method": method,
                "precision_mean": precision[0],
                "precision_ci_low": precision[1],
                "precision_ci_high": precision[2],
                "coverage_mean": coverage[0],
                "coverage_ci_low": coverage[1],
                "coverage_ci_high": coverage[2],
                "bootstrap_unit": "homology_cluster",
                "replicates": 500,
            }
        )
    write_tsv(args.out / "method-cluster-bootstrap.tsv", cluster_rows)

    # 2. Motif similarity and overlap sensitivity (species-held-out OOF).
    sensitivity_rows = []
    for method in ("A", "B", "D", "E"):
        rows = [row for row in species_oof if row["method"] == method and row["status"] == "called"]
        for correlation_threshold in SIMILARITY_THRESHOLDS:
            for overlap_threshold in OVERLAP_THRESHOLDS:
                correct = sum(
                    motif_match_at(
                        row["predicted_motifs"].split(";"),
                        row["query_motifs"].split(";"),
                        pair_similarity,
                        correlation_threshold,
                        overlap_threshold,
                    )
                    for row in rows
                    if row["predicted_motifs"]
                )
                precision, low, high = wilson(correct, len(rows))
                sensitivity_rows.append(
                    {
                        "method": method,
                        "motif_similarity_threshold": correlation_threshold,
                        "min_overlap": overlap_threshold,
                        "called": len(rows),
                        "correct": correct,
                        "precision": f"{precision:.4f}",
                        "precision_ci_low": f"{low:.4f}",
                    }
                )
    write_tsv(args.out / "motif-similarity-sensitivity.tsv", sensitivity_rows)

    # 3. Information-content stratification.
    ic_rows = []
    for method in ("A", "B", "D", "E"):
        values = []
        for row in species_oof:
            if row["method"] != method or row["status"] != "called":
                continue
            predicted = [motif for motif in row["predicted_motifs"].split(";") if motif]
            ics = [float(motif_rows[motif]["information_content"]) for motif in predicted if motif in motif_rows]
            if ics:
                values.append((max(ics), 1 if row["correct"] else 0))
        if len(values) < 10:
            continue
        values.sort()
        quartiles = [values[int(len(values) * q / 4)] [0] for q in (1, 2, 3)]
        for index in range(4):
            if index < 3:
                group = [value for value in values if (index == 0 or value[0] >= quartiles[index - 1]) and value[0] < quartiles[index]]
            else:
                group = [value for value in values if value[0] >= quartiles[-1]]
            if not group:
                continue
            correct = sum(value[1] for value in group)
            precision, low, _ = wilson(correct, len(group))
            ic_rows.append(
                {
                    "method": method,
                    "ic_quartile": index + 1,
                    "ic_range": f"{group[0][0]:.2f}-{group[-1][0]:.2f}",
                    "calls": len(group),
                    "correct": correct,
                    "precision": f"{precision:.4f}",
                    "precision_ci_low": f"{low:.4f}",
                }
            )
    write_tsv(args.out / "information-content-stratification.tsv", ic_rows)

    # 4. Representative-motif selection sensitivity.
    representative = {}
    for donor, motifs in {
        row["donor_id"]: [motif for motif in row["motif_ids"].split(";") if motif]
        for row in donors
    }.items():
        values = [(float(motif_rows[motif]["information_content"]), motif) for motif in motifs]
        representative[donor] = max(values)[1] if values else ""
    representative_rows = []
    for method in ("A", "B", "D", "E"):
        calls = [row for row in species_oof if row["method"] == method and row["status"] == "called"]
        all_correct = 0
        rep_correct = 0
        agree = 0
        for row in calls:
            query = [motif for motif in row["query_motifs"].split(";") if motif]
            predicted = [motif for motif in row["predicted_motifs"].split(";") if motif]
            all_value = motif_match_at(predicted, query, pair_similarity, 0.80, 0.50)
            top_donor = row["top_donor"].split(";")[0] if row["top_donor"] else ""
            rep = representative.get(top_donor, "")
            rep_value = bool(rep) and motif_match_at([rep], query, pair_similarity, 0.80, 0.50)
            all_correct += int(all_value)
            rep_correct += int(rep_value)
            agree += int(all_value == rep_value)
        if calls:
            representative_rows.append(
                {
                    "method": method,
                    "calls": len(calls),
                    "all_motifs_correct": all_correct,
                    "representative_correct": rep_correct,
                    "label_agreement": agree,
                    "label_agreement_rate": f"{agree / len(calls):.4f}",
                }
            )
    write_tsv(args.out / "representative-motif-sensitivity.tsv", representative_rows)

    # 5. Family support grid.
    support_rows = []
    for min_donors in (10, 15, 20, 30):
        for min_species in (2, 3, 4):
            for min_groups in (1, 2, 3):
                supported = [
                    row
                    for row in read_tsv(args.dataset / "calibration-family-summary.tsv")
                    if int(row["donors"]) >= min_donors
                    and int(row["species"]) >= min_species
                    and int(row["exact_groups"]) >= min_groups
                ]
                support_rows.append(
                    {
                        "min_donors": min_donors,
                        "min_species": min_species,
                        "min_exact_groups": min_groups,
                        "supported_families": len(supported),
                        "supported_donors": sum(int(row["donors"]) for row in supported),
                        "families": ";".join(row["family"] for row in supported),
                    }
                )
    write_tsv(args.out / "family-support-grid.tsv", support_rows)

    # 6. Leakage audit.
    sequences = {name: sequence for name, sequence in read_fasta(args.dataset / "calibration-proteins.faa")}
    sequence_hashes = collections.Counter(hashlib.sha256(sequence.encode()).hexdigest() for sequence in sequences.values())
    duplicate_sequences = sum(count - 1 for count in sequence_hashes.values() if count > 1)
    pairwise = [
        row
        for row in read_tsv(args.dataset / "dbd-motif-pairwise.tsv")
        if row["motif_similarity"]
        and row["dbd_identity"]
        and row["protein_identity"]
        and float(row["motif_overlap_fraction"]) >= 0.50
    ]
    same_species = sum(row["species_relation"] == "same_species" for row in pairwise)
    same_genus = sum(row["species_relation"] == "same_genus" for row in pairwise)
    exact_duplicate_pairs = sum(row["exact_motif_duplicate"] == "yes" for row in pairwise)
    homology_pairs = sum(
        1
        for row in pairwise
        if row["donor_a"] in cluster_of
        and cluster_of.get(row["donor_a"]) == cluster_of.get(row["donor_b"])
    )
    motif_sources = collections.Counter(row["source_releases"] for row in read_tsv(args.dataset / "calibration-motifs.tsv"))
    leakage = {
        "calibration_donors": len(donor_ids),
        "unique_protein_sequences": len(sequence_hashes),
        "duplicate_protein_sequences": duplicate_sequences,
        "homology_clusters_70_70": len(set(cluster_of.values())),
        "largest_homology_cluster": max(collections.Counter(cluster_of.values()).values()),
        "pairwise_rows_with_motif_similarity": len(pairwise),
        "same_species_pairs": same_species,
        "same_genus_pairs": same_genus,
        "different_genus_pairs": len(pairwise) - same_species - same_genus,
        "exact_duplicate_motif_pairs": exact_duplicate_pairs,
        "homology_pairs": homology_pairs,
        "motif_source_release_combinations": dict(motif_sources),
        "shared_2024_2026_motifs": motif_sources.get("JASPAR2024;JASPAR2026", 0),
    }
    write_tsv(args.out / "leakage-audit.tsv", [{"metric": key, "value": value} for key, value in leakage.items()])
    (args.out / "leakage-audit.json").write_text(json.dumps(leakage, indent=2) + "\n", encoding="utf-8")

    # 7. Pairwise bootstrap comparison.
    bootstrap_rows = []
    for row in pairwise:
        row["label"] = float(row["motif_similarity"]) >= 0.80
    for score_field in ("dbd_identity", "protein_identity"):
        overall = auc(
            [
                {
                    "score": float(row[score_field]),
                    "label": float(row["motif_similarity"]) >= 0.80,
                }
                for row in pairwise
            ],
            "score",
        )
        for mode in ("pair", "donor", "cluster"):
            low, high, mean = bootstrap_auc(
                pairwise,
                score_field,
                mode,
                cluster_of,
                replicates=300,
            )
            bootstrap_rows.append(
                {
                    "score": score_field,
                    "bootstrap_unit": mode,
                    "auc": f"{overall:.4f}",
                    "bootstrap_ci_low": f"{low:.4f}" if low is not None else "",
                    "bootstrap_ci_high": f"{high:.4f}" if high is not None else "",
                    "bootstrap_mean": f"{mean:.4f}" if mean is not None else "",
                    "ci_width": f"{(high - low):.4f}" if low is not None else "",
                }
            )
    write_tsv(args.out / "pairwise-bootstrap-comparison.tsv", bootstrap_rows)
    print("validation rows", len(validation_rows))
    print("leakage", leakage)
    print("bootstrap", bootstrap_rows)


if __name__ == "__main__":
    main()
