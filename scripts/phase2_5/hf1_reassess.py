#!/usr/bin/env python3
"""Compare HF1 strict, exploratory DBD and calibrated DBD mapping methods."""
import argparse
import collections
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_tsv, write_tsv  # noqa: E402


def load_map(path):
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def strict_by_tf(rows):
    return {row["tf_id"]: row for row in rows}


def dbd_by_tf(rows):
    grouped = collections.defaultdict(list)
    for row in rows:
        grouped[row["tf_id"]].append(row)
    return grouped


def motif_similarity_lookup(path):
    lookup = {}
    for row in read_tsv(path):
        key = tuple(sorted((row["motif_a"], row["motif_b"])))
        lookup[key] = (float(row["correlation"]), float(row["overlap_fraction"]))
    return lookup


def motif_match(a, b, lookup, threshold=0.80):
    if a == b:
        return True
    correlation, overlap = lookup.get(tuple(sorted((a, b))), (-1.0, 0.0))
    return correlation >= threshold and overlap >= 0.50


def base_id(motif):
    return re.sub(r"\.\d+$", "", motif)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict2024", type=Path, required=True)
    parser.add_argument("--strict2026", type=Path, required=True)
    parser.add_argument("--dbd-uncalibrated", type=Path, required=True)
    parser.add_argument("--dbd-calibrated", type=Path, required=True)
    parser.add_argument("--dbd-diagnostics", type=Path, required=True)
    parser.add_argument("--catalogue", type=Path, required=True)
    parser.add_argument("--motif-similarity", type=Path, required=True)
    parser.add_argument("--reference-map-2026", type=Path, required=True)
    parser.add_argument("--motif-library-2026", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    strict2024 = strict_by_tf(load_map(args.strict2024))
    strict2026 = strict_by_tf(load_map(args.strict2026))
    dbd_uncalibrated = dbd_by_tf(load_map(args.dbd_uncalibrated))
    dbd_calibrated = dbd_by_tf(load_map(args.dbd_calibrated))
    diagnostics = load_map(args.dbd_diagnostics)
    similarity = motif_similarity_lookup(args.motif_similarity)
    catalogue = load_map(args.catalogue)
    catalogue_family = {row["gene_id"]: row["family"] for row in catalogue}
    family_counts = collections.Counter(catalogue_family.values())
    donors2026 = {row["tf_id"] for row in read_tsv(args.reference_map_2026)}
    motifs2026 = {
        line[1:].split()[0]
        for line in open(args.motif_library_2026, encoding="utf-8", errors="replace")
        if line.startswith(">")
    }

    calibrated_tfs = {
        tf
        for tf, rows in dbd_calibrated.items()
        if any(row["mapping_status"] == "dbd_calibrated" for row in rows)
    }
    ambiguous_tfs = {
        tf
        for tf, rows in dbd_calibrated.items()
        if all(row["mapping_status"] != "dbd_calibrated" for row in rows)
    }
    uncalibrated_tfs = set(dbd_uncalibrated)
    combined = dict(strict2026)
    for tf in calibrated_tfs:
        if tf not in combined:
            rows = [row for row in dbd_calibrated[tf] if row["mapping_status"] == "dbd_calibrated"]
            if rows:
                combined[tf] = rows[0]

    def recovery(rows_by_tf, cutoff):
        exact_donor = 0
        exact = 0
        similar = 0
        covered = 0
        for tf, rows in rows_by_tf.items():
            if tf not in strict2026:
                continue
            covered += 1
            strict = strict2026[tf]
            donors = {row["reference_tf_id"] for row in rows}
            motifs = {row["motif_id"] for row in rows}
            if strict["reference_tf_id"] in donors:
                exact_donor += 1
            if strict["motif_id"] in motifs:
                exact += 1
            if any(motif_match(strict["motif_id"], motif, similarity, cutoff) for motif in motifs):
                similar += 1
        return covered, exact_donor, exact, similar

    uncal_covered, uncal_donor, uncal_exact, uncal_similar = recovery(dbd_uncalibrated, 0.80)
    cal_covered, cal_donor, cal_exact, cal_similar = recovery(
        {tf: rows for tf, rows in dbd_calibrated.items() if tf in calibrated_tfs}, 0.80
    )
    combined_covered, combined_donor, combined_exact, combined_similar = recovery(
        {tf: [row] for tf, row in combined.items()}, 0.80
    )

    def unique_motifs(rows_by_tf):
        return len({row["motif_id"] for rows in rows_by_tf.values() for row in rows})

    methods = [
        ("strict_2024", len(strict2024), len(strict2024), 0, 0, unique_motifs({tf: [row] for tf, row in strict2024.items()}), "", "", "", ""),
        ("strict_2026", len(strict2026), len(strict2026), 0, 0, unique_motifs({tf: [row] for tf, row in strict2026.items()}), "", "", "", ""),
        ("dbd_uncalibrated", sum(len(rows) for rows in dbd_uncalibrated.values()), len(uncalibrated_tfs), 0, 0, unique_motifs(dbd_uncalibrated), uncal_covered, uncal_donor, uncal_exact, uncal_similar),
        ("dbd_calibrated", sum(len(rows) for rows in dbd_calibrated.values()), len(calibrated_tfs), len(calibrated_tfs), len(ambiguous_tfs), unique_motifs(dbd_calibrated), cal_covered, cal_donor, cal_exact, cal_similar),
        ("strict_2026_plus_calibrated", len(combined), len(combined), len(calibrated_tfs), len(ambiguous_tfs), unique_motifs({tf: [row] for tf, row in combined.items()}), combined_covered, combined_donor, combined_exact, combined_similar),
    ]
    method_rows = []
    for method, rows, tfs, calibrated, ambiguous, motifs, covered, donor, exact, similar in methods:
        method_rows.append(
            {
                "method": method,
                "mapping_rows": rows,
                "mapped_tfs": tfs,
                "calibrated_tfs": calibrated,
                "ambiguous_tfs": ambiguous,
                "unique_motifs": motifs,
                "strict2026_overlap_tfs": covered,
                "strict2026_exact_donor": donor,
                "strict2026_exact_motif": exact,
                "strict2026_motif_similar_080": similar,
            }
        )
    combined_rows = []
    for tf, row in sorted(strict2026.items()):
        combined_rows.append(
            {
                "tf_id": tf,
                "family": row["family"],
                "motif_id": row["motif_id"],
                "reference_tf_id": row["reference_tf_id"],
                "mapping_status": "strict_rbh",
                "mapping_source": "strict_2026",
            }
        )
    for tf in sorted(calibrated_tfs):
        if tf in strict2026:
            continue
        rows = [row for row in dbd_calibrated[tf] if row["mapping_status"] == "dbd_calibrated"]
        if not rows:
            continue
        row = rows[0]
        combined_rows.append(
            {
                "tf_id": tf,
                "family": row["family"],
                "motif_id": row["motif_id"],
                "reference_tf_id": row["reference_tf_id"],
                "mapping_status": "dbd_calibrated",
                "mapping_source": "calibrated_dbd",
            }
        )
    write_tsv(args.out / "hf1-strict-plus-calibrated.tsv", combined_rows)

    write_tsv(args.out / "hf1-method-summary.tsv", method_rows)

    families = sorted(family_counts)
    family_rows = []
    for family in families:
        strict_family = {tf for tf, row in strict2026.items() if row["family"] == family}
        family_rows.append(
            {
                "family": family,
                "catalogue_tfs": family_counts[family],
                "strict_2024_tfs": sum(1 for row in strict2024.values() if row["family"] == family),
                "strict_2026_tfs": len(strict_family),
                "dbd_uncalibrated_tfs": sum(1 for tf in uncalibrated_tfs if catalogue_family.get(tf, "") == family),
                "dbd_calibrated_tfs": sum(1 for tf in calibrated_tfs if catalogue_family.get(tf, "") == family),
                "dbd_ambiguous_tfs": sum(1 for tf in ambiguous_tfs if catalogue_family.get(tf, "") == family),
                "combined_tfs": sum(1 for tf in combined if catalogue_family.get(tf, "") == family),
            }
        )
    write_tsv(args.out / "hf1-family-coverage.tsv", family_rows)

    changes = []
    for tf in sorted(set(strict2024) | set(strict2026)):
        old = strict2024.get(tf)
        new = strict2026.get(tf)
        if old and not new:
            motif = old["motif_id"]
            donor = old["reference_tf_id"]
            if donor not in donors2026:
                reason = "donor_absent_from_2026_pack"
            elif motif not in motifs2026:
                reason = "motif_absent_from_2026_library"
            else:
                reason = "competition_or_threshold_change"
            changes.append(
                {
                    "tf_id": tf,
                    "family": old["family"],
                    "change": "lost",
                    "old_donor": donor,
                    "old_motif": motif,
                    "new_donor": "",
                    "new_motif": "",
                    "donor_present_2026": "yes" if donor in donors2026 else "no",
                    "motif_present_2026": "yes" if motif in motifs2026 else "no",
                    "detail": reason,
                }
            )
        elif new and not old:
            changes.append(
                {
                    "tf_id": tf,
                    "family": new["family"],
                    "change": "new",
                    "old_donor": "",
                    "old_motif": "",
                    "new_donor": new["reference_tf_id"],
                    "new_motif": new["motif_id"],
                    "donor_present_2026": "yes",
                    "motif_present_2026": "yes",
                    "detail": "",
                }
            )
        elif old and new and (
            old["reference_tf_id"] != new["reference_tf_id"] or old["motif_id"] != new["motif_id"]
        ):
            if base_id(old["motif_id"]) == base_id(new["motif_id"]):
                detail = "motif_version_update"
            elif old["reference_tf_id"] == new["reference_tf_id"]:
                detail = "motif_representative_change"
            else:
                detail = "donor_change"
            changes.append(
                {
                    "tf_id": tf,
                    "family": new["family"],
                    "change": "changed",
                    "old_donor": old["reference_tf_id"],
                    "old_motif": old["motif_id"],
                    "new_donor": new["reference_tf_id"],
                    "new_motif": new["motif_id"],
                    "donor_present_2026": "yes" if new["reference_tf_id"] in donors2026 else "no",
                    "motif_present_2026": "yes" if new["motif_id"] in motifs2026 else "no",
                    "detail": detail,
                }
            )
    write_tsv(args.out / "hf1-strict-2024-vs-2026-changes.tsv", changes)

    diagnostic_counts = collections.Counter(row["primary_reason"] for row in diagnostics)
    diagnostic_family = collections.Counter()
    for row in diagnostics:
        if row["primary_reason"] in {"DBD_CALIBRATED", "DBD_AMBIGUOUS", "NOT_RECIPROCAL_DBD_BEST", "BELOW_DBD_THRESHOLD"}:
            diagnostic_family[(row["family"], row["primary_reason"])] += 1
    ambiguity_rows = []
    for tf, rows in sorted(dbd_calibrated.items()):
        ambiguity_rows.append(
            {
                "tf_id": tf,
                "family": rows[0]["family"],
                "mapping_status": "dbd_calibrated" if tf in calibrated_tfs else "dbd_ambiguous",
                "row_count": len(rows),
                "specificity_groups": "|".join(sorted({row["specificity_groups"] for row in rows})),
                "reference_donors": "|".join(sorted({row["reference_tf_id"] for row in rows})),
                "motifs": "|".join(sorted({row["motif_id"] for row in rows})),
            }
        )
    write_tsv(args.out / "hf1-dbd-ambiguity.tsv", ambiguity_rows)

    summary = {
        "methods": method_rows,
        "diagnostic_primary_reasons": dict(diagnostic_counts),
        "diagnostic_family_reasons": [
            {"family": family, "reason": reason, "count": count}
            for (family, reason), count in sorted(diagnostic_family.items())
        ],
        "strict_2024_vs_2026_changes": changes,
        "notes": [
            "strict_2026_plus_calibrated uses strict rows first and adds calibrated DBD TFs only when the TF is absent from strict.",
            "DBD ambiguities are not promoted to calibrated calls.",
            "motif similarity is TOMTOM-style Pearson correlation with at least 50% overlap; it is not biological equivalence.",
        ],
    }
    (args.out / "hf1-reassess-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print("methods:")
    for row in method_rows:
        print(row)
    print("changes", len(changes), "ambiguity rows", len(ambiguity_rows))


if __name__ == "__main__":
    main()
