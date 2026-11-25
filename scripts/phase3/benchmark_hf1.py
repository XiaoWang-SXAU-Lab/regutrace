#!/usr/bin/env python3
"""Summarize the five Phase 3 HF1 benchmark runs.

This script reads only completed run directories. It does not rescan or modify
evidence. Runtime values come from the scan manifest; wall-clock values are
recorded explicitly because the legacy pipeline has no runtime field.
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import statistics
from pathlib import Path


def median(values):
    return statistics.median(values) if values else None


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_lines(path):
    with open(path, "rb") as handle:
        return sum(1 for _ in handle) - 1


def db_metrics(db_path):
    con = sqlite3.connect(str(db_path))
    metrics = {}
    metrics["evidence_rows_db"] = con.execute("select count(*) from hits").fetchone()[0]
    metrics["tf_target_pairs"] = con.execute("select count(*) from pairs").fetchone()[0]
    metrics["mapped_tfs_with_targets"] = con.execute(
        "select count(distinct tf_id) from pairs"
    ).fetchone()[0]
    metrics["target_union"] = con.execute(
        "select count(distinct target_id) from pairs"
    ).fetchone()[0]
    metrics["motif_target_pairs"] = con.execute(
        "select count(*) from (select distinct motif_id,target_id from hits)"
    ).fetchone()[0]
    metrics["evidence_motifs"] = con.execute(
        "select count(distinct motif_id) from hits"
    ).fetchone()[0]
    per_tf = [
        row[0]
        for row in con.execute("select count(*) from pairs group by tf_id").fetchall()
    ]
    metrics["targets_per_tf_min"] = min(per_tf) if per_tf else None
    metrics["targets_per_tf_median"] = median(per_tf)
    metrics["targets_per_tf_mean"] = (
        sum(per_tf) / len(per_tf) if per_tf else None
    )
    metrics["targets_per_tf_max"] = max(per_tf) if per_tf else None
    metrics["tf_with_zero_targets"] = metrics["mapped_tfs_with_targets"]
    # Counted below from the tf-map; stored here so the caller can fill it.
    con.close()
    return metrics


def evidence_set_diff(left_db, right_db, left_filter_tfs=None):
    con = sqlite3.connect(":memory:")
    left_sql_path = str(left_db).replace("'", "''")
    right_sql_path = str(right_db).replace("'", "''")
    con.execute(f"attach database '{left_sql_path}' as l")
    con.execute(f"attach database '{right_sql_path}' as r")
    columns = (
        "tf_id,tf_symbol,target_id,target_symbol,evidence_type,evidence_level,"
        "method,source,pmid,condition,promoter_start,promoter_end,motif_id,"
        "hit_start,hit_end,strand,score,p_value,notes"
    )
    left_where = ""
    right_where = ""
    if left_filter_tfs:
        placeholders = ",".join("?" for _ in left_filter_tfs)
        left_where = f" where tf_id in ({placeholders})"
        right_where = f" where tf_id in ({placeholders})"
        params = tuple(sorted(left_filter_tfs))
        left_sql = (
            f"select {columns} from l.hits{left_where} "
            f"except select {columns} from r.hits"
        )
        right_sql = (
            f"select {columns} from r.hits{right_where} "
            f"except select {columns} from l.hits"
        )
        left_only = con.execute(left_sql, params).fetchall()
        right_only = con.execute(right_sql, params).fetchall()
    else:
        left_only = con.execute(
            f"select {columns} from l.hits except select {columns} from r.hits"
        ).fetchall()
        right_only = con.execute(
            f"select {columns} from r.hits except select {columns} from l.hits"
        ).fetchall()
    con.close()
    return len(left_only), len(right_only)


def site_metrics(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        sites = 0
        tf_count_sum = 0
        missing_coords = 0
        promoter_sites = {}
        for row in reader:
            sites += 1
            tf_count_sum += int(row["tf_count"])
            if not row["genomic_hit_start"]:
                missing_coords += 1
            gene = row["gene_id"]
            promoter_sites[gene] = promoter_sites.get(gene, 0) + 1
    values = list(promoter_sites.values())
    return {
        "unique_sites": sites,
        "site_tf_count_sum": tf_count_sum,
        "sites_missing_genomic_coords": missing_coords,
        "promoters_with_sites": len(values),
        "sites_per_promoter_mean": sum(values) / len(values) if values else 0.0,
        "sites_per_promoter_median": median(values) if values else 0,
        "sites_per_promoter_max": max(values) if values else 0,
    }


def testability_metrics(path):
    rows = read_tsv(path)
    counts = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    mapped = [row for row in rows if row["mapped"] == "true"]
    by_length = {}
    for row in mapped:
        length = int(row["matrix_length"])
        bucket = by_length.setdefault(
            length, {"motifs": 0, "attainable": 0, "with_hits": 0}
        )
        bucket["motifs"] += 1
        bucket["attainable"] += int(row["attainable"] == "true")
        bucket["with_hits"] += int(int(row["motif_id_hit_count"]) > 0)
    return counts, by_length


def tf_map_metrics(path):
    rows = read_tsv(path)
    tfs = {row["tf_id"] for row in rows}
    motifs = {row["motif_id"] for row in rows}
    return {
        "tf_motif_rows": len(rows),
        "mapped_tfs": len(tfs),
        "mapped_motif_ids": len(motifs),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run",
        action="append",
        required=True,
        help="label=directory; repeat for each benchmark run",
    )
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--wall-seconds",
        action="append",
        default=[],
        help="label=seconds; wall-clock time including pipeline stages",
    )
    parser.add_argument("--run1-run2-equality", nargs=2, metavar=("RUN1", "RUN2"))
    parser.add_argument(
        "--run5-strict-subset",
        nargs=3,
        metavar=("RUN5", "RUN3", "RUN3_TF_MAP"),
        help="compare Run 5 rows for Run 3 TFs with the full Run 3 evidence",
    )
    args = parser.parse_args()

    runs = {}
    for item in args.run:
        label, path = item.split("=", 1)
        runs[label] = Path(path)
    walls = {}
    for item in args.wall_seconds:
        label, value = item.split("=", 1)
        walls[label] = float(value)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    per_tf_rows = []
    length_rows = []
    for label, directory in runs.items():
        manifest = {}
        manifest_path = directory / "scan-manifest.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        pipeline = {}
        pipeline_path = directory / "pipeline-report.json"
        if pipeline_path.exists():
            pipeline = json.loads(pipeline_path.read_text(encoding="utf-8"))
        tf_map_path = directory / "tf_map.tsv"
        if not tf_map_path.exists():
            inputs = manifest.get("inputs", {})
            if inputs.get("tf_map"):
                tf_map_path = Path(inputs["tf_map"])
        map_metrics = tf_map_metrics(tf_map_path)
        db_path = directory / "evidence.db"
        metrics = db_metrics(db_path)
        metrics["mapped_tfs_with_targets"] = metrics["mapped_tfs_with_targets"]
        metrics["tfs_without_targets"] = (
            map_metrics["mapped_tfs"] - metrics["mapped_tfs_with_targets"]
        )
        evidence_path = directory / "evidence.tsv"
        metrics["evidence_rows_tsv"] = count_lines(evidence_path)
        sites_path = directory / "motif-site-hits.tsv"
        if not sites_path.exists():
            sites_path = directory / "sites.tsv"
        if sites_path.exists():
            metrics.update(site_metrics(sites_path))
            metrics["sites_file_bytes"] = sites_path.stat().st_size
        else:
            metrics.update(
                {
                    "unique_sites": None,
                    "site_tf_count_sum": None,
                    "sites_missing_genomic_coords": None,
                    "promoters_with_sites": None,
                    "sites_per_promoter_mean": None,
                    "sites_per_promoter_median": None,
                    "sites_per_promoter_max": None,
                    "sites_file_bytes": None,
                }
            )
        testability_path = directory / "motif-testability.tsv"
        if not testability_path.exists():
            testability_path = directory / "testability.tsv"
        if testability_path.exists():
            counts, by_length = testability_metrics(testability_path)
            metrics["testability_counts"] = counts
            metrics["testable_motifs"] = counts.get("significant_hits", 0) + counts.get(
                "no_significant_hits", 0
            )
            metrics["unattainable_motifs"] = counts.get(
                "threshold_unattainable", 0
            )
            metrics["not_scanned_motifs"] = counts.get("not_scanned", 0)
            for length, values in sorted(by_length.items()):
                length_rows.append(
                    {
                        "run": label,
                        "matrix_length": length,
                        **values,
                    }
                )
        else:
            scanned_unreachable = (
                pipeline.get("motifs", {})
                .get("threshold", {})
                .get("scanned_unreachable", {})
                .get("count")
            )
            metrics["testability_counts"] = {}
            metrics["testable_motifs"] = None
            metrics["unattainable_motifs"] = scanned_unreachable
            metrics["not_scanned_motifs"] = None

        scan = manifest.get("scan_jobs", {})
        sites = manifest.get("sites", {})
        motifs = manifest.get("motifs", {})
        background = manifest.get("background", {})
        phase3 = pipeline.get("phase3", {})
        metrics.update(
            {
                "run": label,
                "directory": str(directory),
                "pack": pipeline.get("reference_pack", {}).get(
                    "version", manifest.get("reference_pack", "")
                ),
                "background_mode": background.get(
                    "mode", phase3.get("background_mode", "uniform")
                ),
                "background_acgt": ",".join(
                    f"{value:.6f}" for value in background.get("values_acgt", [])
                ),
                "tf_motif_rows": map_metrics["tf_motif_rows"],
                "mapped_tfs": map_metrics["mapped_tfs"],
                "mapped_motif_ids": map_metrics["mapped_motif_ids"],
                "pwm_groups": motifs.get("pwm_groups", phase3.get("pwm_groups")),
                "unique_pwms": motifs.get("unique_pwms", phase3.get("unique_pwms")),
                "duplicate_pwm_groups": motifs.get(
                    "duplicate_pwm_groups", phase3.get("duplicate_pwm_groups")
                ),
                "exact_pwm_reuse": motifs.get("exact_pwm_reuse", False),
                "scan_runtime_seconds": manifest.get("runtime", {}).get("seconds"),
                "wall_seconds": walls.get(label),
                "cache_hits": scan.get("cache_hits", phase3.get("cache_hits")),
                "cache_misses": scan.get("cache_misses", phase3.get("cache_misses")),
                "executed_group_promoter_pairs": scan.get(
                    "executed_group_promoter_pairs"
                ),
                "group_promoter_pairs": scan.get("group_promoter_pairs"),
                "skipped_unattainable_pairs": scan.get(
                    "skipped_unattainable_group_promoter_pairs"
                ),
                "unique_site_hits": sites.get(
                    "unique_pwm_sites", phase3.get("unique_site_hits")
                ),
                "motif_id_site_hits": sites.get(
                    "motif_id_site_hits", phase3.get("motif_id_hits")
                ),
                "scan_evidence_rows": sites.get(
                    "evidence_rows", phase3.get("evidence_rows")
                ),
                "evidence_tsv_bytes": evidence_path.stat().st_size,
                "evidence_db_bytes": db_path.stat().st_size,
                "input_identified_tfs": pipeline.get("tf_catalogue", {}).get(
                    "identified"
                ),
                "promoters": manifest.get("promoters", {}).get("count")
                or pipeline.get("scan", {}).get("promoters"),
            }
        )
        rows.append(metrics)

        con = sqlite3.connect(str(db_path))
        for tf_id, target_count in con.execute(
            "select tf_id,count(*) from pairs group by tf_id order by tf_id"
        ):
            per_tf_rows.append(
                {"run": label, "tf_id": tf_id, "targets": target_count}
            )
        con.close()

    # Column order: stable, machine-readable.
    preferred = [
        "run",
        "pack",
        "background_mode",
        "background_acgt",
        "input_identified_tfs",
        "tf_motif_rows",
        "mapped_tfs",
        "mapped_motif_ids",
        "pwm_groups",
        "unique_pwms",
        "duplicate_pwm_groups",
        "exact_pwm_reuse",
        "promoters",
        "group_promoter_pairs",
        "executed_group_promoter_pairs",
        "skipped_unattainable_pairs",
        "cache_hits",
        "cache_misses",
        "unique_site_hits",
        "motif_id_site_hits",
        "scan_evidence_rows",
        "evidence_rows_db",
        "evidence_rows_tsv",
        "motif_target_pairs",
        "tf_target_pairs",
        "target_union",
        "mapped_tfs_with_targets",
        "tfs_without_targets",
        "targets_per_tf_min",
        "targets_per_tf_median",
        "targets_per_tf_mean",
        "targets_per_tf_max",
        "promoters_with_sites",
        "sites_per_promoter_mean",
        "sites_per_promoter_median",
        "sites_per_promoter_max",
        "sites_missing_genomic_coords",
        "testable_motifs",
        "unattainable_motifs",
        "not_scanned_motifs",
        "scan_runtime_seconds",
        "wall_seconds",
        "sites_file_bytes",
        "evidence_tsv_bytes",
        "evidence_db_bytes",
    ]
    all_keys = []
    for row in rows:
        for key in row:
            if key not in all_keys:
                all_keys.append(key)
    columns = preferred + [key for key in all_keys if key not in preferred]
    with open(out / "hf1-benchmark.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        json.dumps(row.get(key), sort_keys=True)
                        if isinstance(row.get(key), (dict, list))
                        else row.get(key)
                    )
                    for key in columns
                }
            )
    with open(out / "hf1-benchmark.json", "w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with open(out / "hf1-targets-per-tf.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["run", "tf_id", "targets"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(per_tf_rows)
    with open(
        out / "hf1-testability-by-length.tsv", "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["run", "matrix_length", "motifs", "attainable", "with_hits"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(length_rows)

    if args.run1_run2_equality:
        left, right = args.run1_run2_equality
        left_db = runs[left] / "evidence.db"
        right_db = runs[right] / "evidence.db"
        left_only, right_only = evidence_set_diff(left_db, right_db)
        equality = {
            "comparison": f"{left}_vs_{right}",
            "left_only_rows": left_only,
            "right_only_rows": right_only,
            "identical_row_sets": left_only == 0 and right_only == 0,
        }
        (out / "hf1-run1-run2-equality.json").write_text(
            json.dumps(equality, indent=2) + "\n", encoding="utf-8"
        )

    if args.run5_strict_subset:
        run5, run3, run3_map = args.run5_strict_subset
        strict_tfs = sorted({row["tf_id"] for row in read_tsv(Path(run3_map))})
        left_only, right_only = evidence_set_diff(
            runs[run5] / "evidence.db",
            runs[run3] / "evidence.db",
            left_filter_tfs=strict_tfs,
        )
        equality = {
            "comparison": f"{run5}_strict_subset_vs_{run3}",
            "strict_tfs": len(strict_tfs),
            "left_only_rows": left_only,
            "right_only_rows": right_only,
            "identical_row_sets": left_only == 0 and right_only == 0,
        }
        (out / "hf1-run5-strict-subset-equality.json").write_text(
            json.dumps(equality, indent=2) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    main()
