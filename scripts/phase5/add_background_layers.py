#!/usr/bin/env python3
"""Add normalized scan/testability/promoter layers to a ReguTrace background DB.

This is additive: existing tf_catalogue/tf_motif/motif_target/pairs/hits tables
are not modified. The 8.2M-row site table remains an external optional
attachment; motif_target is its compact aggregate.
"""
import argparse
import csv
import json
import os
import sqlite3
import time


def quote(name):
    return '"' + name.replace('"', '""') + '"'


def load_table(conn, path, table, columns):
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in columns if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path}: missing columns {missing}")
        conn.execute(f"DROP TABLE IF EXISTS {table}")
        conn.execute(f"CREATE TABLE {table} ({', '.join(quote(c) + ' TEXT' for c in columns)})")
        placeholders = ",".join("?" for _ in columns)
        batch = []
        for row in reader:
            batch.append(tuple(row.get(c, "") for c in columns))
            if len(batch) >= 5000:
                conn.executemany(f"INSERT INTO {table} VALUES ({placeholders})", batch)
                batch = []
        if batch:
            conn.executemany(f"INSERT INTO {table} VALUES ({placeholders})", batch)
    return conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--testability", required=True)
    ap.add_argument("--promoter-summary", required=True)
    ap.add_argument("--scan-manifest", required=True)
    ap.add_argument("--meta-out", default="")
    args = ap.parse_args()
    start = time.time()
    conn = sqlite3.connect(args.db, timeout=60)
    conn.execute("PRAGMA user_version = 2")

    n_test = load_table(conn, args.testability, "motif_testability", [
        "motif_id", "pwm_hash", "matrix_length", "information_content", "max_score",
        "min_pvalue", "requested_pvalue", "attainable", "status", "mapped",
        "motif_id_hit_count", "tf_count", "background_mode",
    ])
    n_prom = load_table(conn, args.promoter_summary, "promoter", [
        "gene_id", "chrom", "promoter_start", "promoter_end", "promoter_strand",
        "promoter_length", "scanned_pwms", "testable_pwms", "pwms_with_hits",
        "sites", "best_pvalue", "best_score", "best_pwm_hash", "best_motif_ids",
    ])
    conn.execute("CREATE INDEX IF NOT EXISTS motif_testability_id ON motif_testability(motif_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS promoter_gene ON promoter(gene_id)")

    manifest = json.load(open(args.scan_manifest, encoding="utf-8"))
    conn.execute("DROP TABLE IF EXISTS scan_config")
    conn.execute("CREATE TABLE scan_config (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    motif_input = (manifest.get("inputs") or {}).get("motifs") or ""
    scan_jobs = manifest.get("scan_jobs") or {}
    flat = {
        "scan_manifest_path": os.path.basename(args.scan_manifest),
        "background_mode": (manifest.get("background") or {}).get("mode", ""),
        "background_values_acgt": json.dumps((manifest.get("background") or {}).get("values_acgt", [])),
        "motif_library": os.path.basename(motif_input),
        "exact_pwm_reuse": (manifest.get("motifs") or {}).get("exact_pwm_reuse"),
        "pwm_groups": (manifest.get("motifs") or {}).get("pwm_groups"),
        "mapped_motif_ids": (manifest.get("motifs") or {}).get("mapped_motif_ids"),
        "promoter_count": (manifest.get("promoters") or {}).get("count"),
        "scan_seconds": (manifest.get("runtime") or {}).get("seconds"),
        "cache_hits": scan_jobs.get("cache_hits"),
        "cache_misses": scan_jobs.get("cache_misses"),
        "strands": (manifest.get("scan_algorithm") or {}).get("strands"),
        "sites_only": True,
    }
    conn.executemany("INSERT OR REPLACE INTO scan_config VALUES (?,?)",
                     [(k, json.dumps(v) if not isinstance(v, str) else v) for k, v in flat.items()])

    meta_rows = dict(conn.execute("SELECT key,value FROM meta"))
    meta = {}
    for key, value in meta_rows.items():
        try:
            meta[key] = json.loads(value)
        except Exception:
            meta[key] = value
    meta["schema_version"] = "regutrace-regulatory-background-db-v2"
    meta["tables"] = ["tf_catalogue", "tf_motif", "motif_target", "motif_registry",
                      "motif_testability", "promoter", "scan_config", "pairs", "hits"]
    counts = {name: conn.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
              for name in meta["tables"]}
    meta["counts"] = counts
    meta["normalized_layers"] = {
        "scan_config": "scan parameters and background",
        "motif_testability": "attainable p/status per motif",
        "promoter": "per-promoter scanned/testable/site summary",
        "site_table": "external optional attachment (not materialized in this DB)",
    }
    conn.executemany("INSERT OR REPLACE INTO meta VALUES (?,?)",
                     [(k, json.dumps(v) if not isinstance(v, str) else v) for k, v in meta.items()])
    conn.commit()
    conn.execute("VACUUM")
    conn.close()
    meta["database_bytes"] = os.path.getsize(args.db)
    meta["layer_seconds"] = round(time.time() - start, 2)
    if args.meta_out:
        json.dump(meta, open(args.meta_out, "w"), indent=2)
    print(json.dumps({"testability": n_test, "promoter": n_prom,
                      "database_bytes": meta["database_bytes"],
                      "seconds": meta["layer_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
