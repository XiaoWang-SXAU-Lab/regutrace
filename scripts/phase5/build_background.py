#!/usr/bin/env python3
"""Build a standalone, importable ReguTrace regulatory-background SQLite.

The database contains the complete TF catalogue, TF-motif mapping, aggregated
motif-target associations, registry metadata, an empty canonical hits table,
and a policy-filtered TF-target pairs table. It is independent of the running
service and the original job ID.
"""
import argparse
import csv
import hashlib
import json
import os
import sqlite3
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
from aggregate_audit import (dedup_exact, group_rows, read_coord_map,
                             read_pwm_hash, read_specificity_groups)


HITS_SQL = """
CREATE TABLE IF NOT EXISTS hits (
 tf_id TEXT NOT NULL, tf_symbol TEXT NOT NULL DEFAULT '',
 target_id TEXT NOT NULL, target_symbol TEXT NOT NULL DEFAULT '',
 evidence_type TEXT NOT NULL DEFAULT '', evidence_level TEXT NOT NULL DEFAULT '',
 method TEXT NOT NULL DEFAULT '', source TEXT NOT NULL DEFAULT '',
 pmid TEXT NOT NULL DEFAULT '', condition TEXT NOT NULL DEFAULT '',
 promoter_start INTEGER, promoter_end INTEGER,
 motif_id TEXT NOT NULL DEFAULT '', hit_start INTEGER, hit_end INTEGER,
 strand TEXT NOT NULL DEFAULT '', score REAL, p_value REAL,
 notes TEXT NOT NULL DEFAULT ''
);
"""

PAIRS_SQL = """
CREATE TABLE IF NOT EXISTS pairs (
 tf_id TEXT NOT NULL, tf_symbol TEXT NOT NULL, target_id TEXT NOT NULL,
 target_symbol TEXT NOT NULL, n_sites INTEGER NOT NULL, n_motifs INTEGER NOT NULL,
 best_p REAL, best_score REAL, motifs TEXT, evidence_levels TEXT
);
"""


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def best_motifs_by_tf(coord, pwm_hash, specificity, threshold=0.8):
    selected = {}
    for tf, rows in coord.items():
        chosen = []
        for group in group_rows(rows, threshold):
            unique = dedup_exact(group, pwm_hash)
            if not unique:
                continue
            unique.sort(key=lambda r: (r["_rank"], r["motif_id"]), reverse=True)
            chosen.append(unique[0]["motif_id"])
        selected[tf] = sorted(set(chosen))
    return selected


def create_registry(conn, motifs_path, pwm_hash, specificity):
    conn.execute("CREATE TABLE motif_registry (motif_id TEXT PRIMARY KEY, motif_name TEXT, family TEXT, pwm_hash TEXT, matrix_length INTEGER, specificity_group_080 TEXT, source TEXT)")
    if not motifs_path or not os.path.exists(motifs_path):
        return
    motif_id = None
    name = None
    length = 0
    rows = []
    with open(motifs_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if motif_id:
                    rows.append((motif_id, name, "", pwm_hash.get(motif_id, ""), length, specificity.get(motif_id, ""), "JASPAR2024"))
                parts = line[1:].strip().split(None, 1)
                motif_id = parts[0] if parts else ""
                name = parts[1] if len(parts) > 1 else ""
                length = 0
            elif motif_id and line[:1] in "ACGT":
                length += len(line.split("[")[1].split("]")[0].split()) if "[" in line else 0
        if motif_id:
            rows.append((motif_id, name, "", pwm_hash.get(motif_id, ""), length, specificity.get(motif_id, ""), "JASPAR2024"))
    conn.executemany("INSERT OR REPLACE INTO motif_registry VALUES (?,?,?,?,?,?,?)", rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coord-map", required=True)
    ap.add_argument("--base-map", required=True)
    ap.add_argument("--tf-genes", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--testability", required=True)
    ap.add_argument("--motifs", default="")
    ap.add_argument("--pwm-hash", default="")
    ap.add_argument("--specificity-groups", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--policy", default="best@0.8")
    ap.add_argument("--manifest", default="")
    args = ap.parse_args()
    start = time.time()
    if os.path.exists(args.out):
        os.remove(args.out)
    conn = sqlite3.connect(args.out)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    conn.executescript(HITS_SQL)
    conn.executescript(PAIRS_SQL)
    coord = read_coord_map(args.coord_map)
    base = read_coord_map(args.base_map)
    pwm_hash = read_pwm_hash(args.pwm_hash)
    specificity = read_specificity_groups(args.specificity_groups)

    # TF catalogue.
    conn.execute("CREATE TABLE tf_catalogue (tf_id TEXT PRIMARY KEY, tf_symbol TEXT, family TEXT, identification_status TEXT, source TEXT, motif_count INTEGER, target_count INTEGER)")
    with open(args.tf_genes, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = []
        for row in reader:
            gid = row.get("gene_id", "")
            rows.append((gid, row.get("tf_symbol", ""), row.get("family", ""),
                         row.get("identification_status", ""), row.get("source", ""),
                         len({r["motif_id"] for r in base.get(gid, [])}), 0))
    conn.executemany("INSERT OR REPLACE INTO tf_catalogue VALUES (?,?,?,?,?,?,?)", rows)

    # Full TF-motif mapping (strict/calibrated/exploratory), preserving
    # coordinates when the coordinate-aware rerun provides them.
    with open(args.base_map, encoding="utf-8") as f:
        base_rows = list(csv.DictReader(f, delimiter="\t"))
    coord_by_key = {}
    for tf, rows in coord.items():
        for r in rows:
            coord_by_key[(tf, r["motif_id"], r.get("reference_tf_id", ""))] = r
    fields = list(base_rows[0].keys()) if base_rows else []
    extra = ["query_full_dbd_start", "query_full_dbd_end", "donor_full_dbd_start",
             "donor_full_dbd_end", "query_dbd_start", "query_dbd_end",
             "donor_dbd_start", "donor_dbd_end", "dbd_bitscore", "dbd_evalue",
             "alignment_length"]
    all_fields = fields + [x for x in extra if x not in fields]
    cols = ", ".join(f'"{name}" TEXT' for name in all_fields)
    conn.execute(f"CREATE TABLE tf_motif ({cols})")
    placeholders = ",".join("?" for _ in all_fields)
    batch = []
    for row in base_rows:
        c = coord_by_key.get((row.get("tf_id"), row.get("motif_id"), row.get("reference_tf_id")), {})
        merged = dict(row)
        for x in extra:
            merged[x] = c.get(x, "")
        batch.append(tuple(merged.get(name, "") for name in all_fields))
        if len(batch) >= 5000:
            conn.executemany(f"INSERT INTO tf_motif VALUES ({placeholders})", batch)
            batch = []
    if batch:
        conn.executemany(f"INSERT INTO tf_motif VALUES ({placeholders})", batch)

    # Registry metadata.
    create_registry(conn, args.motifs, pwm_hash, specificity)

    # Motif-target aggregation from the site table.
    conn.execute("CREATE TEMP TABLE raw_motif_target (motif_id TEXT, target_id TEXT, p_value REAL, score REAL)")
    wanted = ({r["motif_id"] for rows in base.values() for r in rows} |
              {r["motif_id"] for rows in coord.values() for r in rows})
    with open(args.sites, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        batch = []
        for line in f:
            fields_row = line.rstrip("\n").split("\t")
            if len(fields_row) <= idx["p_value"]:
                continue
            target = fields_row[idx["gene_id"]]
            try:
                p = float(fields_row[idx["p_value"]])
                score = float(fields_row[idx["score"]]) if "score" in idx and fields_row[idx["score"]] else 0.0
            except ValueError:
                continue
            for motif in fields_row[idx["motif_ids"]].split(";"):
                if motif and motif in wanted:
                    batch.append((motif, target, p, score))
            if len(batch) >= 10000:
                conn.executemany("INSERT INTO raw_motif_target VALUES (?,?,?,?)", batch)
                batch = []
        if batch:
            conn.executemany("INSERT INTO raw_motif_target VALUES (?,?,?,?)", batch)
    conn.execute("CREATE TABLE motif_target AS SELECT motif_id, target_id, COUNT(*) n_sites, MIN(p_value) best_p, MAX(score) best_score FROM raw_motif_target GROUP BY motif_id, target_id")
    conn.execute("DROP TABLE raw_motif_target")
    conn.execute("CREATE INDEX motif_target_motif ON motif_target(motif_id)")
    conn.execute("CREATE INDEX motif_target_target ON motif_target(target_id)")

    # Policy-selected TF-target pairs. This is the import contract; it is not
    # the 25M all-candidate cross-product.
    if args.policy != "none":
        selected = {}
        mode = args.policy.split("@", 1)[0].lower()
        for tf, rows in base.items():
            fixed = [r["motif_id"] for r in rows
                     if r.get("mapping_status") in ("strict_rbh", "dbd_calibrated")]
            crows = [r for r in coord.get(tf, [])
                     if r.get("mapping_status") == "dbd_exploratory"]
            chosen = list(fixed)
            if mode == "all":
                chosen.extend(r["motif_id"] for r in rows
                              if r.get("mapping_status") == "dbd_exploratory")
            elif crows:
                for group in group_rows(crows, 0.8):
                    unique = dedup_exact(group, pwm_hash)
                    if not unique:
                        continue
                    unique.sort(key=lambda r: (r["_rank"], r["motif_id"]), reverse=True)
                    if mode == "reps":
                        seen = set()
                        for row in unique:
                            gid = specificity.get(row["motif_id"], row["motif_id"])
                            if gid not in seen:
                                seen.add(gid)
                                chosen.append(row["motif_id"])
                    else:
                        chosen.append(unique[0]["motif_id"])
            else:
                expl = [r for r in rows if r.get("mapping_status") == "dbd_exploratory"]
                if mode == "reps":
                    seen = set()
                    for row in sorted(expl, key=lambda r: (float(r.get("query_dbd_coverage") or 0), r["motif_id"]), reverse=True):
                        gid = specificity.get(row["motif_id"], row["motif_id"])
                        if gid not in seen:
                            seen.add(gid)
                            chosen.append(row["motif_id"])
                elif expl:
                    chosen.append(sorted(
                        expl,
                        key=lambda r: (float(r.get("query_dbd_coverage") or 0), r["motif_id"]),
                        reverse=True,
                    )[0]["motif_id"])
            selected[tf] = sorted(set(chosen))
        conn.execute("CREATE TABLE tf_motif_selected (tf_id TEXT, motif_id TEXT)")
        conn.executemany("INSERT INTO tf_motif_selected VALUES (?,?)",
                         [(tf, m) for tf, motifs in selected.items() for m in motifs])
        conn.execute("""INSERT INTO pairs
            SELECT s.tf_id, COALESCE(c.tf_symbol, s.tf_id), mt.target_id, '',
                   SUM(mt.n_sites), COUNT(DISTINCT mt.motif_id), MIN(mt.best_p),
                   MAX(mt.best_score), GROUP_CONCAT(DISTINCT mt.motif_id), 'D'
            FROM tf_motif_selected s JOIN motif_target mt ON mt.motif_id=s.motif_id
            LEFT JOIN tf_catalogue c ON c.tf_id=s.tf_id
            GROUP BY s.tf_id, mt.target_id""")
    conn.commit()

    # Catalogue target counts under the selected policy. Aggregate once; a
    # correlated per-catalogue-row subquery is prohibitively slow at this size.
    conn.execute("CREATE INDEX IF NOT EXISTS tms_motif ON tf_motif_selected(motif_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS tms_tf ON tf_motif_selected(tf_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS tfm_motif ON tf_motif(motif_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS tfm_tf ON tf_motif(tf_id)")
    conn.execute("CREATE TEMP TABLE tf_counts AS SELECT s.tf_id, COUNT(DISTINCT mt.target_id) n FROM tf_motif_selected s JOIN motif_target mt ON mt.motif_id=s.motif_id GROUP BY s.tf_id")
    conn.execute("UPDATE tf_catalogue SET target_count=COALESCE((SELECT n FROM tf_counts WHERE tf_counts.tf_id=tf_catalogue.tf_id),0)")
    conn.execute("DROP TABLE tf_counts")
    for table, index, column in [
        ("pairs", "pairs_tf", "tf_id"), ("pairs", "pairs_target", "target_id"),
        ("hits", "hits_tf", "tf_id"), ("hits", "hits_target", "target_id"),
        ("hits", "hits_motif", "motif_id"), ("hits", "hits_p", "p_value")]:
        conn.execute(f"CREATE INDEX IF NOT EXISTS {index} ON {table}({column})")

    stats = {name: conn.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
             for name in ("tf_catalogue", "tf_motif", "motif_target", "pairs", "hits")}
    meta = {
        "schema_version": "regutrace-regulatory-background-db-v1",
        "source_job_mode": "site_only_exploratory",
        "policy": args.policy,
        "promoter_upstream_bp": 1000,
        "promoter_downstream_bp": 0,
        "pvalue": 1e-4,
        "evidence_level": "D",
        "coordinate_system": "protein coordinates for query DBD; donor coordinates per reference protein",
        "tables": ["tf_catalogue", "tf_motif", "motif_target", "motif_registry", "pairs", "hits"],
        "counts": stats,
        "inputs": {
            "coord_map": {"path": os.path.basename(args.coord_map), "sha256": sha256(args.coord_map)},
            "sites": {"path": os.path.basename(args.sites), "sha256": sha256(args.sites)},
            "tf_genes": {"path": os.path.basename(args.tf_genes), "sha256": sha256(args.tf_genes)},
            "testability": {"path": os.path.basename(args.testability), "sha256": sha256(args.testability)},
        },
        "limitations": [
            "All PWM motif-target associations are evidence level D; they are not experimental binding.",
            "The pairs table is policy-selected, not the full exploratory union.",
            "Unobserved genes are not negatives in any downstream benchmark.",
        ],
    }
    conn.executemany("INSERT INTO meta VALUES (?,?)", [(k, json.dumps(v) if not isinstance(v, str) else v) for k, v in meta.items()])
    conn.commit()
    conn.execute("VACUUM")
    conn.close()
    meta["build_seconds"] = round(time.time() - start, 2)
    meta["database_bytes"] = os.path.getsize(args.out)
    with open(os.path.join(os.path.dirname(args.out), "background-build.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
