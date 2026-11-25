#!/usr/bin/env python3
"""Rebuild motif_registry from the exact motif library used by a scan.

Fixes the provenance gap where a 2024-only registry did not cover 2026
mapped motif IDs. Does not change tf_motif, motif_target, pairs, or hits.
"""
import argparse
import csv
import json
import os
import sqlite3


def read_jaspar(path):
    rows = []
    motif = None
    name = ""
    length = 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if motif:
                    rows.append((motif, name, length))
                parts = line[1:].strip().split(None, 1)
                motif = parts[0] if parts else ""
                name = parts[1] if len(parts) > 1 else ""
                length = 0
            elif motif and line[:1] in "ACGT":
                if "[" in line and "]" in line:
                    length = len(line.split("[", 1)[1].split("]", 1)[0].split())
        if motif:
            rows.append((motif, name, length))
    return rows


def read_specificity(path):
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["motif_id"]] = (row.get("pwm_hash", ""), row.get("family", ""), row.get("group_080", ""))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--motifs", required=True)
    ap.add_argument("--specificity-groups", default="")
    args = ap.parse_args()
    motifs = read_jaspar(args.motifs)
    spec = read_specificity(args.specificity_groups)
    conn = sqlite3.connect(args.db, timeout=60)
    conn.execute("DROP TABLE IF EXISTS motif_registry")
    conn.execute("CREATE TABLE motif_registry (motif_id TEXT PRIMARY KEY, motif_name TEXT, family TEXT, pwm_hash TEXT, matrix_length INTEGER, specificity_group_080 TEXT, source TEXT)")
    conn.executemany(
        "INSERT INTO motif_registry VALUES (?,?,?,?,?,?,?)",
        [(m, n, spec.get(m, ("", "", ""))[1], spec.get(m, ("", "", ""))[0], length,
          spec.get(m, ("", "", ""))[2], "ReguTrace job motif library") for m, n, length in motifs],
    )
    conn.execute("CREATE INDEX IF NOT EXISTS motif_registry_id ON motif_registry(motif_id)")
    meta = dict(conn.execute("SELECT key,value FROM meta"))
    m = json.loads(meta.get("counts", "{}"))
    m["motif_registry"] = conn.execute("SELECT COUNT(*) FROM motif_registry").fetchone()[0]
    meta["counts"] = json.dumps(m)
    conn.execute("INSERT OR REPLACE INTO meta VALUES ('counts',?)", (meta["counts"],))
    conn.commit()
    print(json.dumps({"motifs": len(motifs), "registry_rows": m["motif_registry"]}, indent=2))
    conn.close()


if __name__ == "__main__":
    main()
