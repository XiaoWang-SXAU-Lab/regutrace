#!/usr/bin/env python3
"""Add exact policy selections and the coordinate rerun layer to a ReguTrace DB.

The policy table stores the exact TF–motif association selected by the Phase 5
audit. Target sets are reconstructed with relational joins; no policy-specific
TF–target table is materialized.

Inputs:
  --db             existing v2 background
  --coord-map      coordinate rerun map (the data layer used by Phase 5)
  --policy-targets phase5-policy-targets.tsv from the audited run
  --policy-version provenance string for the selection rules

The script is idempotent: it replaces the three policy/coordinate tables.
"""
import argparse
import csv
import hashlib
import json
import os
import sqlite3


POLICY_LABELS = [("all", "all"), ("best", "best@0.8"), ("reps", "reps@0.8")]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_policy_targets(path):
    out = {name: {} for name, _ in POLICY_LABELS}
    by_label = {label: name for name, label in POLICY_LABELS}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            name = by_label.get(row.get("policy", ""))
            if not name:
                continue
            motifs = [m for m in (row.get("motifs") or "").split(";") if m]
            out[name][row["tf_id"]] = motifs
    for name, rows in out.items():
        if not rows:
            raise SystemExit(f"policy {name} is missing from {path}")
    return out


def read_keyed(path, key_fields=("tf_id", "motif_id")):
    out = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = tuple(row.get(field, "") for field in key_fields)
            out.setdefault(key, row)
    return out


def create_coordinate_table(conn, path):
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames or []
        if "tf_id" not in header or "motif_id" not in header:
            raise SystemExit("coordinate map needs tf_id and motif_id")
        cols = ", ".join('"' + name.replace('"', '""') + '" TEXT' for name in header)
        conn.execute("DROP TABLE IF EXISTS tf_motif_coordinate")
        conn.execute(f"CREATE TABLE tf_motif_coordinate ({cols})")
        placeholders = ",".join("?" for _ in header)
        batch = []
        count = 0
        for row in reader:
            batch.append(tuple(row.get(name, "") for name in header))
            count += 1
            if len(batch) >= 5000:
                conn.executemany(f"INSERT INTO tf_motif_coordinate VALUES ({placeholders})", batch)
                batch = []
        if batch:
            conn.executemany(f"INSERT INTO tf_motif_coordinate VALUES ({placeholders})", batch)
    conn.execute("CREATE INDEX IF NOT EXISTS tmc_tf_motif ON tf_motif_coordinate(tf_id, motif_id)")
    return count


def selection_group(row):
    if not row:
        return ""
    start = row.get("query_full_dbd_start", "")
    end = row.get("query_full_dbd_end", "")
    return "|".join([
        row.get("extraction_family", ""),
        row.get("domain_mode", ""),
        row.get("domain_count", ""),
        f"{start}-{end}" if start or end else "",
    ])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--coord-map", required=True)
    ap.add_argument("--policy-targets", required=True)
    ap.add_argument("--policy-version", default="phase5-v1")
    ap.add_argument("--meta-out", default="")
    ap.add_argument("--vacuum", action="store_true")
    args = ap.parse_args()

    policies = read_policy_targets(args.policy_targets)
    base = {}
    conn = sqlite3.connect(args.db, timeout=120)
    for row in conn.execute(
        "SELECT tf_id,motif_id,mapping_status,domain_mode,extraction_family,domain_count,"
        "reference_tf_id,specificity_groups,dbd_identity,query_dbd_coverage,"
        "donor_dbd_coverage,dbd_bitscore,query_full_dbd_start,query_full_dbd_end FROM tf_motif"
    ):
        tf, motif, status, mode, family, domain_count, ref, spec, identity, qcov, dcov, bits, qs, qe = row
        base.setdefault((tf, motif), {
            "tf_id": tf, "motif_id": motif, "mapping_status": status or "",
            "domain_mode": mode or "", "extraction_family": family or "",
            "domain_count": domain_count or "", "reference_tf_id": ref or "",
            "specificity_groups": spec or "", "dbd_identity": identity or "",
            "query_dbd_coverage": qcov or "", "donor_dbd_coverage": dcov or "",
            "dbd_bitscore": bits or "", "query_full_dbd_start": qs or "",
            "query_full_dbd_end": qe or "",
        })
    coord = read_keyed(args.coord_map)
    registry = {}
    for motif, group in conn.execute("SELECT motif_id, specificity_group_080 FROM motif_registry"):
        registry[motif] = group or ""

    coord_count = create_coordinate_table(conn, args.coord_map)

    conn.execute("DROP TABLE IF EXISTS tf_motif_policy")
    conn.execute("""
        CREATE TABLE tf_motif_policy (
            policy TEXT NOT NULL,
            policy_version TEXT NOT NULL,
            tf_id TEXT NOT NULL,
            motif_id TEXT NOT NULL,
            selection_role TEXT NOT NULL,
            selection_source TEXT NOT NULL,
            selection_group TEXT,
            specificity_group TEXT,
            reference_tf_id TEXT,
            selection_order INTEGER NOT NULL,
            PRIMARY KEY (policy, tf_id, motif_id)
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS tmpolicy_tf ON tf_motif_policy(policy, tf_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS tmpolicy_motif ON tf_motif_policy(policy, motif_id)")

    conn.execute("DROP TABLE IF EXISTS tf_motif_policy_config")
    conn.execute("""
        CREATE TABLE tf_motif_policy_config (
            policy TEXT PRIMARY KEY,
            policy_version TEXT NOT NULL,
            parameters_json TEXT NOT NULL,
            source_file TEXT NOT NULL,
            source_sha256 TEXT NOT NULL,
            coordinate_source TEXT NOT NULL,
            coordinate_sha256 TEXT NOT NULL,
            row_count INTEGER NOT NULL
        )
    """)

    base_params = {
        "overlap_threshold": 0.8,
        "ranking": ["full_architecture", "min(query_dbd_coverage,donor_dbd_coverage)",
                    "dbd_identity", "dbd_bitscore", "motif_id"],
        "strict_calibrated_retained": True,
        "ambiguous_included": False,
        "uncoordinated_exploratory": "excluded from best/reps when any coordinate rows exist",
        "no_coordinate_fallback": "single best exploratory motif per TF",
        "coordinate_layer": "tf_motif_coordinate",
    }
    params = {
        "all": {
            "candidate_set": "all base-map rows",
            "ambiguous_included": True,
            "single_domain_fallback_included": True,
            "ranking": None,
        },
        "best": dict(base_params),
        "reps": dict(base_params, representative="one per specificity group; motif_id when absent"),
    }
    source_hash = sha256(args.policy_targets)
    coord_hash = sha256(args.coord_map)
    rows = []
    for name, _ in POLICY_LABELS:
        policy_rows = policies[name]
        for tf, motifs in sorted(policy_rows.items()):
            for order, motif in enumerate(motifs):
                b = base.get((tf, motif))
                c = coord.get((tf, motif))
                fixed = bool(b and b["mapping_status"] in ("strict_rbh", "dbd_calibrated"))
                role = "fixed" if fixed else "exploratory"
                if c:
                    source = "coordinate_rerun"
                    group = selection_group(c)
                    ref = c.get("reference_tf_id", "")
                    spec = c.get("specificity_groups", "")
                elif b:
                    source = "base_map"
                    group = selection_group(b)
                    ref = b.get("reference_tf_id", "")
                    spec = b.get("specificity_groups", "")
                else:
                    raise SystemExit(f"{name} selection {tf}/{motif} is in neither base nor coordinate map")
                if not spec:
                    spec = registry.get(motif, "")
                rows.append((name, args.policy_version, tf, motif, role, source,
                             group, spec, ref, order))
    conn.executemany("INSERT INTO tf_motif_policy VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
    counts = {}
    for name, _ in POLICY_LABELS:
        counts[name] = sum(1 for row in rows if row[0] == name)
        conn.execute(
            "INSERT INTO tf_motif_policy_config VALUES (?,?,?,?,?,?,?,?)",
            (name, args.policy_version, json.dumps(params[name], sort_keys=True),
             os.path.basename(args.policy_targets), source_hash,
             os.path.basename(args.coord_map), coord_hash, counts[name]),
        )

    meta = dict(conn.execute("SELECT key,value FROM meta"))
    table_counts = {name: conn.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
                    for name in ("tf_catalogue", "tf_motif", "tf_motif_coordinate",
                                 "tf_motif_selected", "tf_motif_policy",
                                 "tf_motif_policy_config", "motif_target",
                                 "motif_registry", "motif_testability", "promoter",
                                 "scan_config", "pairs", "hits")}
    tables = json.loads(meta.get("tables", "[]"))
    for name in ("tf_motif_coordinate", "tf_motif_policy", "tf_motif_policy_config"):
        if name not in tables:
            tables.append(name)
    meta["schema_version"] = "regutrace-regulatory-background-db-v3"
    meta["tables"] = json.dumps(tables)
    meta["counts"] = json.dumps(table_counts)
    meta["policy_version"] = args.policy_version
    meta["policy_sources"] = json.dumps({
        "policy_targets": {"file": os.path.basename(args.policy_targets), "sha256": source_hash},
        "coordinate_rerun": {"file": os.path.basename(args.coord_map), "sha256": coord_hash},
    }, sort_keys=True)
    for key, value in meta.items():
        conn.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", (key, value))
    conn.execute("PRAGMA user_version=3")
    conn.commit()
    if args.vacuum:
        conn.execute("VACUUM")
    conn.close()

    summary = {
        "schema_version": "regutrace-regulatory-background-db-v3",
        "policy_version": args.policy_version,
        "coordinate_rows": coord_count,
        "policy_rows": counts,
        "counts": table_counts,
        "policy_targets_sha256": source_hash,
        "coordinate_sha256": coord_hash,
    }
    if args.meta_out:
        with open(args.meta_out, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
