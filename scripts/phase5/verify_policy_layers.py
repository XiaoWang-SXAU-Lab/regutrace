#!/usr/bin/env python3
"""Verify exact policy TF sets and TF-target pair sets in a v3 background.

The reference selection is phase5-policy-targets.tsv; the reference target
space is the DB motif_target table. For each TF we compare the expected target
set with the DB policy join as sets, not only as aggregate counts.
"""
import argparse
import collections
import csv
import hashlib
import json
import sqlite3


POLICIES = [("all", "all"), ("best", "best@0.8"), ("reps", "reps@0.8")]


def read_policy_targets(path):
    out = {name: {} for name, _ in POLICIES}
    by_label = {label: name for name, label in POLICIES}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            name = by_label.get(row.get("policy", ""))
            if not name:
                continue
            out[name][row["tf_id"]] = [m for m in (row.get("motifs") or "").split(";") if m]
    return out


def read_summary(path):
    if not path:
        return {}
    data = json.load(open(path, encoding="utf-8"))
    return data.get("data", data)


def target_query(conn, policy, tf, use_pairs=False):
    if use_pairs:
        rows = conn.execute("SELECT target_id FROM pairs WHERE tf_id=?", (tf,))
    else:
        rows = conn.execute(
            "SELECT DISTINCT mt.target_id FROM tf_motif_policy p "
            "JOIN motif_target mt ON mt.motif_id=p.motif_id "
            "WHERE p.policy=? AND p.tf_id=?",
            (policy, tf),
        )
    return {row[0] for row in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--policy-targets", required=True)
    ap.add_argument("--summary", default="")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    conn = sqlite3.connect(args.db, timeout=120)
    expected = read_policy_targets(args.policy_targets)
    summary = read_summary(args.summary)
    motif_targets = collections.defaultdict(list)
    for motif, target in conn.execute("SELECT motif_id, target_id FROM motif_target ORDER BY motif_id"):
        motif_targets[motif].append(target)

    db_selection = {name: collections.defaultdict(list) for name, _ in POLICIES}
    for policy, tf, motif in conn.execute(
        "SELECT policy, tf_id, motif_id FROM tf_motif_policy ORDER BY policy, tf_id, selection_order"
    ):
        db_selection[policy][tf].append(motif)

    report = {
        "schema_version": conn.execute("SELECT value FROM meta WHERE key='schema_version'").fetchone()[0],
        "policies": {},
    }
    for name, label in POLICIES:
        expected_tf = expected[name]
        db_tf = db_selection[name]
        selection_mismatch = []
        for tf in sorted(set(expected_tf) | set(db_tf)):
            if set(expected_tf.get(tf, [])) != set(db_tf.get(tf, [])):
                selection_mismatch.append(tf)
        covered = 0
        pairs = 0
        pair_mismatch = []
        global_hash = hashlib.sha256()
        pairs_table_mismatch = []
        for tf in sorted(expected_tf):
            targets = set()
            for motif in expected_tf[tf]:
                targets.update(motif_targets.get(motif, ()))
            actual = target_query(conn, name, tf)
            if targets != actual:
                pair_mismatch.append({
                    "tf_id": tf,
                    "expected_targets": len(targets),
                    "actual_targets": len(actual),
                    "missing": len(targets - actual),
                    "extra": len(actual - targets),
                })
            if targets:
                covered += 1
                pairs += len(targets)
                global_hash.update((tf + "\n").encode())
                for target in sorted(targets):
                    global_hash.update((target + "\n").encode())
            if name == "best":
                materialized = target_query(conn, name, tf, use_pairs=True)
                if materialized != targets:
                    pairs_table_mismatch.append(tf)
        expected_summary = summary.get(f"{label}|original", {}) if isinstance(summary, dict) else {}
        expected_covered = expected_summary.get("tfs_with_targets")
        expected_pairs = expected_summary.get("total_tf_target_pairs")
        report["policies"][name] = {
            "selection_tfs": len(expected_tf),
            "selection_motifs": sum(len(v) for v in expected_tf.values()),
            "selection_mismatch_tfs": selection_mismatch[:20],
            "selection_mismatch_count": len(selection_mismatch),
            "covered_tfs": covered,
            "pairs": pairs,
            "pair_set_sha256": global_hash.hexdigest(),
            "pair_mismatch_tfs": pair_mismatch[:20],
            "pair_mismatch_count": len(pair_mismatch),
            "pairs_table_mismatch_count": len(pairs_table_mismatch) if name == "best" else None,
            "expected_summary_covered": expected_covered,
            "expected_summary_pairs": expected_pairs,
            "counts_match": (expected_covered in (None, covered)
                             and expected_pairs in (None, pairs)),
        }
    conn.close()
    ok = all(
        item["selection_mismatch_count"] == 0
        and item["pair_mismatch_count"] == 0
        and (item["pairs_table_mismatch_count"] in (None, 0))
        and item["counts_match"]
        for item in report["policies"].values()
    )
    report["ok"] = ok
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
