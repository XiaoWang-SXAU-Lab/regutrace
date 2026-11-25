#!/usr/bin/env python3
"""Build the final schema-v2 ReguTrace background in one command.

Order matters:
  1. build_background.py       canonical v1 hits/pairs + motif_target
  2. add_background_layers.py  motif_testability, promoter, scan_config
  3. rebuild_motif_registry.py registry from the exact scan motif library

The final manifest records the v2 hash and counts. Existing output files are
replaced by the underlying builders; this wrapper refuses to run if the output
directory does not exist.
"""
import argparse
import datetime
import hashlib
import json
import os
import sqlite3
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(script, args):
    cmd = [sys.executable, os.path.join(HERE, script)] + args
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coord-map", required=True)
    ap.add_argument("--base-map", required=True)
    ap.add_argument("--tf-genes", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--testability", required=True)
    ap.add_argument("--motifs", required=True)
    ap.add_argument("--pwm-hash", default="")
    ap.add_argument("--specificity-groups", default="")
    ap.add_argument("--promoter-summary", required=True)
    ap.add_argument("--scan-manifest", required=True)
    ap.add_argument("--policy-targets", required=True,
                    help="Phase 5 phase5-policy-targets.tsv")
    ap.add_argument("--policy-version", default="phase5-v1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--policy", default="best@0.8")
    ap.add_argument("--manifest", default="")
    args = ap.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if not os.path.isdir(out_dir):
        raise SystemExit(f"output directory does not exist: {out_dir}")
    manifest_path = args.manifest or os.path.join(out_dir, "background-build.json")
    start = time.time()

    build_args = [
        "--coord-map", args.coord_map, "--base-map", args.base_map,
        "--tf-genes", args.tf_genes, "--sites", args.sites,
        "--testability", args.testability, "--motifs", args.motifs,
        "--out", args.out, "--policy", args.policy,
    ]
    if args.pwm_hash:
        build_args += ["--pwm-hash", args.pwm_hash]
    if args.specificity_groups:
        build_args += ["--specificity-groups", args.specificity_groups]
    run("build_background.py", build_args)

    run("add_background_layers.py", [
        "--db", args.out, "--testability", args.testability,
        "--promoter-summary", args.promoter_summary,
        "--scan-manifest", args.scan_manifest,
        "--meta-out", manifest_path,
    ])

    registry_args = ["--db", args.out, "--motifs", args.motifs]
    if args.specificity_groups:
        registry_args += ["--specificity-groups", args.specificity_groups]
    run("rebuild_motif_registry.py", registry_args)

    run("add_policy_layers.py", [
        "--db", args.out, "--coord-map", args.coord_map,
        "--policy-targets", args.policy_targets,
        "--policy-version", args.policy_version, "--vacuum",
    ])

    conn = sqlite3.connect(args.out)
    counts = {
        name: conn.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
        for name in ("tf_catalogue", "tf_motif", "tf_motif_coordinate",
                     "tf_motif_selected", "tf_motif_policy", "tf_motif_policy_config",
                     "motif_target", "motif_registry", "motif_testability",
                     "promoter", "scan_config", "pairs", "hits")
    }
    version = conn.execute("PRAGMA user_version").fetchone()[0]
    conn.close()

    manifest = {}
    if os.path.exists(manifest_path):
        with open(manifest_path, encoding="utf-8") as f:
            manifest = json.load(f)
    manifest.update({
        "schema_version": f"regutrace-regulatory-background-db-v{version}",
        "database_file": os.path.basename(args.out),
        "database_bytes": os.path.getsize(args.out),
        "database_sha256": sha256(args.out),
        "counts": counts,
        "build_order": ["build_background.py", "add_background_layers.py",
                        "rebuild_motif_registry.py", "add_policy_layers.py"],
        "build_seconds": round(time.time() - start, 2),
        "completed_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    })
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps({"database_sha256": manifest["database_sha256"],
                      "database_bytes": manifest["database_bytes"],
                      "counts": counts, "manifest": manifest_path}, indent=2))


if __name__ == "__main__":
    main()
