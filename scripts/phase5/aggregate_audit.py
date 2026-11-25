#!/usr/bin/env python3
"""Coordinate-aware motif grouping, representative policies, and HF1 coverage.

Inputs are existing Phase 3/4 artifacts; no promoter rescan is performed.
This script does not alter mapping status or evidence. All target relations
remain level D.
"""
import argparse
import array
import csv
import json
import math
import os
import pickle
import sys
import traceback
from collections import defaultdict

COMPOSITE_FAMILIES = {"RAV"}


def read_coord_map(path):
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    out = defaultdict(list)
    for row in rows:
        def num(key):
            try:
                return float(row.get(key) or 0)
            except ValueError:
                return 0.0
        row["_qs"] = num("query_full_dbd_start")
        row["_qe"] = num("query_full_dbd_end")
        row["_qcov"] = num("query_dbd_coverage")
        row["_dcov"] = num("donor_dbd_coverage")
        row["_id"] = num("dbd_identity")
        row["_bits"] = num("dbd_bitscore")
        row["_rank"] = (
            1.0 if row.get("domain_mode") == "full_architecture" else 0.0,
            min(row["_qcov"], row["_dcov"]),
            row["_id"],
            row["_bits"],
        )
        out[row["tf_id"]].append(row)
    return out


def read_pwm_hash(path):
    hashes = {}
    if not path or not os.path.exists(path):
        return hashes
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            hashes[row["motif_id"]] = row.get("pwm_hash") or row["motif_id"]
    return hashes


def read_specificity_groups(path):
    groups = {}
    if not path or not os.path.exists(path):
        return groups
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            groups[row["motif_id"]] = row.get("group_080") or row["motif_id"]
    return groups


def overlap_coeff(a, b):
    if not a or not b:
        return 0.0
    sa, ea = a
    sb, eb = b
    inter = min(ea, eb) - max(sa, sb) + 1
    if inter <= 0:
        return 0.0
    la, lb = ea - sa + 1, eb - sb + 1
    return inter / min(la, lb)


def group_rows(rows, threshold):
    """Union-find over compatible architecture + query DBD overlap.

    Composite architectures (e.g. RAV) are not split by non-overlapping local
    alignments; they remain one functional architecture group.
    """
    if not rows:
        return []
    base = defaultdict(list)
    for row in rows:
        family = row.get("extraction_family") or row.get("family") or ""
        mode = row.get("domain_mode") or ""
        key = (family, mode, row.get("domain_count") or "")
        base[key].append(row)
    groups = []
    for key, part in base.items():
        family, mode, _ = key
        if family in COMPOSITE_FAMILIES:
            groups.append(part)
            continue
        parent = list(range(len(part)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(i, j):
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[rj] = ri

        for i in range(len(part)):
            for j in range(i + 1, len(part)):
                if overlap_coeff((part[i]["_qs"], part[i]["_qe"]), (part[j]["_qs"], part[j]["_qe"])) >= threshold:
                    union(i, j)
        by_root = defaultdict(list)
        for i, row in enumerate(part):
            by_root[find(i)].append(row)
        groups.extend(by_root.values())
    return groups


def dedup_exact(rows, pwm_hash):
    """Keep the best row per exact PWM matrix; exact duplicates are redundant."""
    best = {}
    for row in rows:
        key = pwm_hash.get(row["motif_id"], row["motif_id"])
        current = best.get(key)
        if current is None or row["_rank"] > current["_rank"]:
            best[key] = row
    return list(best.values())


def select_policies(rows, pwm_hash, specificity):
    policies = {"all": rows, "best": [], "reps": []}
    for threshold in (0.5, 0.7, 0.8, 0.9):
        groups = group_rows(rows, threshold)
        best = []
        reps = []
        for group in groups:
            unique = dedup_exact(group, pwm_hash)
            if not unique:
                continue
            ordered = sorted(unique, key=lambda r: (r["_rank"], r["motif_id"]), reverse=True)
            top = ordered[0]
            best.append(top)
            # One representative per specificity group, chosen by the same
            # sequence/alignment ranking. Ties remain represented by the
            # deterministic motif-ID order and are not merged into a new call.
            seen = set()
            for row in ordered:
                group_id = specificity.get(row["motif_id"], row["motif_id"])
                if group_id in seen:
                    continue
                seen.add(group_id)
                reps.append(row)
        policies[f"best@{threshold}"] = best
        policies[f"reps@{threshold}"] = reps
    return policies


def motif_info(path):
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                out[row["motif_id"]] = {
                    "length": int(float(row.get("matrix_length") or 0)),
                    "status": row.get("status") or "",
                    "attainable": (row.get("attainable") or "").lower() == "true",
                    "hits": int(float(row.get("motif_id_hit_count") or 0)),
                }
            except ValueError:
                pass
    return out


def bits_set(bits, index):
    bits[index >> 3] |= 1 << (index & 7)


def bits_or(dst, src):
    if len(dst) < len(src):
        dst.extend(b"\x00" * (len(src) - len(dst)))
    for i, value in enumerate(src):
        dst[i] |= value


def bits_count(bits):
    return sum(value.bit_count() for value in bits)


def read_sites(path, motif_info_map):
    """One pass: motif-target bitsets, site p arrays, and per-motif target best p."""
    target_index = {}
    bits = defaultdict(bytearray)
    pvals = defaultdict(lambda: array.array("d"))
    tidx = defaultdict(lambda: array.array("I"))
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= idx["p_value"] or len(fields) <= idx["motif_ids"]:
                continue
            target = fields[idx["gene_id"]]
            try:
                p = float(fields[idx["p_value"]])
            except ValueError:
                continue
            ti = target_index.get(target)
            if ti is None:
                ti = len(target_index)
                target_index[target] = ti
            motif_ids = [m for m in fields[idx["motif_ids"]].split(";") if m]
            for motif in motif_ids:
                b = bits[motif]
                if len(b) <= ti >> 3:
                    b.extend(b"\x00" * ((ti >> 3) + 1 - len(b)))
                bits_set(b, ti)
                pvals[motif].append(p)
                tidx[motif].append(ti)
    return target_index, bits, pvals, tidx


def site_q_for_motif(pvals, motif_length, promoter_count, promoter_len, alpha=0.05):
    """Conservative BH over hits ranked within the complete tested-window universe.

    Non-hit windows have p > requested cutoff; omitting them cannot make a hit
    q smaller than the exact BH q, so this is a conservative (upper-bound) q.
    """
    n = len(pvals)
    if n == 0:
        return array.array("B", b"\x00" * n)
    windows = max(1, 2 * (promoter_len - motif_length + 1))
    total = promoter_count * windows
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        value = min(1.0, total * pvals[i] / rank)
        if value < running:
            running = value
        q[i] = running
    return array.array("B", bytes(1 if value <= alpha else 0 for value in q))


def promoter_q_pass(tidx, pvals, motif_length, promoter_count, promoter_len, alpha=0.05):
    """Bonferroni promoter p from best site p and actual tested windows, then BH
    across all eligible promoters. Non-hit promoters are omitted from the rank
    list (their p is 1), which is conservative for hit q-values.
    """
    best = {}
    for ti, p in zip(tidx, pvals):
        if ti not in best or p < best[ti]:
            best[ti] = p
    if not best:
        return set()
    windows = max(1, 2 * (promoter_len - motif_length + 1))
    values = [(min(1.0, p * windows), ti) for ti, p in best.items()]
    values.sort()
    q = [1.0] * len(values)
    running = 1.0
    for rank in range(len(values), 0, -1):
        p, _ = values[rank - 1]
        value = min(1.0, promoter_count * p / rank)
        if value < running:
            running = value
        q[rank - 1] = running
    return {ti for (p, ti), qv in zip(values, q) if qv <= alpha}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coord-map", required=True)
    ap.add_argument("--base-map", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--testability", required=True)
    ap.add_argument("--pwm-hash", default="")
    ap.add_argument("--specificity-groups", default="")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--promoter-count", type=int, default=37415)
    ap.add_argument("--promoter-length", type=int, default=1000)
    ap.add_argument("--state", default="")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    coord = read_coord_map(args.coord_map)
    base = read_coord_map(args.base_map)
    pwm_hash = read_pwm_hash(args.pwm_hash)
    specificity = read_specificity_groups(args.specificity_groups)
    info = motif_info(args.testability)

    state_path = args.state or os.path.join(args.out_dir, "aggregate-state.pkl")
    if state_path and os.path.exists(state_path):
        with open(state_path, "rb") as f:
            target_index, motif_bits, pvals, tidx = pickle.load(f)
        print(f"loaded state targets={len(target_index)} motifs={len(motif_bits)}", flush=True)
    else:
        target_index, motif_bits, pvals, tidx = read_sites(args.sites, info)
        print(f"site targets={len(target_index)} motifs={len(motif_bits)}", flush=True)
        if state_path:
            with open(state_path, "wb") as f:
                pickle.dump((target_index, dict(motif_bits), dict(pvals), dict(tidx)), f, protocol=4)
            print(f"wrote state {state_path}", flush=True)

    # Site- and promoter-level FDR target bitsets.
    site_pass = {}
    promoter_pass = {}
    for n, (motif, values) in enumerate(pvals.items(), 1):
        meta = info.get(motif, {})
        length = meta.get("length") or 6
        sq = site_q_for_motif(values, length, args.promoter_count, args.promoter_length)
        sb = bytearray((len(target_index) + 7) // 8)
        for i, ok in enumerate(sq):
            if ok:
                bits_set(sb, tidx[motif][i])
        site_pass[motif] = sb
        pr = promoter_q_pass(tidx[motif], values, length, args.promoter_count, args.promoter_length)
        pb = bytearray((len(target_index) + 7) // 8)
        for ti in pr:
            bits_set(pb, ti)
        promoter_pass[motif] = pb
        if n % 25 == 0 or n == len(pvals):
            print(f"fdr motif {n}/{len(pvals)} {motif} len={length} hits={len(values)} site_pass={bits_count(sb)} promoter_pass={bits_count(pb)}", flush=True)

    audit_rows = []
    target_rows = []
    all_tfs = sorted(base)
    for tfn, tf in enumerate(all_tfs, 1):
        if tfn % 100 == 0:
            print(f"policy TF {tfn}/{len(all_tfs)}", flush=True)
        rows = coord.get(tf, [])
        brows = base.get(tf, [])
        for threshold in (0.5, 0.7, 0.8, 0.9):
            groups = group_rows(rows, threshold) if rows else []
            audit_rows.append({
                "tf_id": tf, "tf_symbol": (rows or brows)[0].get("tf_symbol", "") if (rows or brows) else "",
                "family": (rows or brows)[0].get("family", "") if (rows or brows) else "",
                "mapping_rows": len(rows), "unique_motifs": len({r["motif_id"] for r in rows}),
                "unique_pwms": len({pwm_hash.get(r["motif_id"], r["motif_id"]) for r in rows}),
                "specificity_groups_080": len({specificity.get(r["motif_id"], r["motif_id"]) for r in rows}),
                "overlap_threshold": threshold, "architecture_groups": len(groups),
                "largest_group_rows": max((len(g) for g in groups), default=0),
                "query_dbd_min": min((r["_qs"] for r in rows), default=0),
                "query_dbd_max": max((r["_qe"] for r in rows), default=0),
                "composite": bool(rows) and rows[0].get("extraction_family") in COMPOSITE_FAMILIES,
            })
        # Strict and calibrated mappings are not DBD exploratory candidates;
        # keep them unchanged and only apply coordinate grouping to exploratory
        # rows. Ambiguous exploratory rows remain in the all-candidate view.
        fixed = [r for r in brows if r.get("mapping_status") in ("strict_rbh", "dbd_calibrated")]
        base_all = brows
        if not rows:
            # A DBD TF that was not reproduced by the coordinate rerun: use the
            # original rows as opaque singleton groups for best/reps.
            rows = []
        policies = {"all": base_all}
        for threshold in (0.5, 0.7, 0.8, 0.9):
            best = list(fixed)
            reps = list(fixed)
            if rows:
                expl = [r for r in rows if r.get("mapping_status") == "dbd_exploratory"]
                for group in group_rows(expl, threshold):
                    unique = dedup_exact(group, pwm_hash)
                    if not unique:
                        continue
                    ordered = sorted(unique, key=lambda r: (r["_rank"], r["motif_id"]), reverse=True)
                    best.append(ordered[0])
                    seen = set()
                    for row in ordered:
                        gid = specificity.get(row["motif_id"], row["motif_id"])
                        if gid in seen:
                            continue
                        seen.add(gid)
                        reps.append(row)
            else:
                # No coordinate information: conservative single best row from
                # the original exploratory rows, deterministic by motif ID.
                expl = [r for r in brows if r.get("mapping_status") == "dbd_exploratory"]
                if expl:
                    selected = sorted(expl, key=lambda r: (float(r.get("query_dbd_coverage") or 0), r["motif_id"]), reverse=True)[0]
                    best.append(selected); reps.append(selected)
            policies[f"best@{threshold}"] = best
            policies[f"reps@{threshold}"] = reps
        for name, selected in policies.items():
            motifs = [r["motif_id"] for r in selected]
            row = {"tf_id": tf, "policy": name, "selected_motifs": len(set(motifs)),
                   "selected_pwms": len({pwm_hash.get(m, m) for m in motifs}),
                   "selected_specificity_groups": len({specificity.get(m, m) for m in motifs}),
                   "selected_untestable": sum(1 for m in set(motifs) if not info.get(m, {}).get("attainable", True)),
                   "motifs": ";".join(sorted(set(motifs)))}
            for label, source in (("original", motif_bits), ("site_fdr", site_pass), ("promoter_fdr", promoter_pass)):
                bits = bytearray((len(target_index) + 7) // 8)
                for motif in set(motifs):
                    bits_or(bits, source.get(motif, b""))
                count = bits_count(bits)
                row[f"targets_{label}"] = count
                row[f"pairs_{label}"] = count
            target_rows.append(row)

    with open(os.path.join(args.out_dir, "phase5-coord-audit.tsv"), "w", newline="", encoding="utf-8") as f:
        fields = list(audit_rows[0].keys())
        w = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(audit_rows)
    with open(os.path.join(args.out_dir, "phase5-policy-targets.tsv"), "w", newline="", encoding="utf-8") as f:
        fields = list(target_rows[0].keys())
        w = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(target_rows)
    summary = {}
    for name in ("all", "best@0.5", "best@0.7", "best@0.8", "best@0.9",
                 "reps@0.5", "reps@0.7", "reps@0.8", "reps@0.9"):
        vals = [r for r in target_rows if r["policy"] == name]
        for label in ("original", "site_fdr", "promoter_fdr"):
            counts = sorted(r[f"targets_{label}"] for r in vals)
            total = sum(counts)
            covered = sum(1 for v in counts if v > 0)
            q1 = counts[len(counts) // 4] if counts else 0
            q3 = counts[(3 * len(counts)) // 4] if counts else 0
            summary[f"{name}|{label}"] = {
                "tfs": len(counts), "tfs_with_targets": covered,
                "coverage": covered / 1907, "median_targets": (counts[len(counts)//2] if counts else 0),
                "q1_targets": q1, "q3_targets": q3, "total_tf_target_pairs": total,
                "median_motifs": sorted(r["selected_motifs"] for r in vals)[len(vals)//2],
                "median_pwms": sorted(r["selected_pwms"] for r in vals)[len(vals)//2],
                "selected_untestable": sum(r["selected_untestable"] for r in vals),
            }
    json.dump({"data": summary, "site_targets": len(target_index),
               "target_index_size": len(target_index)}, open(os.path.join(args.out_dir, "phase5-policy-summary.json"), "w"), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
