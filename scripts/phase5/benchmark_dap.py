#!/usr/bin/env python3
"""Evaluate ReguTrace promoter predictions against O'Malley 2016 DAP-seq peaks.

Positive set: genes whose 1,000 bp upstream promoter interval overlaps at least
one assay peak. Unpeaked genes are UNKNOWN, not established negatives. The
reported "PPV_assay" and AUPRC therefore use unobserved as non-positive and are
optimistic/ambiguous; they are comparative diagnostics, not biological FDR.
"""
import argparse
import gzip
import json
import math
import os
import re
from collections import defaultdict

TF_TO_FILE = {
    "MYB44": "GSE60141_DAPSeq-MYB_tnt-MYB44_col-chr1-5_GEM_events.narrowPeak.gz",
    "WRKY21": "GSE60141_DAPSeq-WRKY_tnt-WRKY21_col-chr1-5_GEM_events.narrowPeak.gz",
    "ERF15": "GSE60141_DAPSeq-AP2EREBP_tnt-ERF15_col-chr1-5_GEM_events.narrowPeak.gz",
    "GBF3": "GSE60141_DAPSeq-bZIP_tnt-GBF3_col-chr1-5_GEM_events.narrowPeak.gz",
    "TGA1": "GSE60141_DAPSeq-bZIP_tnt-TGA1_col-chr1-5_GEM_events.narrowPeak.gz",
    "TCP24": "GSE60141_DAPSeq-TCP_tnt-TCP24_col-chr1-5_GEM_events.narrowPeak.gz",
    "PDF2": "GSE60141_DAPSeq-Homeobox_tnt-PDF2_col-chr1-5_GEM_events.narrowPeak.gz",
    "RAV1": "GSE60141_DAPSeq-RAV_tnt-RAV1_col-chr1-5_GEM_events.narrowPeak.gz",
    "HSFA6B": "GSE60141_DAPSeq-HSF_tnt-HSFA6B_col-chr1-5_GEM_events.narrowPeak.gz",
}


def norm_chrom(value):
    value = value.strip()
    value = re.sub(r"^chr", "", value, flags=re.I)
    return value.lstrip("0") or "0"


def read_promoters(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for row in f:
            fields = row.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[0] == "gene_id":
                continue
            out[fields[0]] = {
                "chrom": norm_chrom(fields[1]),
                "start": int(fields[2]),
                "end": int(fields[3]),
                "strand": fields[4],
                "length": int(fields[5]),
            }
    return out


def read_peaks(path):
    by_chrom = defaultdict(list)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        for row in f:
            fields = row.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = norm_chrom(fields[0])
            start, end = int(fields[1]), int(fields[2])
            by_chrom[chrom].append((start, end))
    for values in by_chrom.values():
        values.sort()
    return by_chrom


def overlaps(values, start, end):
    # Sorted interval list; linear scan is acceptable for these peak sizes.
    for a, b in values:
        if a > end:
            break
        if b >= start:
            return True
    return False


def read_evidence(path):
    by_tf = defaultdict(lambda: defaultdict(lambda: {"best": 2.0, "sites": 0}))
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for row in f:
            fields = row.rstrip("\n").split("\t")
            if len(fields) <= idx["p_value"]:
                continue
            tf = fields[idx["tf_id"]]
            target = fields[idx["target_id"]].split("|", 1)[0]
            try:
                p = float(fields[idx["p_value"]])
            except ValueError:
                continue
            item = by_tf[tf][target]
            item["sites"] += 1
            if p < item["best"]:
                item["best"] = p
    return by_tf


def auprc(scores, labels):
    # scores descending; labels 1 positive.
    order = sorted(range(len(scores)), key=lambda i: (-scores[i], i))
    total_pos = sum(labels)
    if total_pos == 0:
        return None
    tp = 0
    prev_recall = 0.0
    ap = 0.0
    for rank, i in enumerate(order, 1):
        if labels[i]:
            tp += 1
            recall = tp / total_pos
            precision = tp / rank
            ap += (recall - prev_recall) * precision
            prev_recall = recall
    return ap


def hypergeom_sf(observed, draws, successes, population):
    if draws <= 0 or population <= 0:
        return 1.0
    upper = min(successes, draws)
    total = math.comb(population, draws)

    def log_comb(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    logs = []
    for k in range(observed, upper + 1):
        logs.append(log_comb(successes, k) + log_comb(population - successes, draws - k))
    if not logs:
        return 0.0
    m = max(logs)
    log_num = m + math.log(sum(math.exp(v - m) for v in logs))
    log_total = log_comb(population, draws)
    return min(1.0, math.exp(log_num - log_total))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--promoters", required=True)
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--tf-genes", required=True)
    ap.add_argument("--peaks-dir", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    promoters = read_promoters(args.promoters)
    universe = sorted(promoters)
    evidence = read_evidence(args.evidence)
    tf_symbols = {}
    for tf in evidence:
        # Evidence has tf_symbol in every row; re-read first row cheaply.
        with open(args.evidence, encoding="utf-8") as f:
            next(f)
            for row in f:
                fs = row.split("\t")
                if fs and fs[0] == tf:
                    tf_symbols[tf] = fs[1]
                    break

    results = []
    for tf, symbol in sorted(tf_symbols.items()):
        filename = TF_TO_FILE.get(symbol)
        if not filename:
            results.append({"tf_id": tf, "tf_symbol": symbol, "status": "no_peak_file"})
            continue
        peak_path = os.path.join(args.peaks_dir, filename)
        if not os.path.exists(peak_path):
            results.append({"tf_id": tf, "tf_symbol": symbol, "status": "peak_file_missing"})
            continue
        peaks = read_peaks(peak_path)
        positives = set()
        for gene in universe:
            prom = promoters[gene]
            if overlaps(peaks.get(prom["chrom"], []), prom["start"], prom["end"]):
                positives.add(gene)
        targets = evidence.get(tf, {})
        if not positives:
            results.append({"tf_id": tf, "tf_symbol": symbol, "status": "no_positive_genes",
                            "peaks": sum(len(v) for v in peaks.values())})
            continue
        pred = set(targets)
        tp = len(pred & positives)
        precision = tp / len(pred) if pred else 0.0
        recall = tp / len(positives)
        prevalence = len(positives) / len(universe)
        enrichment = precision / prevalence if prevalence > 0 else None
        # AUPRC: all scanned targets + unobserved universe with score 0.
        scores = [1.0 / targets[g]["best"] if g in targets else 0.0 for g in universe]
        labels = [1 if g in positives else 0 for g in universe]
        ap_value = auprc(scores, labels)
        top = {}
        ranked = sorted(targets, key=lambda g: (targets[g]["best"], g))
        for k in (50, 100, 500, 1000):
            hit = sum(1 for g in ranked[:k] if g in positives)
            top[str(k)] = hit / min(k, len(ranked)) if ranked else 0.0
        # Hypergeometric enrichment for predicted set vs universe.
        hg = hypergeom_sf(tp, len(pred), len(positives), len(universe))
        results.append({
            "tf_id": tf, "tf_symbol": symbol, "status": "ok",
            "peaks": sum(len(v) for v in peaks.values()),
            "positive_genes": len(positives), "universe": len(universe),
            "predicted_targets": len(pred), "tp": tp,
            "ppv_assay": precision, "recall": recall,
            "prevalence": prevalence, "enrichment": enrichment,
            "auprc": ap_value, "hypergeom_p": hg,
            "top50": top["50"], "top100": top["100"], "top500": top["500"], "top1000": top["1000"],
            "best_motif_p_median": None,
        })
    fields = ["tf_id", "tf_symbol", "status", "peaks", "positive_genes", "universe",
              "predicted_targets", "tp", "ppv_assay", "recall", "prevalence", "enrichment",
              "auprc", "hypergeom_p", "top50", "top100", "top500", "top1000"]
    with open(args.out_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(fields) + "\n")
        for r in results:
            f.write("\t".join("" if r.get(k) is None else str(r.get(k)) for k in fields) + "\n")
    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump({"benchmark": "O'Malley 2016 Arabidopsis DAP-seq",
                   "positive_rule": "promoter interval overlaps >=1 narrowPeak",
                   "negative_rule": "unobserved (not established negative)",
                   "results": results}, f, indent=2)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
