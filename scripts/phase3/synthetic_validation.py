#!/usr/bin/env python3
"""Independent synthetic validation of the Phase 3 scanner and background model.

The Python implementation below is written from the scanner contract, not
imported from Rust:

* PWM probabilities are normalized with the same +1e-3 pseudocount used by
  `read_jaspar`.
* log-odds are log2(p / background), quantized with scale 10.
* the null distribution is an exact convolution over quantized scores.
* p(score) is the upper tail P(null >= observed).

It then compares those values with the CLI's site table on fixed windows of
random promoters. Because the CLI scans many overlapping windows, the test
reads the p-value of a pre-specified window rather than the minimum p-value.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import subprocess
from pathlib import Path

import numpy as np

BASES = "ACGT"
SITE_COLUMNS = [
    "pwm_hash",
    "motif_ids",
    "gene_id",
    "chrom",
    "promoter_start",
    "promoter_end",
    "promoter_strand",
    "hit_start",
    "hit_end",
    "hit_strand",
    "genomic_hit_start",
    "genomic_hit_end",
    "genomic_hit_strand",
    "score",
    "p_value",
    "tf_count",
    "scan_key",
]


# ---------------------------------------------------------------------------
# Independent scanner model


def normalize_probability_matrix(matrix):
    """Apply the reader's probability-matrix smoothing."""
    out = []
    for column in matrix:
        total = sum(column)
        out.append([(value + 1e-3) / (total + 4e-3) for value in column])
    return out


def log_odds(matrix, background):
    return [
        [math.log2(max(value, 1e-9) / max(background[base], 1e-9)) for base, value in enumerate(column)]
        for column in matrix
    ]


def quantize(log_odds_matrix, scale=10.0):
    return [[int(round(value * scale)) for value in row] for row in log_odds_matrix]


def null_distribution(log_odds_matrix, background, scale=10.0):
    rows = quantize(log_odds_matrix, scale)
    dist = {0: 1.0}
    for row in rows:
        nxt = {}
        for score, probability in dist.items():
            for base in range(4):
                new_score = score + row[base]
                nxt[new_score] = nxt.get(new_score, 0.0) + probability * background[base]
        dist = nxt
    tail = {}
    running = 0.0
    for score in sorted(dist, reverse=True):
        running += dist[score]
        tail[score] = running
    return rows, tail


def tail_arrays(tail):
    keys = sorted(tail)
    return keys, [tail[key] for key in keys]


def tail_pvalue(tail, score):
    # Smallest null score >= observed.
    keys, values = tail_arrays(tail)
    index = bisect.bisect_left(keys, score)
    return values[index] if index < len(values) else 0.0


def tail_pvalue_sorted(keys, values, score):
    index = bisect.bisect_left(keys, score)
    return values[index] if index < len(values) else 0.0


def score_sequence(seq, rows):
    if len(seq) < len(rows):
        return None
    best = None
    for offset in range(len(seq) - len(rows) + 1):
        score = 0
        valid = True
        for k, row in enumerate(rows):
            base = BASES.find(seq[offset + k].upper())
            if base < 0:
                valid = False
                break
            score += row[base]
        if not valid:
            continue
        if best is None or score > best[0]:
            best = (score, offset)
    return best


def score_window(seq, offset, rows):
    score = 0
    for k, row in enumerate(rows):
        base = BASES.find(seq[offset + k].upper())
        if base < 0:
            return None
        score += row[base]
    return score


def information_content(matrix):
    return sum(
        2.0 + sum(value * math.log2(value) for value in column if value > 0)
        for column in matrix
    )


def base_frequencies_gc(gc):
    at = (1.0 - gc) / 2.0
    return [at, gc / 2.0, gc / 2.0, at]


# ---------------------------------------------------------------------------
# Motif and sequence generation


def motif_columns_from_consensus(consensus, strength, gc_preference=False):
    columns = []
    for letter in consensus:
        column = [0.05, 0.05, 0.05, 0.05]
        index = BASES.index(letter)
        column[index] = strength
        rest = (1.0 - strength) / 3.0
        for base in range(4):
            if base != index:
                column[base] = rest
        columns.append(column)
    return columns


def write_motifs(path, motifs):
    lines = []
    for motif_id, (consensus, strength) in motifs.items():
        columns = motif_columns_from_consensus(consensus, strength)
        lines.append(f">{motif_id} synthetic")
        for base in range(4):
            values = " ".join(f"{column[base]:.8f}" for column in columns)
            lines.append(f"{BASES[base]} [ {values} ]")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        motif_id: normalize_probability_matrix(motif_columns_from_consensus(consensus, strength))
        for motif_id, (consensus, strength) in motifs.items()
    }


def write_fasta(path, sequences):
    with open(path, "w", encoding="utf-8") as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f">s{index:05d}\n{sequence}\n")


def random_sequences(rng, count, length, frequencies):
    draws = rng.choice(4, size=(count, length), p=frequencies)
    return ["".join(BASES[index] for index in row) for row in draws]


# ---------------------------------------------------------------------------
# CLI interaction


def run_cli(binary, work, sequences, motifs, tag, background=None, empirical=False):
    fasta = work / f"{tag}.fa"
    motif_path = work / f"{tag}-motifs.txt"
    map_path = work / f"{tag}-map.tsv"
    sites_path = work / f"{tag}-sites.tsv"
    testability_path = work / f"{tag}-testability.tsv"
    evidence_path = work / f"{tag}-evidence.tsv"
    report_path = work / f"{tag}-report.json"
    manifest_path = work / f"{tag}-manifest.json"
    write_fasta(fasta, sequences)
    matrices = write_motifs(motif_path, motifs)
    with open(map_path, "w", encoding="utf-8") as handle:
        handle.write("tf_id\ttf_symbol\tfamily\tmotif_id\n")
        for motif_id in motifs:
            handle.write(f"{motif_id}_TF\t{motif_id}_TF\tsynthetic\t{motif_id}\n")
    command = [
        str(binary),
        "scan",
        "--promoters",
        str(fasta),
        "--motifs",
        str(motif_path),
        "--tf-map",
        str(map_path),
        "--pvalue",
        "1.0",
        "--one-strand",
        "--reuse-exact-pwm",
        "--sites-out",
        str(sites_path),
        "--motif-testability-out",
        str(testability_path),
        "--scan-manifest",
        str(manifest_path),
        "--report-json",
        str(report_path),
        "--threads",
        "1",
        "--out",
        str(evidence_path),
    ]
    if empirical:
        command += ["--background-mode", "empirical"]
    else:
        command += ["--background", ",".join(f"{value:.12f}" for value in background)]
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"scan failed for {tag}: {completed.stderr[-2000:] or completed.stdout[-2000:]}"
        )
    rows = []
    with open(sites_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    return matrices, rows, manifest


def fixed_window_pvalues(rows, motif_id, offset, length, sequence_ids):
    found = {}
    for row in rows:
        if row["motif_ids"] != motif_id:
            continue
        if int(row["hit_start"]) != offset + 1 or int(row["hit_end"]) != offset + length:
            continue
        gene = row["gene_id"]
        if gene in sequence_ids:
            found[gene] = float(row["p_value"])
    missing = [gene for gene in sequence_ids if gene not in found]
    return [found[gene] for gene in sequence_ids if gene in found], missing


def ks_uniform(values):
    values = sorted(values)
    n = len(values)
    if n == 0:
        return None, None
    d = max(
        max(abs((index + 1) / n - value), abs(index / n - value))
        for index, value in enumerate(values)
    )
    # Asymptotic Kolmogorov tail; adequate for n >= 100.
    scaled = math.sqrt(n) * d
    if scaled < 1e-12:
        return d, 1.0
    total = 0.0
    for k in range(1, 101):
        total += 2.0 * (-1) ** (k - 1) * math.exp(-2.0 * k * k * scaled * scaled)
    return d, min(1.0, max(0.0, total))


def binomial_tail(n, k, probability):
    """P(X >= k) for X ~ Binomial(n, probability)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    total = 0.0
    for value in range(k, n + 1):
        total += math.comb(n, value) * probability**value * (1 - probability) ** (n - value)
    return min(1.0, total)


def calibration_summary(pvalues):
    n = len(pvalues)
    out = {"n": n}
    for threshold in (0.05, 0.01, 1e-4):
        if n == 0:
            out[f"fraction_le_{threshold:g}"] = None
            out[f"binomial_p_le_{threshold:g}"] = None
            continue
        k = sum(1 for value in pvalues if value <= threshold)
        out[f"fraction_le_{threshold:g}"] = k / n
        out[f"binomial_p_le_{threshold:g}"] = binomial_tail(n, k, threshold)
    return out


# ---------------------------------------------------------------------------
# Validation blocks


def exact_reproduction(work, binary):
    """Compare independent Python p-values with CLI site p-values."""
    rng = np.random.default_rng(20261125)
    motifs = {
        "M_AT6": ("ACGTAC", 0.70),
        "M_GC6": ("GCGCGC", 0.70),
        "M_LOW4": ("ACGT", 0.40),
        "M_HIGH10": ("ACGTACGTAC", 0.80),
    }
    sequences = random_sequences(rng, 1000, 100, [0.25] * 4)
    ids = [f"s{index:05d}" for index in range(len(sequences))]
    offset = 20
    results = []
    for background_name, background in [
        ("uniform", [0.25, 0.25, 0.25, 0.25]),
        ("at_rich", base_frequencies_gc(0.30)),
    ]:
        tag = f"reproduce_{background_name}"
        matrices, rows, manifest = run_cli(
            binary, work, sequences, motifs, tag, background=background
        )
        cli_testability = {}
        with open(work / f"{tag}-testability.tsv", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="	"):
                cli_testability[row["motif_id"]] = row
        for motif_id, matrix in matrices.items():
            length = len(matrix)
            log_odds_matrix = log_odds(matrix, background)
            qrows, tail = null_distribution(log_odds_matrix, background)
            tail_keys, tail_values = tail_arrays(tail)
            min_score = sum(max(row) for row in qrows)
            min_p = tail_pvalue_sorted(tail_keys, tail_values, min_score)
            cli_values, missing = fixed_window_pvalues(
                rows, motif_id, offset, length, ids
            )
            python_values = []
            for sequence in sequences:
                score = score_window(sequence, offset, qrows)
                python_values.append(tail_pvalue_sorted(tail_keys, tail_values, score))
            max_abs = max(
                (abs(a - b) for a, b in zip(cli_values, python_values)), default=None
            )
            results.append(
                {
                    "test": "exact_reproduction",
                    "background": background_name,
                    "motif": motif_id,
                    "length": length,
                    "information_content": information_content(matrix),
                    "min_pvalue": min_p,
                    "cli_min_pvalue": float(cli_testability[motif_id]["min_pvalue"]),
                    "n": len(cli_values),
                    "missing_fixed_windows": len(missing),
                    "max_abs_pvalue_difference": max_abs,
                    "passed": (
                        max_abs is not None
                        and max_abs < 1e-12
                        and not missing
                        and abs(min_p - float(cli_testability[motif_id]["min_pvalue"])) < 1e-12
                    ),
                }
            )
    return results, motifs, sequences, ids, offset


def null_calibration(work, binary, motifs, sequences, ids, offset):
    """Fixed-window p-values must be uniform under the generating background."""
    results = []
    for motif_id in ("M_AT6", "M_GC6"):
        length = len(motifs[motif_id][0])
        matrices, rows, manifest = run_cli(
            binary,
            work,
            sequences,
            {motif_id: motifs[motif_id]},
            f"null_{motif_id}",
            background=[0.25] * 4,
        )
        pvalues, missing = fixed_window_pvalues(rows, motif_id, offset, length, ids)
        d, _ = ks_uniform(pvalues)
        results.append(
            {
                "test": "null_calibration",
                "motif": motif_id,
                "missing": len(missing),
                "ks_d_discrete": d,
                **calibration_summary(pvalues),
            }
        )
    return results


def background_calibration(work, binary):
    """GC-rich promoters with a GC-rich motif: uniform vs empirical background."""
    rng = np.random.default_rng(20261126)
    motif = {"M_GC6": ("GCGCGC", 0.70)}
    length = 6
    offset = 20
    results = []
    for gc in (0.20, 0.50, 0.80):
        frequencies = base_frequencies_gc(gc)
        sequences = random_sequences(rng, 1000, 100, frequencies)
        ids = [f"s{index:05d}" for index in range(len(sequences))]
        for mode in ("uniform", "empirical"):
            tag = f"gc{int(gc*100)}_{mode}"
            matrices, rows, manifest = run_cli(
                binary,
                work,
                sequences,
                motif,
                tag,
                background=[0.25] * 4,
                empirical=(mode == "empirical"),
            )
            matrix = matrices["M_GC6"]
            pvalues, missing = fixed_window_pvalues(rows, "M_GC6", offset, length, ids)
            d, _ = ks_uniform(pvalues)
            results.append(
                {
                    "test": "gc_calibration",
                    "gc_content": gc,
                    "background_mode": mode,
                    "estimated_acgt": manifest["background"]["values_acgt"],
                    "missing": len(missing),
                    "ks_d_discrete": d,
                    **calibration_summary(pvalues),
                }
            )
    return results


def attainable_floor(work, binary, motifs, sequences, ids):
    """Minimum attainable p-value by motif length and IC."""
    length_motifs = {
        f"L{length}": ("ACGTACGTACGTACGTACGT"[:length], 0.75)
        for length in (4, 6, 8, 10, 12, 16, 20)
    }
    matrices, rows, manifest = run_cli(
        binary,
        work,
        sequences[:50],
        length_motifs,
        "floor",
        background=[0.25] * 4,
    )
    # The testability table is the authoritative CLI floor.
    testability_path = work / "floor-testability.tsv"
    cli = {}
    with open(testability_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cli[row["motif_id"]] = row
    results = []
    for motif_id, matrix in matrices.items():
        background = [0.25] * 4
        log_odds_matrix = log_odds(matrix, background)
        qrows, tail = null_distribution(log_odds_matrix, background)
        tail_keys, tail_values = tail_arrays(tail)
        min_score = sum(max(row) for row in qrows)
        python_min = tail_pvalue_sorted(tail_keys, tail_values, min_score)
        results.append(
            {
                "test": "attainable_floor",
                "motif": motif_id,
                "length": len(matrix),
                "information_content": information_content(matrix),
                "python_min_pvalue": python_min,
                "cli_min_pvalue": float(cli[motif_id]["min_pvalue"]),
                "attainable_at_1e-4": python_min <= 1e-4,
                "status": cli[motif_id]["status"],
            }
        )
    return results


def overlap_dependence():
    rng = np.random.default_rng(20261127)
    sequence = "".join(BASES[index] for index in rng.choice(4, size=5000))
    background = [0.25] * 4
    results = []
    for consensus in ("ACGTAC", "AAAAAA"):
        matrix = normalize_probability_matrix(
            motif_columns_from_consensus(consensus, 0.70)
        )
        log_odds_matrix = log_odds(matrix, background)
        qrows = quantize(log_odds_matrix)
        scores = [
            score_window(sequence, offset, qrows)
            for offset in range(len(sequence) - len(qrows) + 1)
        ]
        correlations = {}
        for lag in (1, 2, 3, 5, 10):
            a = np.array(scores[:-lag], dtype=float)
            b = np.array(scores[lag:], dtype=float)
            correlations[f"lag_{lag}"] = float(np.corrcoef(a, b)[0, 1])
        results.append(
            {
                "test": "overlap_dependence",
                "motif": consensus,
                "sequence_length": len(sequence),
                "motif_length": len(qrows),
                "windows": len(scores),
                "score_correlations": correlations,
                "max_abs_lag1": abs(correlations["lag_1"]),
            }
        )
    return results


def discretization_error():
    rng = np.random.default_rng(20261128)
    results = []
    for consensus, strength in [("ACGTAC", 0.70), ("ACGT", 0.40), ("ACGTACGT", 0.60)]:
        length = len(consensus)
        matrix = normalize_probability_matrix(motif_columns_from_consensus(consensus, strength))
        background = [0.25] * 4
        continuous_odds = [
            [math.log2(max(value, 1e-9) / 0.25) for value in column]
            for column in matrix
        ]
        qrows, tail = null_distribution(continuous_odds, background)
        tail_keys, tail_values = tail_arrays(tail)
        # Enumerate all sequences of this length under the uniform background.
        sequences = []
        for index in range(4**length):
            digits = []
            value = index
            for _ in range(length):
                digits.append(value % 4)
                value //= 4
            sequences.append(digits[::-1])
        probabilities = []
        scores_discrete = []
        scores_continuous = []
        for digits in sequences:
            probability = 1.0
            discrete = 0
            continuous = 0.0
            for position, base in enumerate(digits):
                probability *= background[base]
                discrete += qrows[position][base]
                continuous += continuous_odds[position][base]
            probabilities.append(probability)
            scores_discrete.append(discrete)
            scores_continuous.append(continuous)
        order = np.argsort(scores_continuous)
        continuous_sorted = [scores_continuous[i] for i in order]
        probability_sorted = [probabilities[i] for i in order]
        tail_continuous = {}
        running = 0.0
        for score, probability in zip(reversed(continuous_sorted), reversed(probability_sorted)):
            running += probability
            tail_continuous[score] = running
        continuous_keys = sorted(tail_continuous)
        continuous_values = [tail_continuous[key] for key in continuous_keys]
        max_difference = 0.0
        for discrete, continuous in zip(scores_discrete, scores_continuous):
            discrete_p = tail_pvalue_sorted(tail_keys, tail_values, discrete)
            index = bisect.bisect_left(continuous_keys, continuous - 1e-12)
            continuous_p = (
                continuous_values[index] if index < len(continuous_values) else 0.0
            )
            max_difference = max(max_difference, abs(discrete_p - continuous_p))
        results.append(
            {
                "test": "discretization",
                "motif": consensus,
                "length": length,
                "max_abs_pvalue_difference": max_difference,
                "enumerated_sequences": 4**length,
            }
        )
    return results


def write_tables(out, results):
    rows = []
    for result in results:
        row = dict(result)
        row.pop("score_correlations", None)
        row.pop("estimated_acgt", None)
        rows.append(row)
    columns = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with open(out / "synthetic-validation.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)

    results = []
    reproduction, motifs, sequences, ids, offset = exact_reproduction(
        args.workdir, args.binary
    )
    results.extend(reproduction)
    results.extend(null_calibration(args.workdir, args.binary, motifs, sequences, ids, offset))
    results.extend(background_calibration(args.workdir, args.binary))
    results.extend(attainable_floor(args.workdir, args.binary, motifs, sequences[:50], ids[:50]))
    results.extend(overlap_dependence())
    results.extend(discretization_error())
    write_tables(args.out, results)
    (args.out / "synthetic-validation.json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
