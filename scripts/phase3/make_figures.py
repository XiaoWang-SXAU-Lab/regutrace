#!/usr/bin/env python3
"""Render Phase 3 figures as standalone SVG without external plotting libraries."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "phase2_5"))
from make_figures import (  # noqa: E402
    GRID,
    INK,
    MUTED,
    circle,
    line,
    rect,
    svg_document,
    text,
)

BLUE = "#2F6FB0"
ORANGE = "#D9822B"
GRAY = "#8A8A8A"


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:

        return list(csv.DictReader(handle, delimiter="\t"))


def bar_panel(body, left, top, width, height, title, labels, values, colors, unit):
    body.append(text(left, top - 26, title, 13, weight="600"))
    maximum = (max(values) if values else 1) * 1.12
    axis_top = top + 18
    plot_h = height - 54
    body.append(line(left, axis_top + plot_h, left + width, axis_top + plot_h, INK, 1))
    for tick in range(5):
        frac = tick / 4
        y = axis_top + plot_h - frac * plot_h
        body.append(line(left, y, left + width, y, GRID, 1))
        body.append(
            text(
                left - 7,
                y + 3,
                f"{frac * maximum:,.0f}",
                9,
                anchor="end",
                fill=MUTED,
            )
        )
    bar_width = width / max(1, len(values)) * 0.55
    step = width / max(1, len(values))
    for index, (label, value, color) in enumerate(zip(labels, values, colors)):
        x = left + step * (index + 0.5)
        h = value / maximum * plot_h if maximum else 0
        body.append(rect(x - bar_width / 2, axis_top + plot_h - h, bar_width, h, color, rx=4))
        body.append(text(x, axis_top + plot_h - h - 7, f"{value:,.0f}", 10, anchor="middle", weight="600"))
        body.append(text(x, axis_top + plot_h + 16, label, 9, anchor="middle", fill=MUTED))


def figure_architecture(benchmark, memory, path):
    width, height = 1000, 520
    body = [text(40, 32, "Scan architecture: same evidence, much lower memory", 16, weight="600")]
    body.append(
        text(
            40,
            51,
            "Run 1 (legacy, sorted in memory) and Run 2 (unique-PWM groups, streamed) use the same 2024 map and thresholds.",
            11,
            fill=MUTED,
        )
    )
    legacy = json.loads(Path(memory["legacy"]).read_text(encoding="utf-8"))
    optimized = json.loads(Path(memory["optimized"]).read_text(encoding="utf-8"))
    bar_panel(
        body,
        90,
        105,
        380,
        330,
        "Peak working set (MiB)",
        ["legacy", "optimized"],
        [legacy["peak_working_set_mib"], optimized["peak_working_set_mib"]],
        [BLUE, ORANGE],
        "MiB",
    )
    rows = {row["run"]: row for row in read_tsv(benchmark)}
    run1 = rows["run1_2024_legacy"]
    run2 = rows["run2_2024_reuse"]
    bar_panel(
        body,
        560,
        105,
        380,
        330,
        "Wall time (same scan, seconds)",
        [
            "legacy cold",
            "optimized cold",
            "optimized warm",
        ],
        [
            float(legacy["wall_seconds"]),
            float(run2["scan_runtime_seconds"]),
            float(optimized["wall_seconds"]),
        ],
        [BLUE, ORANGE, "#E8B87A"],
        "seconds",
    )
    body.append(
        text(
            40,
            485,
            "Evidence rows: Run 1 3,090,137 vs Run 2 3,090,137. Full-row set difference including notes: 0.",
            10,
            fill=MUTED,
        )
    )
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def figure_background(benchmark, path):
    rows = {row["run"]: row for row in read_tsv(benchmark)}
    uniform = rows["run3_2026_uniform"]
    empirical = rows["run4_2026_empirical"]
    panels = [
        ("Unique PWM sites", "unique_site_hits", "sites"),
        ("TF–target pairs", "tf_target_pairs", "pairs"),
        ("Testable mapped motifs", "testable_motifs", "motifs"),
        ("Unattainable motifs", "unattainable_motifs", "motifs"),
    ]
    width, height = 1000, 620
    body = [text(40, 32, "Empirical promoter background changes the scan, not the mapping", 16, weight="600")]
    body.append(
        text(
            40,
            51,
            "Same JASPAR 2026 strict map and thresholds. The estimated background is 65.1% AT (A/T 0.3256, C/G 0.1744).",
            11,
            fill=MUTED,
        )
    )
    for index, (title, key, unit) in enumerate(panels):
        left = 70 + (index % 2) * 480
        top = 115 + (index // 2) * 245
        bar_panel(
            body,
            left,
            top,
            390,
            170,
            title,
            ["uniform", "empirical"],
            [float(uniform[key]), float(empirical[key])],
            [BLUE, ORANGE],
            unit,
        )
    body.append(circle(70, 590, 5, BLUE))
    body.append(text(82, 594, "uniform 0.25 each", 10))
    body.append(circle(210, 590, 5, ORANGE))
    body.append(text(222, 594, "empirical promoter model", 10))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def figure_floor(testability, path):
    rows = [
        row
        for row in read_tsv(testability)
        if row["mapped"] == "true"
    ]
    width, height = 1000, 600
    left, top, plot_w, plot_h = 95, 95, 800, 430
    body = [text(40, 32, "Attainable p-value floor depends on motif length and information content", 16, weight="600")]
    body.append(
        text(
            40,
            51,
            "Each point is one mapped JASPAR 2026 motif. Points below the dashed line cannot reach p ≤ 1e-4 for any promoter.",
            11,
            fill=MUTED,
        )
    )
    min_log = -16.5
    max_log = 0.0
    body.append(line(left, top + plot_h, left + plot_w, top + plot_h, INK, 1))
    body.append(line(left, top, left, top + plot_h, INK, 1))
    for log_value in (-16, -12, -8, -4, 0):
        y = top + plot_h - (log_value - min_log) / (max_log - min_log) * plot_h
        body.append(line(left, y, left + plot_w, y, GRID, 1))
        body.append(text(left - 8, y + 3, f"1e{log_value}", 9, anchor="end", fill=MUTED))
    threshold_y = top + plot_h - (math.log10(1e-4) - min_log) / (max_log - min_log) * plot_h
    body.append(line(left, threshold_y, left + plot_w, threshold_y, GRAY, 2, dash="7 5"))
    body.append(text(left + plot_w - 8, threshold_y - 7, "requested p = 1e-4", 10, anchor="end", fill=MUTED))
    lengths = sorted({int(row["matrix_length"]) for row in rows})
    for length in lengths:
        x = left + (length - min(lengths)) / max(1, max(lengths) - min(lengths)) * plot_w
        body.append(text(x, top + plot_h + 18, str(length), 9, anchor="middle", fill=MUTED))
    body.append(text(left + plot_w / 2, height - 28, "matrix length (bp)", 11, anchor="middle", fill=MUTED))
    body.append(text(24, top + plot_h / 2, "minimum attainable p (log10)", 11, anchor="middle", fill=MUTED, rotate=-90))
    attainable = 0
    unattainable = 0
    for index, row in enumerate(rows):
        length = int(row["matrix_length"])
        value = float(row["min_pvalue"])
        x = left + (length - min(lengths)) / max(1, max(lengths) - min(lengths)) * plot_w
        # Deterministic jitter within each length column.
        offset = ((index * 37) % 100 - 50) / 100.0 * 8.0
        x += offset
        y = top + plot_h - (math.log10(max(value, 1e-300)) - min_log) / (max_log - min_log) * plot_h
        if value <= 1e-4:
            attainable += 1
            body.append(circle(x, y, 2.6, BLUE, opacity=0.72))
        else:
            unattainable += 1
            body.append(circle(x, y, 2.9, ORANGE, opacity=0.9))
    body.append(circle(650, 120, 5, BLUE))
    body.append(text(662, 124, f"attainable at 1e-4 ({attainable})", 10))
    body.append(circle(650, 140, 5, ORANGE))
    body.append(text(662, 144, f"unattainable ({unattainable})", 10))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def figure_calibration(synthetic, path):
    data = json.loads(Path(synthetic).read_text(encoding="utf-8"))
    rows = [row for row in data if row.get("test") == "gc_calibration"]
    uniform = {row["gc_content"]: row["fraction_le_0.05"] for row in rows if row["background_mode"] == "uniform"}
    empirical = {row["gc_content"]: row["fraction_le_0.05"] for row in rows if row["background_mode"] == "empirical"}
    width, height = 1000, 610
    left, top, plot_w, plot_h = 95, 95, 820, 380
    body = [text(40, 32, "GC-rich promoters expose the cost of a fixed uniform background", 16, weight="600")]
    body.append(
        text(
            40,
            51,
            "GC-rich motif (GCGCGC), fixed window, 1,000 iid promoters per GC level. Nominal false-positive rate 5%.",
            11,
            fill=MUTED,
        )
    )
    y_max = max(0.20, max(uniform.values()) * 1.15)
    for tick in range(6):
        frac = tick / 5
        y = top + plot_h - frac * plot_h
        body.append(line(left, y, left + plot_w, y, GRID, 1))
        body.append(text(left - 8, y + 3, f"{frac * y_max:.2f}", 9, anchor="end", fill=MUTED))
    body.append(line(left, top + plot_h, left + plot_w, top + plot_h, INK, 1))
    body.append(line(left, top, left, top + plot_h, INK, 1))
    xs = sorted(uniform)
    for gc in xs:
        x = left + (gc - min(xs)) / max(1e-9, max(xs) - min(xs)) * plot_w
        body.append(text(x, top + plot_h + 18, f"{gc:.1f}", 10, anchor="middle", fill=MUTED))
    body.append(text(left + plot_w / 2, top + plot_h + 42, "promoter GC content", 11, anchor="middle", fill=MUTED))
    body.append(text(24, top + plot_h / 2, "fraction p ≤ 0.05", 11, anchor="middle", fill=MUTED, rotate=-90))
    nominal_y = top + plot_h - 0.05 / y_max * plot_h
    body.append(line(left, nominal_y, left + plot_w, nominal_y, GRAY, 2, dash="7 5"))
    body.append(text(left + 8, nominal_y - 7, "nominal 0.05", 10, fill=MUTED))
    for values, color, label in [
        (uniform, BLUE, "uniform"),
        (empirical, ORANGE, "empirical"),
    ]:
        points = []
        for gc in xs:
            x = left + (gc - min(xs)) / max(1e-9, max(xs) - min(xs)) * plot_w
            y = top + plot_h - values[gc] / y_max * plot_h
            points.append((x, y))
        for (x1, y1), (x2, y2) in zip(points, points[1:]):
            body.append(line(x1, y1, x2, y2, color, 2))
        for index, (x, y) in enumerate(points):
            body.append(circle(x, y, 5, color))
            body.append(text(x, y - 10, f"{values[xs[index]]:.3f}", 9, anchor="middle", fill=MUTED))
    body.append(circle(95, 565, 5, BLUE))
    body.append(text(107, 569, "uniform 0.25 each", 10))
    body.append(circle(260, 565, 5, ORANGE))
    body.append(text(272, 569, "empirical promoter model", 10))
    body.append(
        text(
            40,
            592,
            "At 80% GC, the uniform model produces 16.9% false positives; the empirical model produces 4.7%.",
            10,
            fill=MUTED,
        )
    )
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--legacy-memory", required=True)
    parser.add_argument("--optimized-memory", required=True)
    parser.add_argument("--testability", required=True)
    parser.add_argument("--synthetic", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    figure_architecture(
        Path(args.benchmark),
        {"legacy": args.legacy_memory, "optimized": args.optimized_memory},
        out / "scan-architecture.svg",
    )
    figure_background(Path(args.benchmark), out / "background-effect.svg")
    figure_floor(Path(args.testability), out / "attainable-floor.svg")
    figure_calibration(Path(args.synthetic), out / "pvalue-calibration.svg")


if __name__ == "__main__":
    main()
