#!/usr/bin/env python3
"""Render Phase 2.6 ablation figures as standalone SVG."""
import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "phase2_5"))
from common import read_tsv  # noqa: E402
from make_figures import (  # noqa: E402
    INK,
    MUTED,
    GRID,
    BLUE,
    ORANGE,
    circle,
    line,
    rect,
    square,
    svg_document,
    text,
    triangle,
)

METHOD_COLORS = {
    "A": "#0072B2",
    "B": "#D55E00",
    "C": "#009E73",
    "D": "#CC79A7",
    "E": "#E69F00",
}
METHOD_LABELS = {
    "A": "A full-length strict RBH",
    "B": "B DBD reciprocal",
    "C": "C DBD NN global threshold",
    "D": "D DBD NN family thresholds",
    "E": "E DBD specificity-group",
}


def marker(x, y, method, size=5):
    color = METHOD_COLORS[method]
    if method == "A":
        return circle(x, y, size, color)
    if method == "B":
        return square(x, y, size * 1.7, color)
    if method == "C":
        return triangle(x, y, size, color)
    if method == "D":
        points = (
            f"{x:.1f},{y-size:.1f} {x+size:.1f},{y:.1f} "
            f"{x:.1f},{y+size:.1f} {x-size:.1f},{y:.1f}"
        )
        return f'<polygon points="{points}" fill="{color}"/>'
    points = " ".join(
        f"{x + size * 0.95 * math.cos(2 * math.pi * k / 5 - math.pi / 2):.1f},"
        f"{y + size * 0.95 * math.sin(2 * math.pi * k / 5 - math.pi / 2):.1f}"
        for k in range(5)
    )
    return f'<polygon points="{points}" fill="{color}"/>'


def axes(width, height, margins, x_range, y_range, x_label, y_label, title, subtitle=None):
    left, right, top, bottom = margins
    plot_w = width - left - right
    plot_h = height - top - bottom
    body = [text(left, 28, title, 16, weight="600")]
    if subtitle:
        body.append(text(left, 47, subtitle, 11, fill=MUTED))
    for tick in range(6):
        frac = tick / 5
        x = left + frac * plot_w
        value = x_range[0] + frac * (x_range[1] - x_range[0])
        body.append(line(x, top, x, top + plot_h, GRID, 1))
        body.append(text(x, top + plot_h + 17, f"{value:g}", 10, anchor="middle", fill=MUTED))
        y = top + plot_h - frac * plot_h
        y_value = y_range[0] + frac * (y_range[1] - y_range[0])
        body.append(line(left, y, left + plot_w, y, GRID, 1))
        body.append(text(left - 7, y + 3, f"{y_value:g}", 10, anchor="end", fill=MUTED))
    body.append(line(left, top + plot_h, left + plot_w, top + plot_h, INK, 1))
    body.append(line(left, top, left, top + plot_h, INK, 1))
    body.append(text(left + plot_w / 2, height - 18, x_label, 11, anchor="middle", fill=MUTED))
    body.append(text(14, top + plot_h / 2, y_label, 11, anchor="middle", fill=MUTED, rotate=-90))
    return body, left, top, plot_w, plot_h


def xmap(value, x_range, left, plot_w):
    return left + (value - x_range[0]) / (x_range[1] - x_range[0]) * plot_w


def ymap(value, y_range, top, plot_h):
    return top + plot_h - (value - y_range[0]) / (y_range[1] - y_range[0]) * plot_h


def figure_precision_coverage(methods, bootstrap, path):
    boot = {row["method"]: row for row in bootstrap}
    width, height = 1000, 700
    body, left, top, plot_w, plot_h = axes(
        width,
        height,
        (85, 250, 70, 70),
        (0, 0.8),
        (0, 1),
        "coverage (accepted calls / donors)",
        "precision (motif similarity >= 0.80)",
        "Mapping strategy precision vs coverage",
        "Species-held-out nested evaluation; bars are 95% homology-cluster bootstrap intervals.",
    )
    for row in methods:
        method = row["method"]
        x = xmap(float(row["coverage"]), (0, 0.8), left, plot_w)
        y = ymap(float(row["precision"]), (0, 1), top, plot_h)
        if method in boot and boot[method]["coverage_ci_low"]:
            x1 = xmap(float(boot[method]["coverage_ci_low"]), (0, 0.8), left, plot_w)
            x2 = xmap(float(boot[method]["coverage_ci_high"]), (0, 0.8), left, plot_w)
            y1 = ymap(float(boot[method]["precision_ci_low"]), (0, 1), top, plot_h)
            y2 = ymap(float(boot[method]["precision_ci_high"]), (0, 1), top, plot_h)
            body.append(line(x1, y, x2, y, METHOD_COLORS[method], 2))
            body.append(line(x, y1, x, y2, METHOD_COLORS[method], 2))
        body.append(marker(x, y, method, 6))
        body.append(text(x + 9, y + 4, method, 13, weight="600", fill=METHOD_COLORS[method]))
    legend_x = width - 235
    body.append(text(legend_x, 100, "Method", 12, weight="600"))
    for index, method in enumerate(("A", "B", "C", "D", "E")):
        y = 122 + index * 22
        body.append(marker(legend_x, y, method, 5))
        body.append(text(legend_x + 14, y + 4, METHOD_LABELS[method], 10))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def figure_validation(validation, path):
    methods = ("A", "B", "C", "D", "E")
    schemes = ("donor_held_out", "species_held_out", "homology_cluster_held_out")
    labels = {
        "donor_held_out": "donor",
        "species_held_out": "species",
        "homology_cluster_held_out": "homology cluster",
    }
    lookup = {(row["fold_scheme"], row["method"]): row for row in validation}
    width, height = 1100, 540
    body = [text(30, 28, "Validation setting changes the apparent performance", 16, weight="600")]
    body.append(text(30, 47, "Same calibration data and models; only the held-out grouping changes.", 11, fill=MUTED))
    for panel, metric in enumerate(("precision", "coverage")):
        left = 70 + panel * 520
        top, plot_w, plot_h = 90, 420, 330
        body.append(text(left, top - 15, metric.capitalize(), 13, weight="600"))
        for tick in (0, 0.25, 0.5, 0.75, 1.0):
            y = top + plot_h - tick * plot_h
            body.append(line(left, y, left + plot_w, y, GRID, 1))
            body.append(text(left - 7, y + 3, f"{tick:g}", 9, anchor="end", fill=MUTED))
        body.append(line(left, top + plot_h, left + plot_w, top + plot_h, INK, 1))
        body.append(line(left, top, left, top + plot_h, INK, 1))
        for index, scheme in enumerate(schemes):
            x = left + 55 + index * 125
            body.append(text(x, top + plot_h + 20, labels[scheme], 10, anchor="middle", fill=MUTED))
            for offset, method in enumerate(methods):
                row = lookup.get((scheme, method))
                if not row:
                    continue
                value = float(row[metric])
                bar_x = x - 45 + offset * 17
                body.append(rect(bar_x, top + plot_h - value * plot_h, 14, value * plot_h, METHOD_COLORS[method], rx=2))
        if panel == 0:
            body.append(text(left + plot_w / 2, height - 22, "held-out unit", 10, anchor="middle", fill=MUTED))
    for index, method in enumerate(methods):
        x = 70 + index * 170
        body.append(marker(x, 480, method, 5))
        body.append(text(x + 12, 484, f"{method} {METHOD_LABELS[method].split(' ', 1)[1]}", 9))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def figure_bootstrap(rows, path):
    width, height = 1000, 480
    left, plot_w = 250, 560
    body = [text(30, 28, "Pairwise AUC depends on the resampling unit", 16, weight="600")]
    body.append(text(30, 47, "Same 11,784 pairs; 300 bootstrap replicates per unit. Bars are 95% intervals.", 11, fill=MUTED))
    for tick in (0.55, 0.60, 0.65, 0.70, 0.75, 0.80):
        x = xmap(tick, (0.55, 0.80), left, plot_w)
        body.append(line(x, 80, x, 360, GRID, 1))
        body.append(text(x, 375, f"{tick:.2f}", 9, anchor="middle", fill=MUTED))
    labels = [
        ("dbd_identity", "pair", "DBD / pairs"),
        ("dbd_identity", "donor", "DBD / donors"),
        ("dbd_identity", "cluster", "DBD / homology clusters"),
        ("protein_identity", "pair", "protein / pairs"),
        ("protein_identity", "donor", "protein / donors"),
        ("protein_identity", "cluster", "protein / homology clusters"),
    ]
    for index, (score, unit, label) in enumerate(labels):
        row = next(r for r in rows if r["score"] == score and r["bootstrap_unit"] == unit)
        y = 105 + index * 42
        auc_value = float(row["auc"])
        low = float(row["bootstrap_ci_low"])
        high = float(row["bootstrap_ci_high"])
        color = BLUE if score == "dbd_identity" else ORANGE
        body.append(text(20, y + 4, label, 10))
        body.append(line(xmap(low, (0.55, 0.80), left, plot_w), y, xmap(high, (0.55, 0.80), left, plot_w), y, color, 3))
        body.append(circle(xmap(auc_value, (0.55, 0.80), left, plot_w), y, 4, color))
        body.append(text(xmap(high, (0.55, 0.80), left, plot_w) + 6, y + 4, f"{auc_value:.3f}", 9, fill=MUTED))
    body.append(text(left + plot_w / 2, height - 20, "AUC", 10, anchor="middle", fill=MUTED))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ablation", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    figure_precision_coverage(
        read_tsv(args.ablation / "ablation-methods.tsv"),
        read_tsv(args.audit / "method-cluster-bootstrap.tsv"),
        args.out / "strategy-precision-coverage.svg",
    )
    figure_validation(
        read_tsv(args.audit / "validation-setting-comparison.tsv"),
        args.out / "validation-settings.svg",
    )
    figure_bootstrap(
        read_tsv(args.audit / "pairwise-bootstrap-comparison.tsv"),
        args.out / "pairwise-bootstrap.svg",
    )
    print("wrote 3 SVG figures")


if __name__ == "__main__":
    main()
