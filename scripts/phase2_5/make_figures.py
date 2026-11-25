#!/usr/bin/env python3
"""Render Phase 2.5 diagnostic SVG figures without external plotting libraries."""
import argparse
import collections
import html
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_tsv  # noqa: E402

INK = "#222222"
MUTED = "#5B5B5B"
GRID = "#E4E4E4"
BLUE = "#2F6FB0"
ORANGE = "#D9822B"
GREEN = "#4E8F5A"
GRAY = "#8A8A8A"
LIGHT_GRAY = "#C9C9C9"


def esc(text):
    return html.escape(str(text), quote=True)


def text(x, y, value, size=12, anchor="start", fill=INK, weight="normal", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Segoe UI, Arial, sans-serif" '
        f'font-size="{size}" text-anchor="{anchor}" fill="{fill}" font-weight="{weight}"'
        f"{transform}>{esc(value)}</text>"
    )


def line(x1, y1, x2, y2, stroke=INK, width=1, dash=None, opacity=1.0):
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
        f'stroke="{stroke}" stroke-width="{width}" opacity="{opacity}"{dash_attr}/>'
    )


def rect(x, y, width, height, fill, stroke="none", rx=0, opacity=1.0):
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
        f'rx="{rx}" fill="{fill}" stroke="{stroke}" opacity="{opacity}"/>'
    )


def circle(x, y, radius, fill, opacity=1.0):
    return (
        f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="{fill}" '
        f'opacity="{opacity}"/>'
    )


def square(x, y, size, fill, opacity=1.0):
    return rect(x - size / 2, y - size / 2, size, size, fill, opacity=opacity)


def triangle(x, y, size, fill, opacity=1.0):
    points = [
        (x, y - size),
        (x - size, y + size),
        (x + size, y + size),
    ]
    path = " ".join(f"{px:.2f},{py:.2f}" for px, py in points)
    return f'<polygon points="{path}" fill="{fill}" opacity="{opacity}"/>'


def svg_document(width, height, body):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        f'<rect width="{width}" height="{height}" fill="#FCFCFB"/>'
        + "".join(body)
        + "</svg>\n"
    )


def axes(width, height, margins, x_range, y_range, x_label, y_label, title, subtitle=None):
    left, right, top, bottom = margins
    plot_w = width - left - right
    plot_h = height - top - bottom
    body = [text(left, 26, title, 16, weight="600")]
    if subtitle:
        body.append(text(left, 45, subtitle, 11, fill=MUTED))
    for tick in range(6):
        fx = tick / 5
        x = left + fx * plot_w
        value = x_range[0] + fx * (x_range[1] - x_range[0])
        body.append(line(x, top, x, top + plot_h, GRID, 1))
        body.append(text(x, top + plot_h + 18, f"{value:g}", 10, anchor="middle", fill=MUTED))
    for tick in range(6):
        fy = tick / 5
        y = top + plot_h - fy * plot_h
        value = y_range[0] + fy * (y_range[1] - y_range[0])
        body.append(line(left, y, left + plot_w, y, GRID, 1))
        body.append(text(left - 8, y + 3, f"{value:g}", 10, anchor="end", fill=MUTED))
    body.append(line(left, top + plot_h, left + plot_w, top + plot_h, INK, 1.2))
    body.append(line(left, top, left, top + plot_h, INK, 1.2))
    body.append(text(left + plot_w / 2, height - 18, x_label, 12, anchor="middle", fill=MUTED))
    body.append(text(14, top + plot_h / 2, y_label, 12, anchor="middle", fill=MUTED, rotate=-90))
    return body, left, top, plot_w, plot_h


def xmap(value, x_range, left, plot_w):
    return left + (value - x_range[0]) / (x_range[1] - x_range[0]) * plot_w


def ymap(value, y_range, top, plot_h):
    return top + plot_h - (value - y_range[0]) / (y_range[1] - y_range[0]) * plot_h


def median_curve(pairs, x_range, left, top, plot_w, plot_h, bins=20, min_count=20):
    buckets = collections.defaultdict(list)
    for x, y in pairs:
        index = min(bins - 1, max(0, int((x - x_range[0]) / (x_range[1] - x_range[0]) * bins)))
        buckets[index].append(y)
    points = []
    for index in sorted(buckets):
        if len(buckets[index]) < min_count:
            continue
        x = x_range[0] + (index + 0.5) / bins * (x_range[1] - x_range[0])
        points.append((x, statistics.median(buckets[index])))
    if len(points) < 2:
        return []
    path = " ".join(
        f"{xmap(x, x_range, left, plot_w):.2f},{ymap(y, (0, 1), top, plot_h):.2f}"
        for x, y in points
    )
    return [f'<polyline points="{path}" fill="none" stroke="{INK}" stroke-width="2"/>']


def figure_overall(rows, path):
    data = [
        (
            float(row["dbd_identity"]),
            float(row["motif_similarity"]),
            row["species_relation"],
        )
        for row in rows
        if row.get("motif_similarity")
    ]
    width, height = 1200, 800
    body, left, top, plot_w, plot_h = axes(
        width,
        height,
        (80, 30, 70, 70),
        (0, 100),
        (0, 1),
        "DBD amino-acid identity (%)",
        "PWM similarity (max Pearson r)",
        "DBD identity versus PWM similarity",
        f"{len(data):,} within-family donor pairs; shapes and color distinguish species relationship",
    )
    for x, y, relation in data:
        if relation == "different_genus":
            body.append(triangle(xmap(x, (0, 100), left, plot_w), ymap(y, (0, 1), top, plot_h), 2.0, GREEN, 0.35))
    for x, y, relation in data:
        if relation == "same_genus":
            body.append(square(xmap(x, (0, 100), left, plot_w), ymap(y, (0, 1), top, plot_h), 3.4, ORANGE, 0.45))
    for x, y, relation in data:
        if relation == "same_species":
            body.append(circle(xmap(x, (0, 100), left, plot_w), ymap(y, (0, 1), top, plot_h), 1.8, BLUE, 0.35))
    body.extend(median_curve([(x, y) for x, y, _ in data], (0, 100), left, top, plot_w, plot_h))
    legend_x, legend_y = width - 280, 95
    for index, (label, color, shape) in enumerate(
        [("same species", BLUE, "circle"), ("same genus", ORANGE, "square"), ("different genus", GREEN, "triangle")]
    ):
        y = legend_y + index * 22
        if shape == "circle":
            body.append(circle(legend_x, y, 4, color))
        elif shape == "square":
            body.append(square(legend_x, y, 7, color))
        else:
            body.append(triangle(legend_x, y, 5, color))
        body.append(text(legend_x + 14, y + 4, label, 11))
    body.append(text(legend_x, legend_y + 78, "black line: 5%-bin median", 10, fill=MUTED))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def figure_families(rows, path, families):
    width, height = 1200, 760
    body = [text(30, 30, "Family-specific DBD identity versus PWM similarity", 16, weight="600")]
    body.append(text(30, 49, "Each panel uses within-family donor pairs; line is the 5%-bin median.", 11, fill=MUTED))
    panel_w, panel_h = 270, 300
    for index, family in enumerate(families):
        col = index % 4
        row = index // 4
        panel_left = 60 + col * 285
        panel_top = 80 + row * 330
        panel_body, left, top, plot_w, plot_h = axes(
            panel_left + panel_w,
            panel_top + panel_h,
            (panel_left, 20, panel_top, 45),
            (0, 100),
            (0, 1),
            "",
            "",
            family,
        )
        # axes() uses absolute width/height; crop by shifting is awkward, so draw directly.
        panel_body = []
        data = [
            (float(row["dbd_identity"]), float(row["motif_similarity"]))
            for row in rows
            if row["family"] == family and row.get("motif_similarity")
        ]
        left = panel_left + 42
        top = panel_top + 22
        plot_w = panel_w - 55
        plot_h = panel_h - 70
        for x, y in data:
            panel_body.append(circle(xmap(x, (0, 100), left, plot_w), ymap(y, (0, 1), top, plot_h), 1.4, BLUE, 0.32))
        panel_body.extend(median_curve(data, (0, 100), left, top, plot_w, plot_h, bins=10, min_count=5))
        panel_body.append(line(left, top + plot_h, left + plot_w, top + plot_h, INK, 1))
        panel_body.append(line(left, top, left, top + plot_h, INK, 1))
        for tick in (0, 50, 100):
            panel_body.append(text(xmap(tick, (0, 100), left, plot_w), top + plot_h + 14, str(tick), 9, anchor="middle", fill=MUTED))
        for tick in (0, 0.5, 1):
            panel_body.append(text(left - 6, ymap(tick, (0, 1), top, plot_h) + 3, f"{tick:g}", 9, anchor="end", fill=MUTED))
        panel_body.append(text(left, panel_top + 12, family, 12, weight="600"))
        panel_body.append(text(left, panel_top + 27, f"n={len(data):,}", 10, fill=MUTED))
        if row == 1:
            panel_body.append(text(left + plot_w / 2, panel_top + panel_h - 8, "DBD identity (%)", 9, anchor="middle", fill=MUTED))
        if col == 0:
            panel_body.append(text(left - 32, top + plot_h / 2, "PWM similarity", 9, anchor="middle", fill=MUTED, rotate=-90))
        body.extend(panel_body)
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def figure_family_sizes(summary, path):
    rows = sorted(summary, key=lambda row: (-int(row["donors"]), row["family"]))
    width = 900
    height = 70 + len(rows) * 23
    body = [text(30, 28, "Calibration donors per family", 16, weight="600")]
    body.append(text(30, 47, "Blue: sample-size supported (threshold-tested); gray: abstained before calibration.", 11, fill=MUTED))
    max_donors = max(int(row["donors"]) for row in rows)
    left, right, top = 200, 50, 65
    plot_w = width - left - right
    for index, row in enumerate(rows):
        y = top + index * 23
        donors = int(row["donors"])
        supported = row["calibration_status"] == "supported"
        bar_w = donors / max_donors * plot_w
        body.append(text(left - 8, y + 12, row["family"], 10, anchor="end"))
        body.append(rect(left, y + 2, bar_w, 14, BLUE if supported else LIGHT_GRAY, rx=3))
        body.append(text(left + bar_w + 5, y + 12, str(donors), 9, fill=MUTED))
    body.append(text(left + plot_w / 2, height - 16, "donors with DBD and experimentally supported motif", 10, anchor="middle", fill=MUTED))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def figure_precision_coverage(sweep, thresholds, families, path):
    width, height = 1200, 760
    body = [text(30, 30, "Precision versus coverage for candidate DBD thresholds", 16, weight="600")]
    body.append(text(30, 49, "Leave-one-species-out estimates; dashed line is the 0.90 precision target. Selected thresholds are marked.", 11, fill=MUTED))
    panel_w, panel_h = 370, 300
    for index, family in enumerate(families):
        col = index % 3
        row = index // 3
        left = 70 + col * 380
        top = 85 + row * 330
        plot_w, plot_h = 300, 230
        body.append(text(left, top - 12, family, 13, weight="600"))
        selected = next((row for row in thresholds if row["family"] == family), None)
        status = "calibrated" if selected else "abstained"
        body.append(text(left, top + 3, status, 10, fill=BLUE if selected else MUTED))
        for tick in (0, 0.5, 1):
            y = top + plot_h - tick * plot_h
            body.append(line(left, y, left + plot_w, y, GRID, 1))
            body.append(text(left - 6, y + 3, f"{tick:g}", 9, anchor="end", fill=MUTED))
            x = left + tick * plot_w
            body.append(line(x, top, x, top + plot_h, GRID, 1))
            body.append(text(x, top + plot_h + 14, f"{tick:g}", 9, anchor="middle", fill=MUTED))
        body.append(line(left, top + plot_h, left + plot_w, top + plot_h, INK, 1))
        body.append(line(left, top, left, top + plot_h, INK, 1))
        body.append(line(left, top + plot_h - 0.9 * plot_h, left + plot_w, top + plot_h - 0.9 * plot_h, MUTED, 1, dash="4 4"))
        points = [
            row
            for row in sweep
            if row["family"] == family and float(row["calibrated"]) > 0
        ]
        points.sort(key=lambda row: float(row["coverage"]))
        if len(points) > 1:
            path_points = " ".join(
                f"{left + float(row['coverage']) * plot_w:.1f},{top + plot_h - min(1.0, float(row['precision'])) * plot_h:.1f}"
                for row in points
            )
            body.append(f'<polyline points="{path_points}" fill="none" stroke="{BLUE if selected else GRAY}" stroke-width="2" opacity="0.8"/>')
        for row in points:
            body.append(circle(left + float(row["coverage"]) * plot_w, top + plot_h - min(1.0, float(row["precision"])) * plot_h, 1.2, BLUE if selected else GRAY, 0.35))
        if selected:
            x = left + float(selected["min_coverage"] if "min_coverage" in selected else 0) * plot_w
            # Locate the selected point in the sweep.
            selected_point = next(
                (
                    row
                    for row in points
                    if abs(float(row["min_identity"]) - float(selected["min_identity"])) < 1e-9
                    and abs(float(row["min_coverage"]) - float(selected["min_coverage"])) < 1e-9
                ),
                None,
            )
            if selected_point:
                body.append(circle(left + float(selected_point["coverage"]) * plot_w, top + plot_h - min(1.0, float(selected_point["precision"])) * plot_h, 4.5, BLUE))
                body.append(text(left + float(selected_point["coverage"]) * plot_w + 7, top + plot_h - min(1.0, float(selected_point["precision"])) * plot_h + 4, "selected", 9, fill=BLUE))
        if row == 1:
            body.append(text(left + plot_w / 2, top + plot_h + 32, "coverage", 10, anchor="middle", fill=MUTED))
        if col == 0:
            body.append(text(left - 40, top + plot_h / 2, "precision", 10, anchor="middle", fill=MUTED, rotate=-90))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")



def figure_predictor_auc(rows, path):
    rows = sorted(rows, key=lambda row: -float(row["auc_dbd"]))
    width = 1000
    height = 80 + len(rows) * 26
    body = [text(30, 28, "DBD versus full-length similarity as motif-specificity predictors", 16, weight="600")]
    body.append(text(30, 47, "Rank AUC for PWM similarity >= 0.80; higher is better. Bars are descriptive, not a calibrated model.", 11, fill=MUTED))
    left, right, top = 130, 230, 70
    plot_w = width - left - right
    for tick in (0, 0.25, 0.5, 0.75, 1.0):
        x = left + tick * plot_w
        body.append(line(x, top, x, top + len(rows) * 26, GRID, 1))
        body.append(text(x, top - 8, f"{tick:g}", 9, anchor="middle", fill=MUTED))
    for index, row in enumerate(rows):
        y = top + index * 26
        body.append(text(left - 8, y + 12, row["family"], 10, anchor="end"))
        dbd = float(row["auc_dbd"])
        protein = float(row["auc_protein"])
        body.append(rect(left, y + 2, dbd * plot_w, 9, BLUE, rx=2))
        body.append(rect(left, y + 13, protein * plot_w, 9, ORANGE, rx=2))
        body.append(text(left + max(dbd, protein) * plot_w + 6, y + 11, f"{dbd:.2f} / {protein:.2f}", 9, fill=MUTED))
    body.append(text(left + plot_w / 2, height - 18, "AUC", 10, anchor="middle", fill=MUTED))
    body.append(square(width - 180, 92, 7, BLUE))
    body.append(text(width - 168, 96, "DBD", 10))
    body.append(square(width - 180, 112, 7, ORANGE))
    body.append(text(width - 168, 116, "full-length protein", 10))
    Path(path).write_text(svg_document(width, height, body), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    pairwise = read_tsv(args.dataset / "dbd-motif-pairwise.tsv")
    summary = read_tsv(args.dataset / "calibration-family-summary.tsv")
    sweep = read_tsv(args.dataset / "calibration-threshold-sweep.tsv")
    thresholds = read_tsv(args.dataset / "calibration-thresholds.tsv")
    figure_overall(pairwise, args.out / "dbd-identity-vs-motif-similarity.svg")
    family_order = [
        "ERF",
        "MYB",
        "WRKY",
        "bZIP",
        "NAC",
        "HD-ZIP",
        "TCP",
        "bHLH",
    ]
    figure_families(pairwise, args.out / "family-dbd-motif-similarity.svg", family_order)
    figure_family_sizes(summary, args.out / "family-sample-sizes.svg")
    auc_rows = read_tsv(args.dataset / "pairwise-predictor-auc.tsv")
    figure_predictor_auc(auc_rows, args.out / "predictor-auc.svg")
    figure_precision_coverage(
        sweep,
        thresholds,
        ["HD-ZIP", "TCP", "WRKY", "Dof", "bZIP", "ERF"],
        args.out / "precision-coverage.svg",
    )
    print("wrote 5 SVG figures")


if __name__ == "__main__":
    main()
