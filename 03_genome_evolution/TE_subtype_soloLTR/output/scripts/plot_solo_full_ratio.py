#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


SPECIAL_OFFSETS = {
    "genome_Mar": (-10, 18),
    "genome_Mru": (-12, 22),
    "genome_ZM4": (-16, -18),
    "genome_482": (12, 18),
    "genome_468": (-14, -18),
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot solo-LTR to full-length LTR ratio against genome size."
    )
    parser.add_argument("--input", required=True, help="solo_full_ltr_summary.tsv")
    parser.add_argument("--output-dir", required=True, help="Directory for plots")
    parser.add_argument(
        "--ratio-column",
        default="solo_to_full_count_ratio",
        choices=["solo_to_full_count_ratio", "solo_to_full_bp_ratio"],
        help="Ratio column to plot",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output filename prefix. Defaults to the ratio column name.",
    )
    parser.add_argument("--zoom-ymax", type=float, default=4.0, help="Y max for linear zoom plot")
    parser.add_argument("--point-color", default="#ef6f6a", help="Point color")
    parser.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="Species IDs to exclude, for example genome_Msa_T2T genome_474_T2T",
    )
    parser.add_argument(
        "--plot-mode",
        choices=["both", "linear", "logy"],
        default="both",
        help="Which plot type to write",
    )
    return parser.parse_args()


def load_rows(path, ratio_column, exclude):
    exclude = set(exclude)
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            species = row["species"]
            if species in exclude:
                continue
            ratio_text = row.get(ratio_column, "")
            if ratio_text in {"", "NA", "nan", "NaN"}:
                continue
            ratio = float(ratio_text)
            if ratio <= 0:
                continue
            label = species
            if label.startswith("genome_"):
                label = label[len("genome_") :]
            rows.append(
                {
                    "species": species,
                    "label": label,
                    "genome_gb": float(row["genome_bp"]) / 1_000_000_000,
                    "ratio": ratio,
                }
            )
    return rows


def x_ticks(rows):
    min_x = min(row["genome_gb"] for row in rows)
    max_x = max(row["genome_gb"] for row in rows)
    ticks = [0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6]
    return [tick for tick in ticks if min_x / 1.4 <= tick <= max_x * 1.8]


def label_offsets(rows):
    offsets = {}
    median_x = sorted(row["genome_gb"] for row in rows)[len(rows) // 2]
    ordered = sorted(rows, key=lambda row: (round(row["genome_gb"], 1), row["ratio"]))
    pattern = [8, -8, 15, -15, 22, -22]
    for index, row in enumerate(ordered):
        dx = 12 if row["genome_gb"] <= median_x else -12
        offsets[row["species"]] = (dx, pattern[index % len(pattern)])
    offsets.update(SPECIAL_OFFSETS)
    return offsets


def style_axes(ax, rows, y_label):
    ax.set_xscale("log", base=2)
    ticks = x_ticks(rows)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{tick:g}" for tick in ticks])
    ax.set_xlim(min(row["genome_gb"] for row in rows) / 1.25, max(row["genome_gb"] for row in rows) * 1.25)
    ax.set_xlabel("Genome size (Gb)", fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.tick_params(axis="both", labelsize=10, width=1.5, length=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)


def add_legend(ax, color):
    handle = Line2D(
        [0],
        [0],
        marker="s",
        color="none",
        markerfacecolor=color,
        markeredgecolor=color,
        markersize=10,
        label="Angiosperms",
    )
    ax.legend(handles=[handle], loc="upper right", bbox_to_anchor=(1.0, 1.18), frameon=False, fontsize=11, handlelength=0.8, handletextpad=0.4)


def annotate_rows(ax, rows, offsets, y_cap=None):
    for row in rows:
        x = row["genome_gb"]
        y = row["ratio"] if y_cap is None else min(row["ratio"], y_cap)
        label = row["label"]
        if y_cap is not None and row["ratio"] > y_cap:
            label = f"{label} ({row['ratio']:.1f})"
        dx, dy = offsets[row["species"]]
        ax.annotate(
            label,
            (x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8,
            fontstyle="italic",
            ha="left" if dx >= 0 else "right",
            va="center",
            arrowprops={"arrowstyle": "-", "lw": 0.7, "color": "black"} if row["ratio"] > (y_cap or math.inf) else None,
            clip_on=False,
        )


def plot_logy(rows, out_path, color, y_label):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.scatter([row["genome_gb"] for row in rows], [row["ratio"] for row in rows], s=38, color=color, edgecolors="white", linewidth=0.8, zorder=3)
    style_axes(ax, rows, y_label)
    ax.set_yscale("log")
    ymin = min(row["ratio"] for row in rows) / 1.5
    ymax = max(row["ratio"] for row in rows) * 1.8
    ax.set_ylim(ymin, ymax)
    ax.text(-0.14, 1.02, "g", transform=ax.transAxes, fontsize=20, fontweight="bold")
    add_legend(ax, color)
    annotate_rows(ax, rows, label_offsets(rows))
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.84)
    fig.savefig(f"{out_path}.png", dpi=300)
    fig.savefig(f"{out_path}.pdf")
    plt.close(fig)


def plot_zoom(rows, out_path, color, y_label, y_max):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    xs = [row["genome_gb"] for row in rows]
    ys = [min(row["ratio"], y_max) for row in rows]
    markers = ["^" if row["ratio"] > y_max else "o" for row in rows]
    for x, y, marker in zip(xs, ys, markers):
        ax.scatter([x], [y], s=42, marker=marker, color=color, edgecolors="white", linewidth=0.8, zorder=3)
    style_axes(ax, rows, y_label)
    ax.set_ylim(0, y_max * 1.12)
    ax.text(-0.14, 1.02, "g", transform=ax.transAxes, fontsize=20, fontweight="bold")
    add_legend(ax, color)
    annotate_rows(ax, rows, label_offsets(rows), y_cap=y_max)
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.84)
    fig.savefig(f"{out_path}.png", dpi=300)
    fig.savefig(f"{out_path}.pdf")
    plt.close(fig)


def main():
    args = parse_args()
    rows = load_rows(args.input, args.ratio_column, args.exclude)
    if not rows:
        raise SystemExit("No rows with a valid ratio were found.")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or args.ratio_column
    if args.ratio_column == "solo_to_full_count_ratio":
        y_label = "Ratio of solo-LTR to fl-LTR"
    else:
        y_label = "Ratio of solo-LTR bp to fl-LTR bp"
    if args.plot_mode in {"both", "logy"}:
        plot_logy(rows, out_dir / f"{prefix}.logy", args.point_color, y_label)
    if args.plot_mode in {"both", "linear"}:
        suffix = "zoom" if args.plot_mode == "both" else "linear"
        plot_zoom(rows, out_dir / f"{prefix}.{suffix}", args.point_color, y_label, args.zoom_ymax)


if __name__ == "__main__":
    main()
