#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Plot top TE subtype percentages per species.")
    parser.add_argument("--input", required=True, help="te_subtype_percent.tsv")
    parser.add_argument("--output-dir", required=True, help="Directory for plot files")
    parser.add_argument("--top-n", type=int, default=10, help="Top subtype count per species")
    parser.add_argument("--min-xmax", type=float, default=10.0, help="Minimum x-axis maximum")
    return parser.parse_args()


def read_rows(path):
    by_species = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["percent_of_all_TE"] = float(row["percent_of_all_TE"])
            row["percent_of_genome"] = float(row["percent_of_genome"])
            by_species.setdefault(row["species"], []).append(row)
    return by_species


def fmt_percent(value):
    if value >= 10:
        return f"{value:.1f}"
    if value >= 1:
        return f"{value:.2f}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def bubble_sizes(values):
    if not values:
        return []
    max_value = max(values)
    if max_value <= 0:
        return [20 for _ in values]
    return [35 + 650 * (value / max_value) for value in values]


def plot_species(species, rows, output_dir, top_n, min_xmax):
    rows = sorted(rows, key=lambda row: (-row["percent_of_all_TE"], row["subtype"]))[:top_n]
    labels = [row["subtype"] for row in rows]
    te_values = [row["percent_of_all_TE"] for row in rows]
    genome_values = [row["percent_of_genome"] for row in rows]
    y = list(range(len(rows)))

    height = max(3.8, 0.38 * len(rows) + 1.0)
    fig = plt.figure(figsize=(6.0, height), dpi=300)
    gs = fig.add_gridspec(1, 2, width_ratios=[4.1, 1.35], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    ax_genome = fig.add_subplot(gs[0, 1], sharey=ax)

    ax.barh(y, te_values, color="#625f5c", edgecolor="none", height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    xmax = max(min_xmax, max(te_values) * 1.12 if te_values else min_xmax)
    tick_step = 10 if xmax <= 80 else 20
    xmax = ((int(xmax) + tick_step - 1) // tick_step) * tick_step
    ax.set_xlim(0, xmax)
    ax.set_xlabel("TEs (%)", fontsize=12)
    ax.set_title(species, fontsize=12, pad=8)
    ax.tick_params(axis="x", labelsize=10, width=1.0, length=4)
    ax.tick_params(axis="y", labelsize=10, width=1.0, length=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    sizes = bubble_sizes(genome_values)
    ax_genome.scatter([0.18] * len(y), y, s=sizes, color="black", zorder=3)
    for yi, value in zip(y, genome_values):
        ax_genome.text(0.38, yi, fmt_percent(value), va="center", ha="left", fontsize=10)
    ax_genome.set_xlim(0, 1.0)
    ax_genome.set_title("Genome (%)", fontsize=12, pad=8)
    ax_genome.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    for spine in ax_genome.spines.values():
        spine.set_visible(False)

    fig.tight_layout()
    png = output_dir / f"{species}.top{top_n}_te_subtypes.png"
    pdf = output_dir / f"{species}.top{top_n}_te_subtypes.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    by_species = read_rows(args.input)
    outputs = []
    for species in sorted(by_species):
        outputs.extend(plot_species(species, by_species[species], output_dir, args.top_n, args.min_xmax))
    with open(output_dir / "plot_manifest.tsv", "w", encoding="utf-8") as handle:
        handle.write("file\n")
        for path in outputs:
            handle.write(f"{path}\n")


if __name__ == "__main__":
    main()
