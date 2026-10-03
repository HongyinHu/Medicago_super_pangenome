#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXCLUDE = {"Unknown", "Unclassified LTR-RTs"}


def parse_args():
    parser = argparse.ArgumentParser(description="Plot classified LTR-RT subtype percentages per species.")
    parser.add_argument("--input", required=True, help="ltr_rt_subtype_percent.tsv")
    parser.add_argument("--output-dir", required=True, help="Directory for plot files")
    parser.add_argument("--top-n", type=int, default=10, help="Top classified LTR-RT subtype count per species")
    return parser.parse_args()


def read_rows(path):
    by_species = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["subtype"] in EXCLUDE:
                continue
            row["percent_of_ltr_rt"] = float(row["percent_of_ltr_rt"])
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


def nice_xmax(values):
    xmax = max(values) * 1.12 if values else 10.0
    xmax = max(10.0, xmax)
    step = 10 if xmax <= 80 else 20
    return ((int(xmax) + step - 1) // step) * step


def plot_species(species, rows, output_dir, top_n):
    rows = sorted(rows, key=lambda row: (-row["percent_of_ltr_rt"], row["subtype"]))[:top_n]
    labels = [row["subtype"] for row in rows]
    ltr_values = [row["percent_of_ltr_rt"] for row in rows]
    genome_values = [row["percent_of_genome"] for row in rows]
    y = list(range(len(rows)))

    height = max(3.8, 0.38 * len(rows) + 1.0)
    fig = plt.figure(figsize=(6.0, height), dpi=300)
    gs = fig.add_gridspec(1, 2, width_ratios=[4.1, 1.35], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    ax_genome = fig.add_subplot(gs[0, 1], sharey=ax)

    ax.barh(y, ltr_values, color="#625f5c", edgecolor="none", height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, nice_xmax(ltr_values))
    ax.set_xlabel("LTR-RTs (%)", fontsize=12)
    ax.set_title(species, fontsize=12, pad=8)
    ax.tick_params(axis="x", labelsize=10, width=1.0, length=4)
    ax.tick_params(axis="y", labelsize=10, width=1.0, length=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    ax_genome.scatter([0.18] * len(y), y, s=bubble_sizes(genome_values), color="black", zorder=3)
    for yi, value in zip(y, genome_values):
        ax_genome.text(0.38, yi, fmt_percent(value), va="center", ha="left", fontsize=10)
    ax_genome.set_xlim(0, 1.0)
    ax_genome.set_title("Genome (%)", fontsize=12, pad=8)
    ax_genome.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    for spine in ax_genome.spines.values():
        spine.set_visible(False)

    png = output_dir / f"{species}.classified_ltr_rt_top{top_n}.png"
    pdf = output_dir / f"{species}.classified_ltr_rt_top{top_n}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def write_filtered_table(by_species, output_dir, top_n):
    path = output_dir / f"classified_ltr_rt_top{top_n}.tsv"
    with open(path, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "species",
            "subtype",
            "bp_merged",
            "annotation_count",
            "percent_of_genome",
            "percent_of_ltr_rt",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for species in sorted(by_species):
            rows = sorted(by_species[species], key=lambda row: (-row["percent_of_ltr_rt"], row["subtype"]))[:top_n]
            for row in rows:
                writer.writerow({key: row[key] for key in fieldnames})
    return path


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    by_species = read_rows(args.input)

    outputs = []
    for species in sorted(by_species):
        outputs.extend(plot_species(species, by_species[species], output_dir, args.top_n))
    filtered_table = write_filtered_table(by_species, output_dir, args.top_n)
    with open(output_dir / "plot_manifest.tsv", "w", encoding="utf-8") as handle:
        handle.write("file\n")
        for path in outputs:
            handle.write(f"{path}\n")
        handle.write(f"{filtered_table}\n")


if __name__ == "__main__":
    main()
