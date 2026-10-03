#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


GENOME_ORDER = ["genome_A17", "genome_R108", "genome_Mpo", "genome_Msa", "genome_474"]
LABEL_MAP = {
    "genome_A17": "A17",
    "genome_R108": "R108",
    "genome_Mpo": "Mpo",
    "genome_Msa": "Msa",
    "genome_474": "474",
}
COLORS = {
    "genome_A17": "#008751",
    "genome_R108": "#ff6a00",
    "genome_Mpo": "#3b6fb6",
    "genome_Msa": "#8a4fb3",
    "genome_474": "#c83f49",
}
GROUP_RANK = {
    "Copia": 0,
    "Gypsy": 1,
    "LTR": 2,
    "LINE": 3,
    "SINE": 4,
    "MITE": 5,
    "DNA": 6,
    "TIR": 6,
    "Helitron": 6,
    "unknown": 9,
}


def display_label(label):
    return str(label).replace("_unclassified", "\nunc.").replace("_", "\n")


def plot_grouped(summary, out_prefix, suffix="", symlog=False):
    present = [genome for genome in GENOME_ORDER if genome in set(summary["genome"])]
    labels = (
        summary[["plot_group", "fine_type"]]
        .drop_duplicates()
        .assign(rank=lambda d: d["plot_group"].map(GROUP_RANK).fillna(8))
        .sort_values(["rank", "plot_group", "fine_type"])
    )
    label_order = labels["fine_type"].tolist()
    x = np.arange(len(label_order))
    width = min(0.78 / max(len(present), 1), 0.15)
    offsets = (np.arange(len(present)) - (len(present) - 1) / 2) * width

    fig, ax = plt.subplots(figsize=(max(11, len(label_order) * 0.7), 5.2))
    for i, genome in enumerate(present):
        vals = []
        for label in label_order:
            val = summary[(summary["genome"] == genome) & (summary["fine_type"] == label)]["copy_number"].sum()
            vals.append(int(val))
        ax.bar(
            x + offsets[i],
            vals,
            width=width,
            color=COLORS.get(genome, "#777777"),
            label=LABEL_MAP.get(genome, genome),
            edgecolor="white",
            linewidth=0.2,
        )

    ax.set_ylabel("Copy number")
    if symlog:
        ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(x)
    ax.set_xticklabels([display_label(label) for label in label_order], rotation=90)
    ax.text(0.01, 0.98, "Functional centromere", ha="left", va="top", transform=ax.transAxes, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, ncol=min(5, len(present)), loc="upper center", bbox_to_anchor=(0.5, -0.25))

    trans = ax.get_xaxis_transform()
    for group, sublabels in labels.groupby("plot_group", sort=False):
        idxs = [label_order.index(v) for v in sublabels["fine_type"]]
        mid = (min(idxs) + max(idxs)) / 2
        ax.text(mid, -0.37, group, ha="center", va="top", transform=trans, fontsize=10)
        ax.plot(
            [min(idxs) - 0.45, max(idxs) + 0.45],
            [-0.30, -0.30],
            color="black",
            linewidth=0.6,
            transform=trans,
            clip_on=False,
        )

    fig.subplots_adjust(bottom=0.45, left=0.08, right=0.995, top=0.96)
    fig.savefig(f"{out_prefix}.TEsorter_clade_copy_number{suffix}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.TEsorter_clade_copy_number{suffix}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("summaries", nargs="+")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    frames = [pd.read_csv(path, sep="\t") for path in args.summaries]
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(outdir / "all_species.functional_centromere.TE_type_summary.TEsorter_clade.tsv", sep="\t", index=False)

    wide = combined.pivot_table(
        index=["plot_group", "fine_type"],
        columns="genome",
        values="copy_number",
        aggfunc="sum",
        fill_value=0,
    ).reset_index()
    wide.to_csv(outdir / "all_species.functional_centromere.TE_type_copy_number_matrix.TEsorter_clade.tsv", sep="\t", index=False)
    out_prefix = str(outdir / "all_species.functional_centromere.TE_type")
    plot_grouped(combined, out_prefix)
    plot_grouped(combined, out_prefix, suffix=".symlog", symlog=True)
    no_unknown = combined[combined["plot_group"] != "unknown"].copy()
    plot_grouped(no_unknown, out_prefix, suffix=".without_unknown.symlog", symlog=True)
    ltr_major = combined[combined["plot_group"].isin(["Copia", "Gypsy", "LTR"])].copy()
    plot_grouped(ltr_major, out_prefix, suffix=".LTR_Copia_Gypsy_LTRunknown", symlog=False)
    plot_grouped(ltr_major, out_prefix, suffix=".LTR_Copia_Gypsy_LTRunknown.symlog", symlog=True)
    print(f"Combined {len(args.summaries)} summary files")
    print(f"Output directory: {outdir}")


if __name__ == "__main__":
    main()
