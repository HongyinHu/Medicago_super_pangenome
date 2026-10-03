#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def pvalue_to_stars(pvalue):
    if pvalue is None or np.isnan(pvalue):
        return "NA"
    if pvalue < 0.001:
        return "***"
    if pvalue < 0.01:
        return "**"
    if pvalue < 0.05:
        return "*"
    return "ns"


def mann_whitney(non_cen, cen):
    try:
        from scipy.stats import mannwhitneyu
    except Exception:
        return mann_whitney_normal_approx(non_cen, cen)
    if len(non_cen) == 0 or len(cen) == 0:
        return np.nan
    return mannwhitneyu(non_cen, cen, alternative="two-sided").pvalue


def mann_whitney_normal_approx(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    n1 = len(x)
    n2 = len(y)
    if n1 == 0 or n2 == 0:
        return np.nan
    combined = pd.Series(np.concatenate([x, y]))
    ranks = combined.rank(method="average").to_numpy()
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2
    mean_u = n1 * n2 / 2

    _, counts = np.unique(combined.to_numpy(), return_counts=True)
    tie_term = np.sum(counts**3 - counts)
    n = n1 + n2
    var_u = n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else np.nan
    if not np.isfinite(var_u) or var_u <= 0:
        return np.nan

    z = (abs(u1 - mean_u) - 0.5) / math.sqrt(var_u)
    return math.erfc(z / math.sqrt(2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--out_prefix", required=True)
    parser.add_argument(
        "--genomes",
        nargs="+",
        default=["genome_A17", "genome_R108", "genome_Mpo", "genome_Msa", "genome_474"],
    )
    parser.add_argument("--ylim", type=float, default=6.0)
    args = parser.parse_args()

    root = Path(args.root)
    rows = []
    for genome in args.genomes:
        label = genome.replace("genome_", "")
        path = root / genome / "CENH3_functional_centromere" / f"genome_{label}.LTRRT.cent_noncent.tsv"
        if not path.exists():
            raise FileNotFoundError(path)
        df = pd.read_csv(path, sep="\t")
        df["genome"] = genome
        df["label"] = label
        rows.append(df)
    all_df = pd.concat(rows, ignore_index=True)
    all_df.to_csv(f"{args.out_prefix}.data.tsv", sep="\t", index=False)

    palette = {
        "A17": "#008751",
        "R108": "#ff6a00",
        "Mpo": "#3b6fb6",
        "Msa": "#8a4fb3",
        "474": "#c83f49",
    }

    positions = []
    data = []
    colors = []
    tick_labels = []
    group_centers = []
    stats_rows = []
    pos = 1.0
    for genome in args.genomes:
        label = genome.replace("genome_", "")
        sub = all_df[all_df["label"] == label]
        non = sub[sub["region"] == "Non-cen"]["insertion_mya"].dropna()
        cen = sub[sub["region"] == "Cen"]["insertion_mya"].dropna()
        positions.extend([pos, pos + 0.78])
        data.extend([non, cen])
        colors.extend([palette.get(label, "black"), palette.get(label, "black")])
        tick_labels.extend(["Non-cen", "Cen"])
        group_centers.append((pos + 0.39, label))
        pvalue = mann_whitney(non, cen)
        stats_rows.append(
            {
                "genome": genome,
                "non_cen_n": len(non),
                "cen_n": len(cen),
                "non_cen_median_mya": non.median() if len(non) else np.nan,
                "cen_median_mya": cen.median() if len(cen) else np.nan,
                "mann_whitney_p": pvalue,
                "significance": pvalue_to_stars(pvalue),
            }
        )
        pos += 2.15

    stats = pd.DataFrame(stats_rows)
    stats.to_csv(f"{args.out_prefix}.stats.tsv", sep="\t", index=False)

    fig_width = max(7.5, 1.35 * len(data))
    fig, ax = plt.subplots(figsize=(fig_width, 4.2))
    box = ax.boxplot(
        data,
        positions=positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
    )
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor("white")
        patch.set_edgecolor(color)
        patch.set_linewidth(1.4)
    for i, median in enumerate(box["medians"]):
        median.set_color(colors[i])
        median.set_linewidth(1.4)
    for i, whisker in enumerate(box["whiskers"]):
        whisker.set_color(colors[i // 2])
        whisker.set_linewidth(1.4)
    for i, cap in enumerate(box["caps"]):
        cap.set_color(colors[i // 2])
        cap.set_linewidth(1.4)

    ax.set_ylabel("Insertion time\n(million years ago)", fontsize=13)
    ax.set_ylim(0, args.ylim)
    ax.set_xticks(positions)
    ax.set_xticklabels(tick_labels, rotation=35, ha="right", fontsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for center, label in group_centers:
        ax.text(center, -0.24, label, ha="center", va="top", fontsize=11, transform=ax.get_xaxis_transform())

    legend_handles = []
    for label in [g.replace("genome_", "") for g in args.genomes]:
        color = palette.get(label, "black")
        handle = plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=color, linewidth=1.4)
        legend_handles.append(handle)
    ax.legend(
        legend_handles,
        [g.replace("genome_", "") for g in args.genomes],
        frameon=False,
        loc="upper right",
        ncol=min(5, len(args.genomes)),
        fontsize=11,
    )

    for row, genome in zip(stats_rows, args.genomes):
        label = genome.replace("genome_", "")
        if row["significance"] in {"NA", "ns"}:
            continue
        idx = args.genomes.index(genome)
        x1 = positions[idx * 2]
        x2 = positions[idx * 2 + 1]
        sub = all_df[all_df["label"] == label]
        y = min(args.ylim * 0.92, max(0.35, np.nanpercentile(sub["insertion_mya"], 95) + 0.25))
        ax.plot([x1, x1, x2, x2], [y - 0.08, y, y, y - 0.08], color="black", linewidth=1.0)
        ax.text((x1 + x2) / 2, y + 0.05, row["significance"], ha="center", va="bottom", fontsize=11)

    fig.subplots_adjust(bottom=0.27, left=0.08, right=0.98, top=0.88)
    fig.savefig(f"{args.out_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{args.out_prefix}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"Wrote {args.out_prefix}.pdf")
    print(f"Wrote {args.out_prefix}.png")
    print(f"Wrote {args.out_prefix}.stats.tsv")


if __name__ == "__main__":
    main()
