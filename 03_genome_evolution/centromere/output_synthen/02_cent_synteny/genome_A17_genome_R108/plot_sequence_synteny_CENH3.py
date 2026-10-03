#!/usr/bin/env python3
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def read_chr_pairs(path):
    return pd.read_csv(path, sep="\t", header=None, names=["query_chr", "target_chr"])


def read_bedgraph(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "value"],
        usecols=[0, 1, 2, 3],
        comment="#",
    )


def read_bed3(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end"],
        usecols=[0, 1, 2],
        comment="#",
    )


def scale_pos(x, chrom_len, x0, x1):
    if chrom_len <= 0:
        chrom_len = 1
    return x0 + x / chrom_len * (x1 - x0)


def chr_len(bg, blocks, chrom, side):
    vals = []
    sub_bg = bg[bg["chr"] == chrom]
    if len(sub_bg):
        vals.append(int(sub_bg["end"].max()))
    if side == "query":
        sub = blocks[blocks["query_chr"] == chrom]
        if len(sub):
            vals.append(int(sub["query_end"].max()))
    else:
        sub = blocks[blocks["target_chr"] == chrom]
        if len(sub):
            vals.append(int(sub["target_end"].max()))
    return max(vals) if vals else 1


def draw_chromosome(ax, x0, x1, y, color, lw=4):
    ax.plot([x0, x1], [y, y], color=color, linewidth=lw, solid_capstyle="butt")


def draw_signal(ax, bg, chrom, chrom_len, x0, x1, y, height, color, mode):
    sub = bg[bg["chr"] == chrom].copy()
    if sub.empty:
        return
    vals = sub["value"].replace([np.inf, -np.inf], np.nan).dropna()
    if vals.empty:
        return
    vmax = np.nanpercentile(vals, 99)
    vmin = np.nanpercentile(vals, 1)
    max_abs = max(abs(vmin), abs(vmax), 1e-9)

    xs = []
    ys = []
    for _, row in sub.iterrows():
        mid = (row["start"] + row["end"]) / 2
        x = scale_pos(mid, chrom_len, x0, x1)
        value = row["value"]
        if pd.isna(value):
            value = 0
        value = max(min(value, max_abs), -max_abs)
        if mode == "up":
            yy = y + max(value, 0) / max_abs * height
        else:
            yy = y - max(value, 0) / max_abs * height
        xs.append(x)
        ys.append(yy)

    ax.fill_between(xs, y, ys, color=color, alpha=0.95, linewidth=0)
    ax.plot(xs, ys, color=color, linewidth=0.35)


def draw_regions(ax, bed, chrom, chrom_len, x0, x1, y, height=0.018):
    sub = bed[bed["chr"] == chrom]
    for _, row in sub.iterrows():
        xs = scale_pos(row["start"], chrom_len, x0, x1)
        xe = scale_pos(row["end"], chrom_len, x0, x1)
        ax.add_patch(
            Rectangle(
                (xs, y - height / 2),
                xe - xs,
                height,
                facecolor="#E6550D",
                edgecolor="none",
                alpha=0.95,
                zorder=6,
            )
        )


def draw_links(ax, blocks, query_chr, target_chr, query_len, target_len, x0, x1, y_query, y_target, max_links):
    sub = blocks[(blocks["query_chr"] == query_chr) & (blocks["target_chr"] == target_chr)].copy()
    if sub.empty:
        return 0
    if max_links and len(sub) > max_links:
        sub = sub.sample(max_links, random_state=1)
    for _, row in sub.iterrows():
        q_mid = (row["query_start"] + row["query_end"]) / 2
        t_mid = (row["target_start"] + row["target_end"]) / 2
        xq = scale_pos(q_mid, query_len, x0, x1)
        xt = scale_pos(t_mid, target_len, x0, x1)
        ax.plot([xq, xt], [y_query, y_target], color="lightgray", linewidth=0.35, alpha=0.42, zorder=0)
    return len(sub)


def draw_legend(ax, query_label, target_label):
    y = 0.95
    ax.text(0.18, y, "Log2(ChIP/input)", fontsize=10, ha="center", va="center")
    ax.plot([0.33, 0.36], [y, y], color="#123C69", linewidth=6, solid_capstyle="butt")
    ax.text(0.375, y, query_label, fontsize=10, va="center", ha="left")
    ax.plot([0.47, 0.50], [y, y], color="#A6B7D3", linewidth=6, solid_capstyle="butt")
    ax.text(0.515, y, target_label, fontsize=10, va="center", ha="left")
    ax.add_patch(Rectangle((0.62, y - 0.01), 0.035, 0.02, facecolor="#E6550D", edgecolor="none"))
    ax.text(0.67, y, "Pericentromeric 6 Mb", fontsize=10, va="center", ha="left")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks", required=True)
    parser.add_argument("--chr_pairs", required=True)
    parser.add_argument("--query_signal", required=True)
    parser.add_argument("--target_signal", required=True)
    parser.add_argument("--query_cen", required=True)
    parser.add_argument("--target_cen", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min_len", type=int, default=10000)
    parser.add_argument("--max_links", type=int, default=1000)
    parser.add_argument("--fig_width", type=float, default=10)
    parser.add_argument("--fig_height", type=float, default=5)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--query_label", default="A17")
    parser.add_argument("--target_label", default="R108")
    args = parser.parse_args()

    blocks = pd.read_csv(args.blocks, sep="\t")
    blocks = blocks[(blocks["query_len"] >= args.min_len) & (blocks["target_len"] >= args.min_len)].copy()
    pairs = read_chr_pairs(args.chr_pairs)
    query_bg = read_bedgraph(args.query_signal)
    target_bg = read_bedgraph(args.target_signal)
    query_cen = read_bed3(args.query_cen)
    target_cen = read_bed3(args.target_cen)

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
    x0, x1 = 0.12, 0.92
    y_start = 0.82
    row_gap = 0.22 if len(pairs) <= 4 else 0.17
    q_color = "#123C69"
    t_color = "#A6B7D3"

    for i, row in pairs.iterrows():
        q_chr = row["query_chr"]
        t_chr = row["target_chr"]
        yq = y_start - i * row_gap
        yt = yq - 0.075
        q_len = chr_len(query_bg, blocks, q_chr, "query")
        t_len = chr_len(target_bg, blocks, t_chr, "target")
        ax.text(0.065, (yq + yt) / 2, q_chr.replace("Chr", "Chr"), fontsize=10, ha="right", va="center")
        draw_links(ax, blocks, q_chr, t_chr, q_len, t_len, x0, x1, yq - 0.015, yt + 0.015, args.max_links)
        draw_chromosome(ax, x0, x1, yq, q_color)
        draw_chromosome(ax, x0, x1, yt, t_color)
        draw_signal(ax, query_bg, q_chr, q_len, x0, x1, yq + 0.005, 0.050, q_color, "up")
        draw_signal(ax, target_bg, t_chr, t_len, x0, x1, yt - 0.005, 0.050, q_color, "down")
        draw_regions(ax, query_cen, q_chr, q_len, x0, x1, yq)
        draw_regions(ax, target_cen, t_chr, t_len, x0, x1, yt)

    draw_legend(ax, args.query_label, args.target_label)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    plt.savefig(args.out, bbox_inches="tight", dpi=args.dpi)
    plt.close()


if __name__ == "__main__":
    main()
