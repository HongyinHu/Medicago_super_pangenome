#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import numpy as np


def read_chr_pairs(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["sa_chr", "sm_chr"]
    )


def read_links(path):
    return pd.read_csv(path, sep="\t")


def read_bedgraph(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "value"],
        comment="#"
    )


def read_bed3(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end"],
        comment="#"
    )


def get_chr_len_from_bedgraph(bg, chrom):
    sub = bg[bg["chr"] == chrom]
    if len(sub) == 0:
        return 1
    return int(sub["end"].max())


def get_chr_len_from_links(links, chrom, side):
    if side == "sa":
        sub = links[links["sa_chr"] == chrom]
        if len(sub) == 0:
            return 1
        return int(sub["sa_end"].max())
    else:
        sub = links[links["sm_chr"] == chrom]
        if len(sub) == 0:
            return 1
        return int(sub["sm_end"].max())


def scale_pos(x, chr_len, x0, x1):
    return x0 + x / chr_len * (x1 - x0)


def draw_chromosome(ax, x0, x1, y, color, lw=4):
    ax.plot([x0, x1], [y, y], color=color, linewidth=lw, solid_capstyle="butt")


def draw_cenh3_signal(ax, bg, chrom, chr_len, x0, x1, y, height, mode="up"):
    sub = bg[bg["chr"] == chrom].copy()

    if len(sub) == 0:
        return

    vals = sub["value"].replace([np.inf, -np.inf], np.nan).dropna()

    if len(vals) == 0:
        return

    vmax = np.nanpercentile(vals, 99)
    vmin = np.nanpercentile(vals, 1)

    max_abs = max(abs(vmin), abs(vmax))
    if max_abs <= 0:
        max_abs = 1

    xs = []
    ys = []

    for _, row in sub.iterrows():
        mid = (row["start"] + row["end"]) / 2
        x = scale_pos(mid, chr_len, x0, x1)

        value = row["value"]
        if pd.isna(value):
            value = 0

        value = max(min(value, max_abs), -max_abs)

        if mode == "up":
            yy = y + max(value, 0) / max_abs * height
        elif mode == "down":
            yy = y - max(value, 0) / max_abs * height
        else:
            yy = y + value / max_abs * height

        xs.append(x)
        ys.append(yy)

    baseline = y

    if mode == "up":
        ax.fill_between(xs, baseline, ys, color="#123C69", alpha=0.95, linewidth=0)
    elif mode == "down":
        ax.fill_between(xs, baseline, ys, color="#123C69", alpha=0.95, linewidth=0)
    else:
        ax.fill_between(xs, baseline, ys, color="#123C69", alpha=0.95, linewidth=0)

    ax.plot(xs, ys, color="#123C69", linewidth=0.35)


def draw_pericentromere(ax, cen, chrom, chr_len, x0, x1, y, height=0.018):
    sub = cen[cen["chr"] == chrom]

    for _, row in sub.iterrows():
        xs = scale_pos(row["start"], chr_len, x0, x1)
        xe = scale_pos(row["end"], chr_len, x0, x1)

        ax.add_patch(
            Rectangle(
                (xs, y - height / 2),
                xe - xs,
                height,
                facecolor="#E6550D",
                edgecolor="none",
                alpha=0.95,
                zorder=5,
            )
        )


def draw_links(ax, links, sa_chr, sm_chr, sa_len, sm_len, x0, x1, y_sa, y_sm, max_links=None):
    sub = links[
        (links["sa_chr"] == sa_chr) &
        (links["sm_chr"] == sm_chr)
    ].copy()

    if len(sub) == 0:
        return

    if max_links is not None and len(sub) > max_links:
        sub = sub.sample(max_links, random_state=1)

    for _, row in sub.iterrows():
        sa_mid = (row["sa_start"] + row["sa_end"]) / 2
        sm_mid = (row["sm_start"] + row["sm_end"]) / 2

        x_sa = scale_pos(sa_mid, sa_len, x0, x1)
        x_sm = scale_pos(sm_mid, sm_len, x0, x1)

        ax.plot(
            [x_sa, x_sm],
            [y_sa, y_sm],
            color="lightgray",
            linewidth=0.4,
            alpha=0.45,
            zorder=0,
        )


def draw_legend(ax):
    y = 0.95

    ax.text(0.20, y, "Log2(ChIP/input)", fontsize=10, ha="center", va="center")

    ax.plot([0.33, 0.36], [y, y], color="#123C69", linewidth=6, solid_capstyle="butt")
    ax.text(0.375, y, "Sa01", fontsize=10, va="center", ha="left")

    ax.plot([0.47, 0.50], [y, y], color="#A6B7D3", linewidth=6, solid_capstyle="butt")
    ax.text(0.515, y, "Sm01", fontsize=10, va="center", ha="left")

    ax.add_patch(
        Rectangle(
            (0.62, y - 0.01),
            0.035,
            0.02,
            facecolor="#E6550D",
            edgecolor="none"
        )
    )
    ax.text(0.67, y, "Pericentromeric 6 Mb", fontsize=10, va="center", ha="left")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chr_pairs", required=True)
    parser.add_argument("--links", required=True)
    parser.add_argument("--sa_signal", required=True)
    parser.add_argument("--sm_signal", required=True)
    parser.add_argument("--sa_cen", required=True)
    parser.add_argument("--sm_cen", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max_links", type=int, default=None)
    parser.add_argument("--fig_width", type=float, default=10)
    parser.add_argument("--fig_height", type=float, default=5)
    args = parser.parse_args()

    chr_pairs = read_chr_pairs(args.chr_pairs)
    links = read_links(args.links)
    sa_bg = read_bedgraph(args.sa_signal)
    sm_bg = read_bedgraph(args.sm_signal)
    sa_cen = read_bed3(args.sa_cen)
    sm_cen = read_bed3(args.sm_cen)

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

    x0, x1 = 0.12, 0.92

    y_start = 0.82
    row_gap = 0.22

    sa_color = "#123C69"
    sm_color = "#A6B7D3"

    for i, row in chr_pairs.iterrows():
        sa_chr = row["sa_chr"]
        sm_chr = row["sm_chr"]

        y_sa = y_start - i * row_gap
        y_sm = y_sa - 0.075

        sa_len = max(
            get_chr_len_from_bedgraph(sa_bg, sa_chr),
            get_chr_len_from_links(links, sa_chr, "sa")
        )

        sm_len = max(
            get_chr_len_from_bedgraph(sm_bg, sm_chr),
            get_chr_len_from_links(links, sm_chr, "sm")
        )

        ax.text(
            0.065,
            (y_sa + y_sm) / 2,
            sa_chr.replace("Sa", "Chr"),
            fontsize=10,
            ha="right",
            va="center"
        )

        draw_links(
            ax,
            links,
            sa_chr,
            sm_chr,
            sa_len,
            sm_len,
            x0,
            x1,
            y_sa - 0.015,
            y_sm + 0.015,
            max_links=args.max_links
        )

        draw_chromosome(ax, x0, x1, y_sa, sa_color, lw=4)
        draw_chromosome(ax, x0, x1, y_sm, sm_color, lw=4)

        draw_cenh3_signal(
            ax,
            sa_bg,
            sa_chr,
            sa_len,
            x0,
            x1,
            y_sa + 0.005,
            height=0.050,
            mode="up"
        )

        draw_cenh3_signal(
            ax,
            sm_bg,
            sm_chr,
            sm_len,
            x0,
            x1,
            y_sm - 0.005,
            height=0.050,
            mode="down"
        )

        draw_pericentromere(
            ax,
            sa_cen,
            sa_chr,
            sa_len,
            x0,
            x1,
            y_sa,
            height=0.020
        )

        draw_pericentromere(
            ax,
            sm_cen,
            sm_chr,
            sm_len,
            x0,
            x1,
            y_sm,
            height=0.020
        )

    draw_legend(ax)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    plt.savefig(args.out, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()