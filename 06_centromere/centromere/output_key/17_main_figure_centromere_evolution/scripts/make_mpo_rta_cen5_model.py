#!/usr/bin/env python3
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.patches import FancyArrowPatch


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})


COL = {
    "chr3": "#7477b8",
    "chr5": "#2f8f83",
    "chr6": "#d99058",
    "relic": "#34a985",
    "active": "#161616",
    "inactive": "#d95f5f",
    "gray": "#333333",
    "light_gray": "#eeeeee",
    "text": "#202020",
}


def save_pub(fig, prefix):
    fig.savefig(prefix + ".svg", bbox_inches="tight")
    fig.savefig(prefix + ".pdf", bbox_inches="tight")
    fig.savefig(prefix + ".png", dpi=450, bbox_inches="tight")
    fig.savefig(prefix + ".tiff", dpi=600, bbox_inches="tight")


def draw_ancestral_module(ax, x, y, w, label, color, cen_label=None, highlight=False):
    ax.plot([x, x + w], [y, y], color=color, lw=10, solid_capstyle="round", zorder=1)
    ax.text(x + w / 2, y + 0.33, label, ha="center", va="bottom",
            fontsize=7.0, fontweight="bold", color=color)
    if cen_label:
        cen_pos = x + w * 0.50
        ax.add_patch(patches.Ellipse((cen_pos, y), 0.24, 0.16,
                                     fc=COL["active"], ec="none", zorder=4))
        if highlight:
            ax.add_patch(patches.Ellipse((cen_pos, y), 0.38, 0.28,
                                         fc="none", ec=COL["inactive"], lw=1.3, zorder=5))
        ax.text(cen_pos, y - 0.31, cen_label, ha="center", va="top",
                fontsize=6.4, color=COL["inactive"] if highlight else COL["text"])


def draw_current_chr(ax, y, name, length_mb, blocks, active_centers, relics):
    x0, x1 = 1.55, 10.85
    ax.plot([x0, x1], [y, y], color=COL["gray"], lw=7, solid_capstyle="round", zorder=0)
    ax.text(0.35, y, name, ha="left", va="center", fontweight="bold")

    for start, end, color, label, label_y in blocks:
        xs = x0 + (start / length_mb) * (x1 - x0)
        xe = x0 + (end / length_mb) * (x1 - x0)
        ax.plot([xs, xe], [y, y], color=color, lw=9, solid_capstyle="butt", zorder=2)
        if label:
            ax.text((xs + xe) / 2, y + label_y, label, ha="center",
                    va="bottom" if label_y > 0 else "top", fontsize=6.2, color=color)

    for pos, label in relics:
        xr = x0 + (pos / length_mb) * (x1 - x0)
        ax.add_patch(patches.Rectangle((xr - 0.10, y - 0.26), 0.20, 0.52,
                                       fc=COL["relic"], ec="white", lw=0.4, zorder=5))
        ax.add_patch(patches.Rectangle((xr - 0.18, y - 0.34), 0.36, 0.68,
                                       fill=False, ec=COL["relic"], lw=1.0, ls="--", zorder=5))
        ax.text(xr, y - 0.43, label, ha="center", va="top", fontsize=6.1, color=COL["relic"])

    for pos, label in active_centers:
        xa = x0 + (pos / length_mb) * (x1 - x0)
        ax.add_patch(patches.Ellipse((xa, y), 0.27, 0.22, fc=COL["active"], ec="none", zorder=6))
        if label:
            ax.text(xa, y + 0.36, label, ha="center", va="bottom", fontsize=5.8)


def curved_arrow(ax, xy1, xy2, color, rad=0.18, lw=1.4, alpha=0.9):
    arrow = FancyArrowPatch(
        xy1, xy2,
        connectionstyle=f"arc3,rad={rad}",
        arrowstyle="-|>",
        mutation_scale=10,
        lw=lw,
        color=color,
        alpha=alpha,
    )
    ax.add_patch(arrow)


def main():
    outdir = Path(__file__).resolve().parents[1] / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.3, 4.15), constrained_layout=True)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5.7)
    ax.axis("off")

    ax.set_title("Mpo 8-to-7 karyotype transition through RTA-like Chr3/Chr5/Chr6 rearrangement",
                 loc="left", fontweight="bold", fontsize=9)

    ax.text(0.42, 5.02, "Ancestral x=8 modules", ha="left", va="center",
            fontsize=7.2, fontweight="bold")
    draw_ancestral_module(ax, 2.1, 4.62, 1.8, "Chr3", COL["chr3"], "CEN3")
    draw_ancestral_module(ax, 5.0, 4.62, 1.8, "Chr5", COL["chr5"], "CEN5", highlight=True)
    draw_ancestral_module(ax, 7.9, 4.62, 1.8, "Chr6", COL["chr6"], "CEN6")

    ax.annotate("RTA-like / complex reciprocal\nrearrangement",
                xy=(6.95, 3.34), xytext=(8.40, 4.02),
                ha="center", va="center", fontsize=7.0,
                arrowprops=dict(arrowstyle="-|>", lw=1.2, color="#555555"))

    ax.text(0.42, 3.25, "Current Mpo x=7 chromosomes", ha="left", va="center",
            fontsize=7.2, fontweight="bold")

    draw_current_chr(
        ax, 2.78, "Mpo Chr3", 62,
        blocks=[
            (18.0, 22.2, COL["chr5"], "Chr5-derived block", 0.34),
            (22.6, 39.5, COL["chr6"], "Chr6-derived block", 0.34),
            (40.4, 50.6, COL["chr3"], "Chr3-derived block", 0.34),
        ],
        active_centers=[(46.6, "")],
        relics=[(21.4, "CEN5 relic")]
    )

    draw_current_chr(
        ax, 1.68, "Mpo Chr5", 60,
        blocks=[
            (25.3, 25.7, COL["chr5"], "", 0.0),
            (26.8, 55.7, COL["chr6"], "Chr6-derived block", 0.34),
            (56.9, 60.0, COL["chr3"], "Chr3-derived block", 0.34),
        ],
        active_centers=[(32.0, "")],
        relics=[(25.45, "CEN5 relic")]
    )

    ax.annotate("", xy=(5.08, 2.62), xytext=(5.90, 4.42),
                arrowprops=dict(arrowstyle="->", color=COL["inactive"], lw=1.1,
                                linestyle="--", shrinkA=2, shrinkB=2))
    ax.annotate("", xy=(5.45, 1.88), xytext=(5.90, 4.42),
                arrowprops=dict(arrowstyle="->", color=COL["inactive"], lw=1.1,
                                linestyle="--", shrinkA=2, shrinkB=2))
    ax.text(6.58, 3.82, "CEN5 fate:\nrelics, not active CENH3",
            ha="center", va="center", fontsize=6.8, color=COL["inactive"])

    ax.add_patch(patches.FancyBboxPatch(
        (0.65, 0.34), 10.55, 0.48,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        fc=COL["light_gray"], ec="none", alpha=0.9, zorder=-1
    ))
    ax.text(
        5.95, 0.58,
        "CEN5-derived relics remain near Chr5-to-Chr6 transition regions, whereas the current active CENH3 domains are elsewhere.",
        ha="center", va="center", fontsize=6.5, color=COL["text"]
    )

    handles = [
        patches.Patch(color=COL["chr5"], label="ancestral Chr5 / CEN5-derived sequence"),
        patches.Patch(color=COL["chr6"], label="ancestral Chr6-derived sequence"),
        patches.Patch(color=COL["chr3"], label="ancestral Chr3-derived sequence"),
        patches.Patch(color=COL["relic"], label="CEN5-derived relic"),
        patches.Ellipse((0, 0), 0.28, 0.18, color=COL["active"], label="current active CENH3"),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.05),
              ncol=3, fontsize=6.2, handlelength=1.3, columnspacing=1.3)

    prefix = str(outdir / "Mpo_RTA_like_CEN5_relic_model_redraw")
    save_pub(fig, prefix)


if __name__ == "__main__":
    main()
