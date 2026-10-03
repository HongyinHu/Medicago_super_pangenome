#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


ROOT = Path("path/to/project/10.centromere_analysis/output_key/12_recall_Mpo_CENH3")
BEDGRAPH = ROOT / "bedgraph"
RESULTS = ROOT / "results"
OUTDIR = ROOT / "figures"

COL = {
    "q20": "#2f5f9e",
    "raw": "#9a9a9a",
    "relaxed": "#4d9b5f",
    "cen": "#bf4b55",
    "c5": "#d97941",
    "grey": "#707070",
    "light": "#eef1f4",
    "axis": "#2a2a2a",
}


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.7,
        "legend.frameon": False,
    }
)


def read_sizes(path: Path) -> dict[str, int]:
    out = {}
    with path.open() as fh:
        for line in fh:
            if line.strip():
                chrom, size = line.rstrip("\n").split("\t")[:2]
                if chrom.startswith("Chr"):
                    out[chrom] = int(size)
    return out


def read_consensus(path: Path) -> dict[str, dict]:
    out = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["chrom"]] = row
    return out


def read_bedgraph(path: Path, chrom: str | None = None, start: int | None = None, end: int | None = None) -> dict[str, list[tuple[float, float]]]:
    data: dict[str, list[tuple[float, float]]] = {}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith(("track", "#")):
                continue
            p = raw.rstrip("\n").split()
            c, s, e, v = p[0], int(float(p[1])), int(float(p[2])), float(p[3])
            if not c.startswith("Chr"):
                continue
            if chrom is not None and c != chrom:
                continue
            if start is not None and e <= start:
                continue
            if end is not None and s >= end:
                continue
            mid = (s + e) / 2 / 1e6
            data.setdefault(c, []).append((mid, v))
    return data


def save(fig, stem: str):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.png", dpi=450, bbox_inches="tight")
    plt.close(fig)


def plot_genomewide():
    sizes = read_sizes(RESULTS / "genome_Mpo.genome.sizes")
    consensus = read_consensus(RESULTS / "Mpo_CENH3_recall_consensus.tsv")
    q20 = read_bedgraph(BEDGRAPH / "genome_Mpo.CENH3_vs_Input.q20.10k.smooth50k.CPM.log2.bedGraph")
    relaxed = read_bedgraph(BEDGRAPH / "genome_Mpo.CENH3_vs_Input.relaxed.10k.smooth50k.CPM.log2.bedGraph")

    chroms = [f"Chr{i}" for i in range(1, 8)]
    fig, axes = plt.subplots(len(chroms), 1, figsize=(7.2, 7.0), sharey=True)
    fig.suptitle("Recalled Mpo CENH3/Input enrichment across seven chromosomes", y=0.995, fontsize=10, fontweight="bold")

    for ax, chrom in zip(axes, chroms):
        for mode, data, color, alpha, lw in [
            ("relaxed", relaxed, COL["relaxed"], 0.35, 0.55),
            ("q20", q20, COL["q20"], 0.95, 0.75),
        ]:
            xy = data.get(chrom, [])
            if xy:
                x = [a for a, _ in xy]
                y = [b for _, b in xy]
                ax.plot(x, y, color=color, lw=lw, alpha=alpha, label=mode if ax is axes[0] else None)
        ax.axhline(0, color="#bfc5cc", lw=0.45)
        ax.axhline(0.5, color="#c9a9ae", lw=0.45, ls=(0, (3, 2)))
        if chrom in consensus:
            r = consensus[chrom]
            s, e = int(r["consensus_start"]) / 1e6, int(r["consensus_end"]) / 1e6
            ax.axvspan(s, e, color=COL["cen"], alpha=0.22, lw=0)
            ax.text((s + e) / 2, 2.40, "recalled CENH3", ha="center", va="top", fontsize=6.2, color=COL["cen"])
        ax.set_xlim(0, sizes[chrom] / 1e6)
        ax.set_ylim(-1.2, 2.6)
        ax.set_ylabel(chrom, rotation=0, ha="right", va="center", fontsize=7.5)
        ax.tick_params(axis="both", labelsize=6.5, length=2)
        ax.spines["left"].set_color("#aeb5bd")
        ax.spines["bottom"].set_color("#aeb5bd")
    axes[-1].set_xlabel("Position on current Mpo assembly (Mb)", fontsize=8)
    axes[0].legend(loc="upper right", ncol=2, fontsize=7, handlelength=1.4)
    fig.text(0.01, 0.50, "log2(CENH3/Input), 10-kb bins, 50-kb smoothing", rotation=90, va="center", fontsize=8)
    save(fig, "Mpo_CENH3_recall_signal_genomewide")


def plot_chr4_focus():
    consensus = read_consensus(RESULTS / "Mpo_CENH3_recall_consensus.tsv")
    chrom = "Chr4"
    start, end = 16_000_000, 26_500_000
    c5_start, c5_end = 18_316_439, 18_432_188
    data = {
        "raw": read_bedgraph(BEDGRAPH / "genome_Mpo.CENH3_vs_Input.raw.10k.smooth50k.CPM.log2.bedGraph", chrom, start, end).get(chrom, []),
        "relaxed": read_bedgraph(BEDGRAPH / "genome_Mpo.CENH3_vs_Input.relaxed.10k.smooth50k.CPM.log2.bedGraph", chrom, start, end).get(chrom, []),
        "q20": read_bedgraph(BEDGRAPH / "genome_Mpo.CENH3_vs_Input.q20.10k.smooth50k.CPM.log2.bedGraph", chrom, start, end).get(chrom, []),
    }

    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 0.9], hspace=0.14)
    ax = fig.add_subplot(gs[0, 0])
    for mode, color, lw, alpha in [("raw", COL["raw"], 0.55, 0.55), ("relaxed", COL["relaxed"], 0.65, 0.75), ("q20", COL["q20"], 0.9, 0.95)]:
        xy = data[mode]
        ax.plot([x for x, _ in xy], [y for _, y in xy], color=color, lw=lw, alpha=alpha, label=mode)

    ax.axhline(0, color="#bfc5cc", lw=0.5)
    ax.axhline(0.5, color="#c9a9ae", lw=0.5, ls=(0, (3, 2)))
    cen = consensus[chrom]
    cen_s, cen_e = int(cen["consensus_start"]), int(cen["consensus_end"])
    ax.axvspan(cen_s / 1e6, cen_e / 1e6, color=COL["cen"], alpha=0.22, lw=0)
    ax.axvspan(c5_start / 1e6, c5_end / 1e6, color=COL["c5"], alpha=0.28, lw=0)
    ax.text((c5_start + c5_end) / 2 / 1e6, 2.55, "CEN5-like\nprojection", ha="center", va="top", fontsize=7, color=COL["c5"])
    ax.text((cen_s + cen_e) / 2 / 1e6, 2.55, "recalled primary\nCENH3 domain", ha="center", va="top", fontsize=7, color=COL["cen"])
    ax.set_title("Mpo Chr4 CENH3 signal: CEN5-like projection is separate from the recalled primary CENH3 domain", loc="left", fontsize=9, fontweight="bold")
    ax.set_xlim(start / 1e6, end / 1e6)
    ax.set_ylim(-1.2, 2.8)
    ax.set_ylabel("log2(CENH3/Input)", fontsize=8)
    ax.tick_params(labelsize=7, length=2)
    ax.legend(loc="upper right", ncol=3, fontsize=7, handlelength=1.3)
    ax.spines["left"].set_color("#aeb5bd")
    ax.spines["bottom"].set_color("#aeb5bd")

    ax2 = fig.add_subplot(gs[1, 0], sharex=ax)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.hlines(0.48, start / 1e6, end / 1e6, color="#c5cbd0", lw=6)
    ax2.add_patch(Rectangle((c5_start / 1e6, 0.27), (c5_end - c5_start) / 1e6, 0.42, facecolor=COL["c5"], edgecolor="none"))
    ax2.add_patch(Rectangle((cen_s / 1e6, 0.27), (cen_e - cen_s) / 1e6, 0.42, facecolor=COL["cen"], edgecolor="none", alpha=0.85))
    ax2.text((c5_start + c5_end) / 2 / 1e6, 0.05, "18.32-18.43 Mb", ha="center", va="bottom", fontsize=6.8, color=COL["c5"])
    ax2.text((cen_s + cen_e) / 2 / 1e6, 0.05, "21.96-24.57 Mb", ha="center", va="bottom", fontsize=6.8, color=COL["cen"])
    ax2.set_xlabel("Position on Mpo Chr4 (Mb)", fontsize=8)
    ax2.tick_params(axis="x", labelsize=7, length=2)
    ax2.spines["left"].set_visible(False)
    ax2.spines["bottom"].set_color("#aeb5bd")
    save(fig, "Mpo_CENH3_recall_signal_Chr4_focus")


def main():
    plot_genomewide()
    plot_chr4_focus()


if __name__ == "__main__":
    main()
