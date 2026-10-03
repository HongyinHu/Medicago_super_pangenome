#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt


CHR_RE = re.compile(r"^Chr(\d+)$")


COL = {
    "q20": "#2f5f9e",
    "raw": "#9a9a9a",
    "relaxed": "#4d9b5f",
    "cen": "#bf4b55",
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


def chr_key(chrom: str) -> tuple[int, str]:
    m = CHR_RE.match(chrom)
    return (int(m.group(1)), chrom) if m else (10**9, chrom)


def read_sizes(path: Path) -> dict[str, int]:
    out = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, size = line.rstrip("\n").split("\t")[:2]
            if CHR_RE.match(chrom):
                out[chrom] = int(float(size))
    return dict(sorted(out.items(), key=lambda kv: chr_key(kv[0])))


def read_consensus(path: Path) -> dict[str, dict]:
    out = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["chrom"]] = row
    return out


def read_bedgraph(path: Path, keep_chroms: set[str]) -> dict[str, list[tuple[float, float]]]:
    data: dict[str, list[tuple[float, float]]] = {}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith(("track", "#")):
                continue
            p = raw.rstrip("\n").split()
            chrom, start, end, value = p[0], int(float(p[1])), int(float(p[2])), float(p[3])
            if chrom not in keep_chroms:
                continue
            mid = (start + end) / 2 / 1e6
            data.setdefault(chrom, []).append((mid, value))
    return data


def save(fig, outdir: Path, stem: str):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=450, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--root", required=True)
    args = parser.parse_args()

    species = args.species
    root = Path(args.root)
    bedgraph = root / species / "bedgraph"
    results = root / species / "results"
    outdir = root / species / "figures"
    sizes = read_sizes(results / f"{species}.genome.sizes")
    chroms = list(sizes)
    keep = set(chroms)
    curated_path = results / f"{species}_CENH3_recall_curated_consensus.tsv"
    if curated_path.exists():
        raw_consensus = read_consensus(curated_path)
        consensus = {
            chrom: {
                "consensus_start": row["curated_start"],
                "consensus_end": row["curated_end"],
                "q20_status": row.get("q20_status", ""),
            }
            for chrom, row in raw_consensus.items()
        }
        consensus_label = "high-confidence CENH3"
    else:
        consensus = read_consensus(results / f"{species}_CENH3_recall_consensus.tsv")
        consensus_label = "recalled CENH3"
    q20 = read_bedgraph(bedgraph / f"{species}.CENH3_vs_Input.q20.10k.smooth50k.CPM.log2.bedGraph", keep)
    relaxed = read_bedgraph(bedgraph / f"{species}.CENH3_vs_Input.relaxed.10k.smooth50k.CPM.log2.bedGraph", keep)
    raw = read_bedgraph(bedgraph / f"{species}.CENH3_vs_Input.raw.10k.smooth50k.CPM.log2.bedGraph", keep)

    fig_h = max(6.4, 0.86 * len(chroms))
    fig, axes = plt.subplots(len(chroms), 1, figsize=(7.2, fig_h), sharey=True)
    if len(chroms) == 1:
        axes = [axes]
    fig.suptitle(f"Recalled {species} CENH3/Input enrichment across main chromosomes", y=0.995, fontsize=10, fontweight="bold")

    for ax, chrom in zip(axes, chroms):
        for mode, data, color, alpha, lw in [
            ("raw", raw, COL["raw"], 0.26, 0.45),
            ("relaxed", relaxed, COL["relaxed"], 0.48, 0.55),
            ("q20", q20, COL["q20"], 0.58, 0.62),
        ]:
            xy = data.get(chrom, [])
            if xy:
                ax.plot([a for a, _ in xy], [b for _, b in xy], color=color, lw=lw, alpha=alpha, label=mode if ax is axes[0] else None)
        ax.axhline(0, color="#bfc5cc", lw=0.45)
        ax.axhline(0.5, color="#c9a9ae", lw=0.45, ls=(0, (3, 2)))
        if chrom in consensus:
            r = consensus[chrom]
            s, e = int(r["consensus_start"]) / 1e6, int(r["consensus_end"]) / 1e6
            ax.axvspan(s, e, color=COL["cen"], alpha=0.22, lw=0)
            ax.text(
                (s + e) / 2,
                2.45,
                f"{consensus_label}\n{s:.2f}-{e:.2f} Mb",
                ha="center",
                va="top",
                fontsize=6.0,
                color=COL["cen"],
                linespacing=0.95,
            )
        ax.set_xlim(0, sizes[chrom] / 1e6)
        ax.set_ylim(-1.2, 2.6)
        ax.set_ylabel(chrom, rotation=0, ha="right", va="center", fontsize=7.5)
        ax.tick_params(axis="both", labelsize=6.5, length=2)
        ax.spines["left"].set_color("#aeb5bd")
        ax.spines["bottom"].set_color("#aeb5bd")
    axes[-1].set_xlabel(f"Position on current {species} assembly (Mb)", fontsize=8)
    axes[0].legend(loc="upper right", ncol=3, fontsize=7, handlelength=1.4)
    fig.text(0.01, 0.50, "log2(CENH3/Input), 10-kb bins, 50-kb smoothing", rotation=90, va="center", fontsize=8)
    save(fig, outdir, f"{species}_CENH3_recall_signal_genomewide")


if __name__ == "__main__":
    main()
