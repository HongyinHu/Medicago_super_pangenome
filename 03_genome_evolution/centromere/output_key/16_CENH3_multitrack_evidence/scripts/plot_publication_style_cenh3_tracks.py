#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle


SPECIES = ["genome_Msa", "genome_474"]
BIN = 100_000
LOCAL_FLANK = 10_000_000

COL = {
    "raw": "#9a9a9a",
    "relaxed": "#4d9b5f",
    "q20": "#2f5f9e",
    "cen": "#bf4b55",
    "sat": "#7463a9",
    "repeat": "#c9c9c9",
    "gene": "#d95f59",
    "A": "#c64b3b",
    "B": "#3c7da0",
    "tad": "#606060",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 6.5,
        "axes.linewidth": 0.55,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "legend.frameon": False,
    }
)


def read_bed(path: Path) -> list[dict]:
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            rows.append({"chrom": p[0], "start": int(float(p[1])), "end": int(float(p[2])), "name": p[3] if len(p) > 3 else ""})
    return rows


def read_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            c, s = line.rstrip("\n").split("\t")[:2]
            if c.startswith("Chr"):
                sizes[c] = int(float(s))
    return dict(sorted(sizes.items(), key=lambda kv: int(kv[0].replace("Chr", ""))))


def read_bedgraph(path: Path, chroms: set[str]) -> dict[str, list[tuple[int, int, float]]]:
    out = {c: [] for c in chroms}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith(("track", "#")):
                continue
            p = raw.rstrip("\n").split()
            if len(p) < 4 or p[0] not in chroms:
                continue
            try:
                value = float(p[3])
            except ValueError:
                continue
            if math.isfinite(value):
                out[p[0]].append((int(float(p[1])), int(float(p[2])), value))
    return out


def read_intervals_from_gff(path: Path, chroms: set[str], feature: str | None = None) -> dict[str, list[tuple[int, int]]]:
    out = {c: [] for c in chroms}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] not in chroms:
                continue
            if feature is not None and p[2] != feature:
                continue
            out[p[0]].append((int(p[3]) - 1, int(p[4])))
    return out


def read_intervals_from_bed(path: Path, chroms: set[str]) -> dict[str, list[tuple[int, int]]]:
    out = {c: [] for c in chroms}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) >= 3 and p[0] in chroms:
                out[p[0]].append((int(float(p[1])), int(float(p[2]))))
    return out


def read_ab(path: Path, chroms: set[str]) -> dict[str, list[tuple[int, int, float]]]:
    out = {c: [] for c in chroms}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            chrom = row.get("chrom", "")
            if chrom not in chroms or not row.get("E1"):
                continue
            try:
                out[chrom].append((int(row["start"]), int(row["end"]), float(row["E1"])))
            except ValueError:
                pass
    return out


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    intervals = sorted(intervals)
    if not intervals:
        return []
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def ov(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def coverage_track(intervals: list[tuple[int, int]], start: int, end: int, bin_size: int = BIN) -> tuple[list[float], list[float]]:
    n = max(1, math.ceil((end - start) / bin_size))
    cov = np.zeros(n, dtype=float)
    lengths = np.array([min(start + (i + 1) * bin_size, end) - (start + i * bin_size) for i in range(n)], dtype=float)
    for s, e in merge_intervals([(max(start, s), min(end, e)) for s, e in intervals if e > start and s < end]):
        b0 = max(0, (s - start) // bin_size)
        b1 = min(n - 1, (max(s, e - 1) - start) // bin_size)
        for b in range(int(b0), int(b1) + 1):
            bs = start + b * bin_size
            be = min(bs + bin_size, end)
            cov[b] += ov((bs, be), (s, e))
    xs = [(start + i * bin_size + min(bin_size, end - (start + i * bin_size)) / 2) / 1e6 for i in range(n)]
    ys = (cov / np.maximum(lengths, 1)).tolist()
    return xs, ys


def count_track(intervals: list[tuple[int, int]], start: int, end: int, bin_size: int = BIN) -> tuple[list[float], list[float]]:
    n = max(1, math.ceil((end - start) / bin_size))
    counts = np.zeros(n, dtype=float)
    for s, e in intervals:
        if e <= start or s >= end:
            continue
        ss, ee = max(start, s), min(end, e)
        b0 = max(0, (ss - start) // bin_size)
        b1 = min(n - 1, (max(ss, ee - 1) - start) // bin_size)
        counts[int(b0) : int(b1) + 1] += 1
    xs = [(start + i * bin_size + min(bin_size, end - (start + i * bin_size)) / 2) / 1e6 for i in range(n)]
    return xs, counts.tolist()


def line_track(rows: list[tuple[int, int, float]], start: int, end: int) -> tuple[list[float], list[float]]:
    xs, ys = [], []
    for s, e, v in rows:
        if e <= start or s >= end:
            continue
        xs.append((s + e) / 2 / 1e6)
        ys.append(v)
    return xs, ys


def add_common(ax, cen: dict, xlim: tuple[float, float], show_cen=True):
    ax.set_xlim(*xlim)
    if show_cen:
        ax.axvspan(cen["start"] / 1e6, cen["end"] / 1e6, color=COL["cen"], alpha=0.15, lw=0)
    ax.spines["left"].set_color("#aeb5bd")
    ax.spines["bottom"].set_color("#aeb5bd")
    ax.tick_params(axis="both", labelsize=5.6, length=2, pad=1)


def plot_region_tracks(axes, chrom, start, end, cen, data, include_hic=False):
    xlim = (start / 1e6, end / 1e6)
    for ax in axes:
        add_common(ax, cen, xlim)
    for ax in axes[:-1]:
        ax.set_xticklabels([])

    ax = axes[0]
    for mode, alpha, lw in [("raw", 0.26, 0.45), ("relaxed", 0.48, 0.55), ("q20", 0.58, 0.62)]:
        xs, ys = line_track(data["cenh3"][mode].get(chrom, []), start, end)
        ax.plot(xs, ys, color=COL[mode], alpha=alpha, lw=lw, label=mode)
    ax.axhline(0, color="#c8cdd2", lw=0.4)
    ax.axhline(0.5, color="#d8a2a8", lw=0.45, ls=(0, (3, 2)))
    ax.set_ylim(-1.1, 2.6)
    ax.set_ylabel("CENH3", rotation=0, ha="right", va="center")

    ax = axes[1]
    sx, sy = coverage_track(data["sat"].get(chrom, []), start, end)
    ax.fill_between(sx, sy, color=COL["sat"], alpha=0.78, step="mid")
    ax.set_ylim(0, 1)
    ax.set_ylabel("sat", rotation=0, ha="right", va="center")

    ax = axes[2]
    rx, ry = coverage_track(data["repeat"].get(chrom, []), start, end)
    ax.fill_between(rx, ry, color=COL["repeat"], alpha=0.92, step="mid")
    ax.set_ylim(0, 1)
    ax.set_ylabel("repeat", rotation=0, ha="right", va="center")

    ax = axes[3]
    gx, gy = count_track(data["gene"].get(chrom, []), start, end)
    ax.plot(gx, gy, color=COL["gene"], lw=0.68, alpha=0.9)
    ax.set_ylim(0, max(2, max(gy + [1]) * 1.15))
    ax.set_ylabel("gene", rotation=0, ha="right", va="center")

    ax = axes[4]
    ab_vals = [(s, e, v) for s, e, v in data["ab"].get(chrom, []) if e > start and s < end]
    for s, e, v in ab_vals:
        ax.add_patch(
            Rectangle(
                (max(s, start) / 1e6, min(0, v)),
                (min(e, end) - max(s, start)) / 1e6,
                abs(v),
                facecolor=COL["A"] if v >= 0 else COL["B"],
                edgecolor="none",
                alpha=0.82,
            )
        )
    ax.axhline(0, color="#c8cdd2", lw=0.45)
    m = max([abs(v) for _, _, v in ab_vals] + [0.2])
    ax.set_ylim(-m * 1.08, m * 1.08)
    ax.set_ylabel("A/B", rotation=0, ha="right", va="center")

    ax = axes[5]
    tx, ty = line_track(data["tad_score"].get(chrom, []), start, end)
    ax.plot(tx, ty, color=COL["tad"], lw=0.55, alpha=0.85)
    for s, e in data["tad_boundary"].get(chrom, []):
        if e > start and s < end:
            ax.axvline((s + e) / 2 / 1e6, color="#222222", lw=0.25, alpha=0.55)
    if ty:
        ymin, ymax = min(ty), max(ty)
        if ymin == ymax:
            ymin, ymax = ymin - 1, ymax + 1
        ax.set_ylim(ymin - 0.08 * abs(ymax - ymin), ymax + 0.08 * abs(ymax - ymin))
    ax.set_ylabel("TAD", rotation=0, ha="right", va="center")


def hic_triangle_array(cool_path: Path, chrom: str, start: int, end: int, max_distance_bins: int = 100):
    try:
        import cooler
    except Exception:
        return None
    try:
        clr = cooler.Cooler(str(cool_path))
        mat = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    except Exception:
        return None
    mat = np.asarray(mat, dtype=float)
    mat[~np.isfinite(mat)] = 0
    mat = np.log10(mat + 1)
    n = mat.shape[0]
    maxd = min(max_distance_bins, n)
    tri = np.full((maxd, n), np.nan)
    for d in range(maxd):
        vals = np.diag(mat, k=d)
        if len(vals):
            tri[d, : len(vals)] = vals
    return tri


def load_data(project: Path, species: str):
    recall = project / "output_key" / "15_recall_Msa_474_CENH3" / species
    hic = project / "output_all" / "08_hic_AB_TAD" / species
    repeat_root = project / "output_all" / "03_repeat"
    cen = read_bed(recall / "results" / f"{species}.CENH3.functional_centromere.recall.high_confidence.bed")
    cen = sorted(cen, key=lambda r: int(r["chrom"].replace("Chr", "")))
    chroms = {x["chrom"] for x in cen}
    data = {
        "cen": cen,
        "sizes": read_sizes(hic / f"{species}.chrom.sizes"),
        "sat": read_intervals_from_bed(repeat_root / f"{species}.TRASH" / f"{species}.TRASH_arrays.bed", chroms),
        "repeat": read_intervals_from_gff(repeat_root / f"{species}.EDTA" / f"{species}.fa.mod.EDTA.TEanno.gff3", chroms),
        "gene": read_intervals_from_gff(project / "output_all" / "01_gff" / f"{species}.gff3", chroms, feature="gene"),
        "ab": read_ab(hic / f"{species}.AB.100kb.cis.vecs.tsv", chroms),
        "tad_score": read_bedgraph(hic / f"{species}.TAD.100kb_score.bedgraph", chroms),
        "tad_boundary": read_intervals_from_bed(hic / f"{species}.TAD.100kb_boundaries.bed", chroms),
        "cool100": hic / f"{species}.100000.cool",
        "cenh3": {},
    }
    for mode in ["raw", "relaxed", "q20"]:
        data["cenh3"][mode] = read_bedgraph(recall / "bedgraph" / f"{species}.CENH3_vs_Input.{mode}.10k.smooth50k.CPM.log2.bedGraph", chroms)
    return data


def plot_global(species: str, data, outdir: Path):
    tracks = 6
    fig = plt.figure(figsize=(8.4, 0.86 * tracks * len(data["cen"])))
    gs = GridSpec(tracks * len(data["cen"]), 1, figure=fig, hspace=0.06, left=0.085, right=0.995, top=0.985, bottom=0.018)
    fig.suptitle(f"{species} whole-chromosome centromere evidence overview", fontsize=11, fontweight="bold", y=0.997)
    for i, c in enumerate(data["cen"]):
        chrom = c["chrom"]
        start, end = 0, data["sizes"][chrom]
        axes = [fig.add_subplot(gs[i * tracks + j, 0]) for j in range(tracks)]
        plot_region_tracks(axes, chrom, start, end, c, data)
        axes[0].text(0, 2.45, f"{chrom}  CENH3 {c['start']/1e6:.2f}-{c['end']/1e6:.2f} Mb", color=COL["cen"], fontsize=6.2, fontweight="bold", va="top")
        if i == 0:
            axes[0].legend(loc="upper right", ncol=3, fontsize=5.8, handlelength=1.2)
        axes[-1].set_xlabel("Position (Mb)")
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"{species}.publication_style_global.{ext}", dpi=420 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def plot_local(species: str, data, outdir: Path):
    tracks = 7
    fig = plt.figure(figsize=(8.4, 1.02 * tracks * len(data["cen"])))
    gs = GridSpec(tracks * len(data["cen"]), 1, figure=fig, height_ratios=[1, 0.72, 0.72, 0.72, 0.62, 0.62, 1.45] * len(data["cen"]), hspace=0.06, left=0.085, right=0.995, top=0.985, bottom=0.018)
    fig.suptitle(f"{species} local centromere evidence with Hi-C contact triangles", fontsize=11, fontweight="bold", y=0.997)
    for i, c in enumerate(data["cen"]):
        chrom = c["chrom"]
        mid = (c["start"] + c["end"]) // 2
        start = max(0, mid - LOCAL_FLANK)
        end = min(data["sizes"][chrom], mid + LOCAL_FLANK)
        axes = [fig.add_subplot(gs[i * tracks + j, 0]) for j in range(tracks)]
        plot_region_tracks(axes[:6], chrom, start, end, c, data)
        axes[0].text(start / 1e6, 2.45, f"{chrom}  CENH3 {c['start']/1e6:.2f}-{c['end']/1e6:.2f} Mb", color=COL["cen"], fontsize=6.2, fontweight="bold", va="top")
        if i == 0:
            axes[0].legend(loc="upper right", ncol=3, fontsize=5.8, handlelength=1.2)
        ax = axes[6]
        add_common(ax, c, (start / 1e6, end / 1e6))
        tri = hic_triangle_array(data["cool100"], chrom, start, end, max_distance_bins=100)
        if tri is not None:
            vmax = np.nanpercentile(tri, 98)
            ax.imshow(
                tri,
                aspect="auto",
                interpolation="nearest",
                origin="lower",
                extent=(start / 1e6, end / 1e6, 0, min(10, (end - start) / 2e6)),
                cmap="YlOrBr",
                vmin=0,
                vmax=max(vmax, 0.1),
            )
        ax.set_yticks([])
        ax.set_ylabel("Hi-C", rotation=0, ha="right", va="center")
        ax.set_xlabel("Position (Mb)")
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"{species}.publication_style_local.{ext}", dpi=420 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", default="path/to/project/10.centromere_analysis")
    ap.add_argument("--outdir", default="path/to/project/10.centromere_analysis/output_key/16_CENH3_multitrack_evidence/publication_style")
    ap.add_argument("--species", nargs="+", default=SPECIES)
    args = ap.parse_args()
    project = Path(args.project)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for species in args.species:
        sp_out = outdir / species
        sp_out.mkdir(parents=True, exist_ok=True)
        data = load_data(project, species)
        plot_global(species, data, sp_out)
        plot_local(species, data, sp_out)


if __name__ == "__main__":
    main()
