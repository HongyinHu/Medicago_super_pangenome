#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle


SPECIES = ["genome_Msa", "genome_474"]
BIN = 100_000
LOCAL_FLANK = 10_000_000

COL = {
    "cen_raw": "#8f8f8f",
    "cen_relaxed": "#e24a3b",
    "cen_q20": "#2f6f9f",
    "cen_box": "#8199bd",
    "sat": "#8ea9cf",
    "repeat": "#cfcfcf",
    "gene": "#d4d4d4",
    "gene_line": "#e24a3b",
    "A": "#c94737",
    "B": "#2f7f9f",
    "tad": "#333333",
    "chrom": "#bfbfbf",
    "chrom_edge": "#202020",
    "cen_span": "#efd8d8",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 8,
        "axes.linewidth": 0.5,
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
            rows.append(
                {
                    "chrom": p[0],
                    "start": int(float(p[1])),
                    "end": int(float(p[2])),
                    "name": p[3] if len(p) > 3 else "",
                    "monomer": p[4] if len(p) > 4 else "",
                }
            )
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
                v = float(p[3])
            except ValueError:
                continue
            if math.isfinite(v):
                out[p[0]].append((int(float(p[1])), int(float(p[2])), v))
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


def read_sat_bed(path: Path, chroms: set[str]) -> dict[str, list[dict]]:
    out = {c: [] for c in chroms}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] not in chroms:
                continue
            monomer = p[4] if len(p) > 4 else motif_from_name(p[3] if len(p) > 3 else "")
            out[p[0]].append({"start": int(float(p[1])), "end": int(float(p[2])), "monomer": monomer})
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


def motif_from_name(name: str) -> str:
    m = re.search(r"TRASH_(\d+)bp", name)
    return m.group(1) if m else "sat"


def overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    intervals = sorted(intervals)
    if not intervals:
        return []
    out = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = out[-1]
        if s <= le:
            out[-1] = (ls, max(le, e))
        else:
            out.append((s, e))
    return out


def coverage_track(intervals: list[tuple[int, int]], start: int, end: int, bin_size: int = BIN):
    n = max(1, math.ceil((end - start) / bin_size))
    cov = np.zeros(n)
    lengths = np.array([min(start + (i + 1) * bin_size, end) - (start + i * bin_size) for i in range(n)], dtype=float)
    for s, e in merge_intervals([(max(start, s), min(end, e)) for s, e in intervals if e > start and s < end]):
        b0 = max(0, (s - start) // bin_size)
        b1 = min(n - 1, (max(s, e - 1) - start) // bin_size)
        for b in range(int(b0), int(b1) + 1):
            bs = start + b * bin_size
            be = min(bs + bin_size, end)
            cov[b] += overlap((bs, be), (s, e))
    xs = np.array([(start + i * bin_size + min(bin_size, end - (start + i * bin_size)) / 2) / 1e6 for i in range(n)])
    return xs, cov / np.maximum(lengths, 1)


def count_track(intervals: list[tuple[int, int]], start: int, end: int, bin_size: int = BIN):
    n = max(1, math.ceil((end - start) / bin_size))
    counts = np.zeros(n)
    for s, e in intervals:
        if e <= start or s >= end:
            continue
        ss, ee = max(start, s), min(end, e)
        b0 = max(0, (ss - start) // bin_size)
        b1 = min(n - 1, (max(ss, ee - 1) - start) // bin_size)
        counts[int(b0) : int(b1) + 1] += 1
    xs = np.array([(start + i * bin_size + min(bin_size, end - (start + i * bin_size)) / 2) / 1e6 for i in range(n)])
    return xs, counts


def line_track(rows: list[tuple[int, int, float]], start: int, end: int):
    xs, ys = [], []
    for s, e, v in rows:
        if e <= start or s >= end:
            continue
        xs.append((s + e) / 2 / 1e6)
        ys.append(v)
    return np.array(xs), np.array(ys)


def region_satellite_families(rows: list[dict], start: int, end: int, n: int = 4) -> list[str]:
    cov = {}
    for r in rows:
        ov = overlap((start, end), (r["start"], r["end"]))
        if ov:
            cov[str(r["monomer"])] = cov.get(str(r["monomer"]), 0) + ov
    return [k for k, _ in sorted(cov.items(), key=lambda kv: (-kv[1], kv[0]))[:n]]


def setup_axis(ax, xlim, show_x=False):
    ax.set_xlim(*xlim)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(show_x)
    ax.tick_params(axis="x", length=2 if show_x else 0, labelsize=8, pad=1, labelbottom=show_x)


def label_axis(ax, text, x=-0.035, size=8):
    ax.text(x, 0.5, text, transform=ax.transAxes, ha="right", va="center", fontsize=size, clip_on=False)


def add_cen_span(ax, cen):
    ax.axvspan(cen["start"] / 1e6, cen["end"] / 1e6, color=COL["cen_span"], alpha=0.75, lw=0, zorder=0)


def nice_ticks(start_mb: float, end_mb: float):
    span = end_mb - start_mb
    step = 5 if span <= 60 else 20 if span <= 140 else 50
    first = math.ceil(start_mb / step) * step
    vals = np.arange(first, end_mb + 0.01, step)
    return vals


def draw_cenh3(ax, chrom, start, end, cen, data):
    setup_axis(ax, (start / 1e6, end / 1e6))
    add_cen_span(ax, cen)
    for mode, color, alpha, lw in [
        ("raw", COL["cen_raw"], 0.38, 0.48),
        ("relaxed", COL["cen_relaxed"], 0.72, 0.62),
        ("q20", COL["cen_q20"], 0.58, 0.58),
    ]:
        xs, ys = line_track(data["cenh3"][mode].get(chrom, []), start, end)
        ax.plot(xs, ys, color=color, alpha=alpha, lw=lw, label=mode, zorder=3)
    ax.axhline(0, color="#cccccc", lw=0.35)
    vals = []
    for mode in ["raw", "relaxed", "q20"]:
        vals.extend(line_track(data["cenh3"][mode].get(chrom, []), start, end)[1].tolist())
    ymax = np.nanpercentile(vals, 99.5) if vals else 2
    ax.set_ylim(-0.4, max(1.8, ymax * 1.18))
    label_axis(ax, "CENH3", size=9)
    ax.legend(loc="upper right", ncol=3, fontsize=6.5, handlelength=1.2, columnspacing=0.8)


def draw_sat_rows(axs, chrom, start, end, cen, data, family_count=4):
    families = region_satellite_families(data["sat"].get(chrom, []), start, end, family_count)
    while len(families) < family_count:
        families.append("")
    for ax, fam in zip(axs, families):
        setup_axis(ax, (start / 1e6, end / 1e6))
        add_cen_span(ax, cen)
        if fam:
            for r in data["sat"].get(chrom, []):
                if str(r["monomer"]) != fam or r["end"] <= start or r["start"] >= end:
                    continue
                xs = max(r["start"], start) / 1e6
                xe = min(r["end"], end) / 1e6
                ax.add_patch(Rectangle((xs, 0.26), xe - xs, 0.48, facecolor=COL["sat"], edgecolor="none", alpha=0.95))
            label_axis(ax, f"sat{fam}bp", size=8)
        else:
            label_axis(ax, "sat", size=8)
        ax.set_ylim(0, 1)


def draw_chrom_bar(ax, chrom, start, end, cen, species):
    xlim = (start / 1e6, end / 1e6)
    setup_axis(ax, xlim, show_x=True)
    ax.set_ylim(0, 1)
    ax.add_patch(Rectangle((xlim[0], 0.38), xlim[1] - xlim[0], 0.24, facecolor=COL["chrom"], edgecolor=COL["chrom_edge"], lw=0.7))
    ax.add_patch(Rectangle((cen["start"] / 1e6, 0.38), (cen["end"] - cen["start"]) / 1e6, 0.24, facecolor=COL["cen_box"], edgecolor="none", alpha=0.85))
    ax.set_xticks(nice_ticks(*xlim))
    ax.set_xticklabels([f"{x:g}Mb" for x in nice_ticks(*xlim)])
    label_axis(ax, chrom, size=9)
    ax.text(xlim[0], 0.92, f"{species} {chrom}: CENH3 {cen['start']/1e6:.2f}-{cen['end']/1e6:.2f} Mb", ha="left", va="top", fontsize=7, color="#b33b3b")


def draw_repeat(ax, chrom, start, end, cen, data):
    setup_axis(ax, (start / 1e6, end / 1e6))
    add_cen_span(ax, cen)
    xs, ys = coverage_track(data["repeat"].get(chrom, []), start, end)
    ax.fill_between(xs, 0, ys, color=COL["repeat"], step="mid", lw=0)
    ax.set_ylim(0, 1)
    label_axis(ax, "Repeat Density", size=9)


def draw_gene(ax, chrom, start, end, cen, data):
    setup_axis(ax, (start / 1e6, end / 1e6))
    add_cen_span(ax, cen)
    xs, ys = count_track(data["gene"].get(chrom, []), start, end)
    ymax = max(1, float(np.nanpercentile(ys, 98)) if len(ys) else 1)
    yscaled = np.minimum(ys / ymax, 1.3)
    ax.fill_between(xs, 0, yscaled, color=COL["gene"], step="mid", lw=0)
    ax.plot(xs, yscaled, color=COL["gene_line"], lw=0.35, alpha=0.55)
    ax.set_ylim(0, 1.35)
    label_axis(ax, "Gene Density", size=9)


def draw_ab(ax, chrom, start, end, cen, data):
    setup_axis(ax, (start / 1e6, end / 1e6))
    add_cen_span(ax, cen)
    vals = [(s, e, v) for s, e, v in data["ab"].get(chrom, []) if e > start and s < end]
    m = max([abs(v) for _, _, v in vals] + [0.2])
    for s, e, v in vals:
        x0 = max(s, start) / 1e6
        w = (min(e, end) - max(s, start)) / 1e6
        ax.add_patch(Rectangle((x0, min(0, v)), w, abs(v), facecolor=COL["A"] if v >= 0 else COL["B"], edgecolor="none", alpha=0.92))
    ax.axhline(0, color="#b8b8b8", lw=0.35)
    ax.set_ylim(-m * 1.05, m * 1.05)
    label_axis(ax, "Eigv", size=9)
    ax.text(-0.095, -0.23, "A", transform=ax.transAxes, ha="left", va="center", fontsize=7, clip_on=False)
    ax.add_patch(Rectangle((-0.128, -0.275), 0.026, 0.08, transform=ax.transAxes, facecolor=COL["A"], edgecolor="none", clip_on=False))
    ax.text(-0.047, -0.23, "B", transform=ax.transAxes, ha="left", va="center", fontsize=7, clip_on=False)
    ax.add_patch(Rectangle((-0.08, -0.275), 0.026, 0.08, transform=ax.transAxes, facecolor=COL["B"], edgecolor="none", clip_on=False))


def downsample_matrix(mat: np.ndarray, target_bins: int = 620):
    n = mat.shape[0]
    if n <= target_bins:
        return mat
    f = math.ceil(n / target_bins)
    nn = (n // f) * f
    mat = mat[:nn, :nn]
    with np.errstate(invalid="ignore"):
        small = np.nanmean(mat.reshape(nn // f, f, nn // f, f), axis=(1, 3))
    return small


def fetch_hic(cool_path: Path, chrom: str, start: int, end: int, target_bins: int = 620):
    try:
        import cooler
    except Exception:
        return None
    try:
        clr = cooler.Cooler(str(cool_path))
        mat = clr.matrix(balance=False).fetch(f"{chrom}:{start}-{end}")
    except Exception:
        return None
    mat = np.asarray(mat, dtype=float)
    mat[~np.isfinite(mat)] = np.nan
    return downsample_matrix(mat, target_bins=target_bins)


def draw_hic_triangle(ax, chrom, start, end, cen, data, target_bins=620):
    setup_axis(ax, (start / 1e6, end / 1e6))
    ax.spines["bottom"].set_visible(False)
    mat = fetch_hic(data["cool100"], chrom, start, end, target_bins=target_bins)
    if mat is None or mat.size == 0:
        label_axis(ax, "TAD", size=9)
        ax.text(0.5, 0.5, "Hi-C matrix unavailable", transform=ax.transAxes, ha="center", va="center", fontsize=8, color="#777777")
        return
    positive = mat[np.isfinite(mat) & (mat > 0)]
    if positive.size == 0:
        label_axis(ax, "TAD", size=9)
        ax.text(0.5, 0.5, "Hi-C matrix unavailable", transform=ax.transAxes, ha="center", va="center", fontsize=8, color="#777777")
        return
    vmax = max(5, float(np.nanpercentile(positive, 99.2)))
    vmin = max(1, float(np.nanpercentile(positive, 8)))
    n = mat.shape[0]
    span = (end - start) / 1e6
    bin_mb = span / n
    xs, ys, cs = [], [], []
    for d in range(n):
        vals = np.diag(mat, k=d)
        if not len(vals):
            continue
        i = np.arange(len(vals))
        x = start / 1e6 + (i + d / 2 + 0.5) * bin_mb
        y = d * bin_mb * 0.5
        good = np.isfinite(vals) & (vals > 0)
        xs.extend(x[good].tolist())
        ys.extend([y] * int(np.sum(good)))
        cs.extend(vals[good].tolist())
    size = max(0.5, min(7.0, 14500 / max(n, 1)))
    sc = ax.scatter(
        xs,
        ys,
        c=cs,
        cmap="YlOrBr",
        s=size,
        marker="s",
        linewidths=0,
        norm=LogNorm(vmin=vmin, vmax=max(vmax, vmin * 1.2)),
        rasterized=True,
    )
    ax.set_ylim((n - 1) * bin_mb * 0.5, 0)
    ax.set_yticks([])
    label_axis(ax, "TAD", size=9)
    x0, x1 = start / 1e6, end / 1e6
    ax.plot([x0, (x0 + x1) / 2, x1], [0, (n - 1) * bin_mb * 0.5, 0], color="white", lw=0.25, alpha=0.65)
    cb = ax.inset_axes([1.012, 0.22, 0.012, 0.56])
    cbar = plt.colorbar(sc, cax=cb)
    desired_ticks = [1, 5, 20, 50, 100]
    ticks = [t for t in desired_ticks if vmin <= t <= vmax]
    if ticks:
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([str(t) for t in ticks])
    cbar.ax.tick_params(labelsize=6, length=2)


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
        "sat": read_sat_bed(repeat_root / f"{species}.TRASH" / f"{species}.TRASH_arrays.bed", chroms),
        "repeat": read_intervals_from_gff(repeat_root / f"{species}.EDTA" / f"{species}.fa.mod.EDTA.TEanno.gff3", chroms),
        "gene": read_intervals_from_gff(project / "output_all" / "01_gff" / f"{species}.gff3", chroms, feature="gene"),
        "ab": read_ab(hic / f"{species}.AB.100kb.cis.vecs.tsv", chroms),
        "cool100": hic / f"{species}.100000.cool",
        "cenh3": {},
    }
    for mode in ["raw", "relaxed", "q20"]:
        data["cenh3"][mode] = read_bedgraph(recall / "bedgraph" / f"{species}.CENH3_vs_Input.{mode}.10k.smooth50k.CPM.log2.bedGraph", chroms)
    return data


def draw_brapa_style(species: str, data, cen: dict, out_prefix: Path, start: int, end: int, panel_letter: str):
    chrom = cen["chrom"]
    start = max(0, start)
    end = min(data["sizes"][chrom], end)
    fig = plt.figure(figsize=(10.2, 4.25))
    gs = GridSpec(
        10,
        1,
        figure=fig,
        height_ratios=[0.62, 0.22, 0.22, 0.22, 0.22, 0.42, 0.42, 0.42, 0.58, 2.62],
        hspace=0.02,
        left=0.16,
        right=0.93,
        top=0.94,
        bottom=0.08,
    )
    axes = [fig.add_subplot(gs[i, 0]) for i in range(10)]
    fig.text(0.055, 0.955, panel_letter, fontsize=16, fontweight="bold", ha="left", va="top")
    draw_cenh3(axes[0], chrom, start, end, cen, data)
    draw_sat_rows(axes[1:5], chrom, start, end, cen, data, family_count=4)
    draw_chrom_bar(axes[5], chrom, start, end, cen, species)
    draw_repeat(axes[6], chrom, start, end, cen, data)
    draw_gene(axes[7], chrom, start, end, cen, data)
    draw_ab(axes[8], chrom, start, end, cen, data)
    draw_hic_triangle(axes[9], chrom, start, end, cen, data, target_bins=620)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{out_prefix}.{ext}", dpi=520 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", default="path/to/project/10.centromere_analysis")
    ap.add_argument("--outdir", default="path/to/project/10.centromere_analysis/output_key/16_CENH3_multitrack_evidence/brapa_style")
    ap.add_argument("--species", nargs="+", default=SPECIES)
    ap.add_argument("--chrom", default=None, help="Optional chromosome, for example Chr1. By default all CENH3 chromosomes are drawn.")
    args = ap.parse_args()

    project = Path(args.project)
    outdir = Path(args.outdir)
    for species in args.species:
        data = load_data(project, species)
        sp_out = outdir / species
        sp_out.mkdir(parents=True, exist_ok=True)
        cens = [c for c in data["cen"] if args.chrom is None or c["chrom"] == args.chrom]
        for idx, cen in enumerate(cens, 1):
            chrom = cen["chrom"]
            draw_brapa_style(
                species,
                data,
                cen,
                sp_out / f"{species}.{chrom}.brapa_style_global",
                0,
                data["sizes"][chrom],
                "A",
            )
            mid = (cen["start"] + cen["end"]) // 2
            draw_brapa_style(
                species,
                data,
                cen,
                sp_out / f"{species}.{chrom}.brapa_style_local",
                mid - LOCAL_FLANK,
                mid + LOCAL_FLANK,
                "B",
            )


if __name__ == "__main__":
    main()
