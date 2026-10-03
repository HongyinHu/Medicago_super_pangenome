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
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle


SPECIES = ["genome_Msa", "genome_474"]
WINDOW = 100_000
FLANK = 10_000_000

COL = {
    "raw": "#9a9a9a",
    "relaxed": "#4d9b5f",
    "q20": "#2f5f9e",
    "cen": "#bf4b55",
    "sat": "#7463a9",
    "repeat": "#b8b8b8",
    "gene": "#d84b4b",
    "A": "#c64b3b",
    "B": "#3c7da0",
    "tad": "#5f5f5f",
    "axis": "#9aa3ad",
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
            rows.append(
                {
                    "chrom": p[0],
                    "start": int(float(p[1])),
                    "end": int(float(p[2])),
                    "name": p[3] if len(p) > 3 else "",
                    "extra": p[4:] if len(p) > 4 else [],
                }
            )
    return rows


def read_gene_intervals(path: Path, chroms: set[str]) -> dict[str, list[tuple[int, int]]]:
    out = {c: [] for c in chroms}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 5 or p[2] != "gene" or p[0] not in chroms:
                continue
            out[p[0]].append((int(p[3]) - 1, int(p[4])))
    return out


def read_repeat_intervals(path: Path, chroms: set[str]) -> dict[str, list[tuple[int, int]]]:
    out = {c: [] for c in chroms}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] not in chroms:
                continue
            out[p[0]].append((int(p[3]) - 1, int(p[4])))
    return out


def read_bed_intervals(path: Path, chroms: set[str]) -> dict[str, list[tuple[int, int]]]:
    out = {c: [] for c in chroms}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] not in chroms:
                continue
            out[p[0]].append((int(float(p[1])), int(float(p[2]))))
    return out


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


def read_ab(path: Path, chroms: set[str]) -> dict[str, list[tuple[int, int, float]]]:
    out = {c: [] for c in chroms}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row.get("chrom", "")
            if chrom not in chroms:
                continue
            value_s = row.get("E1", "")
            if not value_s:
                continue
            try:
                value = float(value_s)
            except ValueError:
                continue
            out[chrom].append((int(row["start"]), int(row["end"]), value))
    return out


def overlap(a: tuple[int, int], b: tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def coverage_track(intervals: list[tuple[int, int]], start: int, end: int, window: int = WINDOW) -> list[tuple[float, float]]:
    intervals = merge_intervals([(max(start, s), min(end, e)) for s, e in intervals if e > start and s < end])
    values = []
    pos = start
    while pos < end:
        wend = min(pos + window, end)
        cov = sum(overlap((pos, wend), iv) for iv in intervals)
        values.append(((pos + wend) / 2 / 1e6, cov / max(1, wend - pos)))
        pos = wend
    return values


def count_track(intervals: list[tuple[int, int]], start: int, end: int, window: int = WINDOW) -> list[tuple[float, float]]:
    values = []
    pos = start
    while pos < end:
        wend = min(pos + window, end)
        count = sum(1 for s, e in intervals if e > pos and s < wend)
        values.append(((pos + wend) / 2 / 1e6, count))
        pos = wend
    return values


def bedgraph_track(rows: list[tuple[int, int, float]], start: int, end: int) -> tuple[list[float], list[float]]:
    xs, ys = [], []
    for s, e, v in rows:
        if e <= start or s >= end:
            continue
        xs.append((s + e) / 2 / 1e6)
        ys.append(v)
    return xs, ys


def write_track_table(path: Path, species: str, cen_rows: list[dict], sat, repeat, genes, ab, tad_score, tad_boundaries) -> None:
    with path.open("w", newline="") as fh:
        fields = [
            "species",
            "chrom",
            "cen_start",
            "cen_end",
            "focus_start",
            "focus_end",
            "bin_start",
            "bin_end",
            "satellite_fraction",
            "repeat_fraction",
            "gene_count",
            "AB_E1",
            "TAD_score",
            "TAD_boundary_count",
        ]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for c in cen_rows:
            chrom = c["chrom"]
            focus_start = max(0, (c["start"] + c["end"]) // 2 - FLANK)
            focus_end = (c["start"] + c["end"]) // 2 + FLANK
            sat_t = dict(coverage_track(sat.get(chrom, []), focus_start, focus_end))
            rep_t = dict(coverage_track(repeat.get(chrom, []), focus_start, focus_end))
            gene_t = dict(count_track(genes.get(chrom, []), focus_start, focus_end))
            ab_t = {((s + e) / 2 / 1e6): v for s, e, v in ab.get(chrom, []) if e > focus_start and s < focus_end}
            tad_t = {((s + e) / 2 / 1e6): v for s, e, v in tad_score.get(chrom, []) if e > focus_start and s < focus_end}
            boundary_intervals = tad_boundaries.get(chrom, [])
            pos = focus_start
            while pos < focus_end:
                wend = min(pos + WINDOW, focus_end)
                mid = (pos + wend) / 2 / 1e6
                boundary_count = sum(1 for s, e in boundary_intervals if e > pos and s < wend)
                writer.writerow(
                    {
                        "species": species,
                        "chrom": chrom,
                        "cen_start": c["start"],
                        "cen_end": c["end"],
                        "focus_start": focus_start,
                        "focus_end": focus_end,
                        "bin_start": pos,
                        "bin_end": wend,
                        "satellite_fraction": sat_t.get(mid, ""),
                        "repeat_fraction": rep_t.get(mid, ""),
                        "gene_count": gene_t.get(mid, ""),
                        "AB_E1": ab_t.get(mid, ""),
                        "TAD_score": tad_t.get(mid, ""),
                        "TAD_boundary_count": boundary_count,
                    }
                )
                pos = wend


def plot_species(species: str, root: Path, out_root: Path) -> None:
    project = root
    recall = project / "output_key" / "15_recall_Msa_474_CENH3" / species
    hic = project / "output_all" / "08_hic_AB_TAD" / species
    repeat_dir = project / "output_all" / "03_repeat"
    gff_path = project / "output_all" / "01_gff" / f"{species}.gff3"
    edta_path = repeat_dir / f"{species}.EDTA" / f"{species}.fa.mod.EDTA.TEanno.gff3"
    trash_path = repeat_dir / f"{species}.TRASH" / f"{species}.TRASH_arrays.bed"

    cen_rows = read_bed(recall / "results" / f"{species}.CENH3.functional_centromere.recall.high_confidence.bed")
    cen_rows = sorted(cen_rows, key=lambda r: int(r["chrom"].replace("Chr", "")))
    chroms = {r["chrom"] for r in cen_rows}

    genes = read_gene_intervals(gff_path, chroms)
    repeats = read_repeat_intervals(edta_path, chroms)
    sats = read_bed_intervals(trash_path, chroms)
    ab = read_ab(hic / f"{species}.AB.100kb.cis.vecs.tsv", chroms)
    tad_score = read_bedgraph(hic / f"{species}.TAD.100kb_score.bedgraph", chroms)
    tad_boundaries = read_bed_intervals(hic / f"{species}.TAD.100kb_boundaries.bed", chroms)
    cenh3 = {
        mode: read_bedgraph(
            recall / "bedgraph" / f"{species}.CENH3_vs_Input.{mode}.10k.smooth50k.CPM.log2.bedGraph",
            chroms,
        )
        for mode in ["raw", "relaxed", "q20"]
    }

    outdir = out_root / species
    outdir.mkdir(parents=True, exist_ok=True)
    write_track_table(outdir / f"{species}.CENH3_multitrack_focus_100kb.tsv", species, cen_rows, sats, repeats, genes, ab, tad_score, tad_boundaries)

    tracks_per_chr = 6
    fig = plt.figure(figsize=(8.2, 0.94 * tracks_per_chr * len(cen_rows)))
    gs = GridSpec(
        tracks_per_chr * len(cen_rows),
        1,
        figure=fig,
        hspace=0.08,
        left=0.085,
        right=0.995,
        top=0.985,
        bottom=0.018,
    )
    fig.suptitle(
        f"{species} centromere multi-track evidence around high-confidence CENH3 domains",
        fontsize=11,
        fontweight="bold",
        y=0.997,
    )

    legend_handles = []
    for i, c in enumerate(cen_rows):
        chrom = c["chrom"]
        cen_mid = (c["start"] + c["end"]) // 2
        focus_start = max(0, cen_mid - FLANK)
        focus_end = cen_mid + FLANK
        xlim = (focus_start / 1e6, focus_end / 1e6)
        axes = [fig.add_subplot(gs[i * tracks_per_chr + j, 0]) for j in range(tracks_per_chr)]
        for ax in axes:
            ax.set_xlim(*xlim)
            ax.axvspan(c["start"] / 1e6, c["end"] / 1e6, color=COL["cen"], alpha=0.16, lw=0)
            ax.spines["left"].set_color(COL["axis"])
            ax.spines["bottom"].set_color(COL["axis"])
            ax.tick_params(axis="both", labelsize=5.7, length=2, pad=1)
        for ax in axes[:-1]:
            ax.set_xticklabels([])

        # CENH3
        ax = axes[0]
        for mode, alpha, lw in [("raw", 0.25, 0.45), ("relaxed", 0.45, 0.55), ("q20", 0.58, 0.62)]:
            xs, ys = bedgraph_track(cenh3[mode].get(chrom, []), focus_start, focus_end)
            ax.plot(xs, ys, color=COL[mode], alpha=alpha, lw=lw, label=mode)
        ax.axhline(0, color="#c8cdd2", lw=0.4)
        ax.axhline(0.5, color="#d8a2a8", lw=0.45, ls=(0, (3, 2)))
        ax.set_ylim(-1.1, 2.6)
        ax.set_ylabel("CENH3", rotation=0, ha="right", va="center", fontsize=6.2)
        ax.text(
            xlim[0],
            2.45,
            f"{chrom}  CENH3 {c['start']/1e6:.2f}-{c['end']/1e6:.2f} Mb",
            ha="left",
            va="top",
            fontsize=6.4,
            color=COL["cen"],
            fontweight="bold",
        )
        if i == 0:
            ax.legend(loc="upper right", ncol=3, fontsize=5.8, handlelength=1.2)

        # Satellite
        ax = axes[1]
        sat_vals = coverage_track(sats.get(chrom, []), focus_start, focus_end)
        ax.fill_between([x for x, _ in sat_vals], [v for _, v in sat_vals], color=COL["sat"], alpha=0.72, step="mid")
        ax.set_ylim(0, 1)
        ax.set_ylabel("sat", rotation=0, ha="right", va="center", fontsize=6.2)

        # Repeat
        ax = axes[2]
        rep_vals = coverage_track(repeats.get(chrom, []), focus_start, focus_end)
        ax.fill_between([x for x, _ in rep_vals], [v for _, v in rep_vals], color=COL["repeat"], alpha=0.85, step="mid")
        ax.set_ylim(0, 1)
        ax.set_ylabel("repeat", rotation=0, ha="right", va="center", fontsize=6.2)

        # Gene
        ax = axes[3]
        gene_vals = count_track(genes.get(chrom, []), focus_start, focus_end)
        ax.plot([x for x, _ in gene_vals], [v for _, v in gene_vals], color=COL["gene"], lw=0.75, alpha=0.85)
        max_gene = max([v for _, v in gene_vals] + [1])
        ax.set_ylim(0, max(2, max_gene * 1.15))
        ax.set_ylabel("gene", rotation=0, ha="right", va="center", fontsize=6.2)

        # AB
        ax = axes[4]
        ab_rows = [(s, e, v) for s, e, v in ab.get(chrom, []) if e > focus_start and s < focus_end]
        for s, e, v in ab_rows:
            ax.add_patch(
                Rectangle(
                    (max(s, focus_start) / 1e6, 0 if v >= 0 else v),
                    (min(e, focus_end) - max(s, focus_start)) / 1e6,
                    abs(v),
                    facecolor=COL["A"] if v >= 0 else COL["B"],
                    edgecolor="none",
                    alpha=0.82,
                )
            )
        ax.axhline(0, color="#c8cdd2", lw=0.45)
        vals = [v for _, _, v in ab_rows]
        m = max([abs(v) for v in vals] + [0.2])
        ax.set_ylim(-m * 1.08, m * 1.08)
        ax.set_ylabel("A/B", rotation=0, ha="right", va="center", fontsize=6.2)

        # TAD score and boundaries
        ax = axes[5]
        tx, ty = bedgraph_track(tad_score.get(chrom, []), focus_start, focus_end)
        if ty:
            ax.plot(tx, ty, color=COL["tad"], lw=0.55, alpha=0.85)
        ymin, ymax = (min(ty), max(ty)) if ty else (-1, 1)
        if ymin == ymax:
            ymin, ymax = ymin - 1, ymax + 1
        for s, e in tad_boundaries.get(chrom, []):
            if e > focus_start and s < focus_end:
                ax.axvline((s + e) / 2 / 1e6, color="#222222", lw=0.28, alpha=0.65)
        ax.set_ylim(ymin - 0.08 * abs(ymax - ymin), ymax + 0.08 * abs(ymax - ymin))
        ax.set_ylabel("TAD", rotation=0, ha="right", va="center", fontsize=6.2)
        ax.set_xlabel("Position (Mb)", fontsize=6.3)

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(outdir / f"{species}.CENH3_multitrack_focus.{ext}", dpi=420 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default="path/to/project/10.centromere_analysis")
    parser.add_argument("--outdir", default="path/to/project/10.centromere_analysis/output_key/16_CENH3_multitrack_evidence")
    parser.add_argument("--species", nargs="+", default=SPECIES)
    args = parser.parse_args()

    root = Path(args.project)
    out_root = Path(args.outdir)
    out_root.mkdir(parents=True, exist_ok=True)
    for species in args.species:
        plot_species(species, root, out_root)


if __name__ == "__main__":
    main()
