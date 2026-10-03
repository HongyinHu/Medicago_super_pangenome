#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

try:
    import cooler
except Exception as exc:  # pragma: no cover
    cooler = None
    COOLER_IMPORT_ERROR = exc
else:
    COOLER_IMPORT_ERROR = None


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    }
)


ROOT = Path("path/to/project/10.centromere_analysis/output_key")
HIC_ROOT = Path("path/to/project/10.centromere_analysis/output_all/08_hic_AB_TAD")
OUTDIR = ROOT / "09_karyotype_model"

FATE_TSV = ROOT / "05_fusion_chr_CEN_fate/Mpo_ancestral_centromere_fate_integrated.tsv"
CONTEXT_TSV = ROOT / "05_fusion_chr_CEN_fate/Mpo_Chr4_CEN4_CEN5_fusion_context.tsv"
LOCAL_BLOCKS_TSV = ROOT / "05_fusion_chr_CEN_fate/Mpo_Chr4_16_26Mb_R108_primary_chr_sequence_blocks_by_chr_summary.tsv"
MULTI_TSV = ROOT / "08_multispecies_validation/multispecies_CEN5_species_summary.tsv"

COL = {
    "cen4": "#3b6fb6",
    "cen5": "#e6862d",
    "active": "#4d9b5f",
    "missing": "#d95b5b",
    "grey": "#6f7680",
    "lightgrey": "#e8ebef",
    "text": "#222222",
    "chr4": "#4e79a7",
    "chr5": "#f28e2b",
    "matrix": "#b84c4c",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_bed(path: Path, cols: list[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=cols)
    return pd.read_csv(path, sep="\t", header=None, names=cols, comment="#")


def fetch_cool_matrix(cool_path: Path, chrom: str, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
    if cooler is None:
        raise RuntimeError(f"cooler import failed: {COOLER_IMPORT_ERROR}")
    c = cooler.Cooler(str(cool_path))
    region = f"{chrom}:{start}-{end}"
    # Raw contacts give a more stable visual signal for these pericentromeric
    # windows, where balanced matrices can be very sparse after filtering.
    mat = c.matrix(balance=False).fetch(region)
    if np.isnan(mat).all():
        mat = c.matrix(balance=False).fetch(region)
    mat = np.asarray(mat, dtype=float)
    mat[~np.isfinite(mat)] = 0
    if np.nanmax(mat) > 0:
        mat = np.log1p(mat)
        hi = np.nanpercentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
        if hi > 0:
            mat = np.clip(mat, 0, hi)
    bins = c.bins().fetch(region)
    centers = ((bins["start"].to_numpy() + bins["end"].to_numpy()) / 2.0) / 1e6
    return mat, centers


def subset_intervals(df: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    if df.empty:
        return df
    sub = df[(df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)].copy()
    return sub


def plot_heatmap(ax, genome: str, chrom: str, start: int, end: int, title: str, domains: pd.DataFrame, highlight: list[tuple[int, int, str, str]]):
    cool_path = HIC_ROOT / genome / f"{genome}.40000.cool"
    mat, centers = fetch_cool_matrix(cool_path, chrom, start, end)
    extent = [start / 1e6, end / 1e6, end / 1e6, start / 1e6]
    ax.imshow(mat, extent=extent, cmap="Reds", aspect="auto", interpolation="nearest")
    ax.set_title(title, loc="left", fontsize=8.5, fontweight="bold", pad=2)
    ax.set_ylabel(f"{chrom} (Mb)", fontsize=7)
    ax.tick_params(labelsize=7, length=2)
    ax.spines["left"].set_color("#b6bbc2")
    ax.spines["bottom"].set_color("#b6bbc2")

    sub = subset_intervals(domains, chrom, start, end)
    for _, r in sub.iterrows():
        s = max(int(r["start"]), start) / 1e6
        e = min(int(r["end"]), end) / 1e6
        w = e - s
        if w <= 0:
            continue
        ax.add_patch(Rectangle((s, s), w, w, fill=False, edgecolor="#57935f", lw=0.55, alpha=0.70))

    y_top = start / 1e6 - 0.08 * ((end - start) / 1e6)
    for hs, he, color, label in highlight:
        ax.add_patch(Rectangle((hs / 1e6, y_top), (he - hs) / 1e6, 0.05 * ((end - start) / 1e6), facecolor=color, edgecolor="none", clip_on=False))
        ax.text((hs + he) / 2 / 1e6, y_top - 0.05 * ((end - start) / 1e6), label, ha="center", va="top", fontsize=7.1, color=color, clip_on=False)


def plot_score_track(ax, score_df: pd.DataFrame, chrom: str, start: int, end: int, label: str):
    sub = subset_intervals(score_df, chrom, start, end)
    if sub.empty:
        ax.text(0.5, 0.5, "TAD score unavailable", transform=ax.transAxes, ha="center", va="center", color=COL["grey"], fontsize=7)
        return
    x = ((sub["start"].to_numpy() + sub["end"].to_numpy()) / 2.0) / 1e6
    y = sub["score"].astype(float).to_numpy()
    ax.plot(x, y, color="#334f74", lw=0.8)
    ax.axhline(0, color="#c2c7ce", lw=0.5)
    ax.set_xlim(start / 1e6, end / 1e6)
    lim = np.nanmax(np.abs(y)) if np.any(np.isfinite(y)) else 1
    ax.set_ylim(-lim * 1.12, lim * 1.12)
    ax.set_ylabel(label, fontsize=7)
    ax.tick_params(labelsize=7, length=2)
    ax.spines["left"].set_color("#b6bbc2")
    ax.spines["bottom"].set_color("#b6bbc2")


def plot_linear_context(ax, context, local_blocks):
    ax.set_xlim(16.0, 26.5)
    ax.set_ylim(0, 5)
    ax.set_yticks([])
    ax.set_xlabel("Mpo Chr4 position (Mb)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7, length=2)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#b6bbc2")

    ctx = {r["feature"]: r for r in context}
    local = {r["R108_chr"]: r for r in local_blocks}
    active = ctx["active_Mpo_Chr4_CENH3_core"]
    c5 = ctx["R108_CEN5_top_projection_cluster_on_Mpo_Chr4"]

    c5_s, c5_e = int(c5["start"]) / 1e6, int(c5["end"]) / 1e6
    active_s, active_e = int(active["start"]) / 1e6, int(active["end"]) / 1e6

    ax.hlines(4.0, 16.0, 26.5, color="#c7ccd3", lw=7)
    ax.add_patch(Rectangle((c5_s, 3.72), c5_e - c5_s, 0.56, facecolor=COL["cen5"], edgecolor="none"))
    ax.add_patch(Rectangle((active_s, 3.72), active_e - active_s, 0.56, facecolor=COL["cen4"], edgecolor="none"))
    ax.text((c5_s + c5_e) / 2, 4.55, "CEN5-associated\nremnant", ha="center", va="bottom", fontsize=7, color=COL["cen5"])
    ax.text((active_s + active_e) / 2, 4.55, "active CEN4-like\nCENH3 core", ha="center", va="bottom", fontsize=7, color=COL["cen4"])

    for y, key, color in [(2.65, "Chr5", COL["chr5"]), (1.85, "Chr4", COL["chr4"])]:
        row = local[key]
        s, e = int(row["Mpo_span_start"]) / 1e6, int(row["Mpo_span_end"]) / 1e6
        ax.hlines(y, s, e, color=color, lw=5)
        ax.text(
            16.05,
            y + 0.20,
            f"R108 {key} sequence blocks: {row['block_count']} blocks, {int(row['sum_query_len'])/1000:.0f} kb",
            fontsize=7,
            color=color,
        )
    ax.text(16.05, 0.65, "All coordinates use the current R108-renamed Mpo assembly.", fontsize=7.2, color=COL["grey"])


def plot_cen_fate(ax, fate, multi):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.95, "Centromere projection evidence", fontsize=10, fontweight="bold")
    cens = [f"CEN{i}" for i in range(1, 9)]
    fate_by_cen = {r["ancestral_centromere_id"]: r for r in fate}
    x0, dx = 0.16, 0.095
    for i, cen in enumerate(cens):
        x = x0 + i * dx
        retained = fate_by_cen[cen]["classification"] == "retained_active_centromere"
        color = COL["active"] if retained else COL["missing"]
        ax.text(x, 0.79, cen, ha="center", fontsize=7.4)
        ax.add_patch(Circle((x, 0.64), 0.034, facecolor=COL["active"], edgecolor="white", lw=0.5))
        ax.add_patch(Circle((x, 0.48), 0.034, facecolor=color, edgecolor="white", lw=0.5))
        ax.text(x, 0.48, "+" if retained else "-", ha="center", va="center", fontsize=7.4, color="white", fontweight="bold")
        if cen == "CEN5":
            ax.add_patch(Rectangle((x - 0.045, 0.40), 0.09, 0.42, fill=False, edgecolor=COL["cen5"], lw=1))
    ax.text(0.01, 0.64, "x=8 refs", ha="left", va="center", fontsize=7.4)
    ax.text(0.01, 0.48, "Mpo active", ha="left", va="center", fontsize=7.4)

    multi_by_species = {r["species"]: r for r in multi}
    ax.text(0.01, 0.30, "CEN5 validation", ha="left", fontsize=7.5, fontweight="bold")
    y = 0.22
    for sp in ["R108", "A17", "474", "Msa", "Mpo"]:
        if sp not in multi_by_species:
            continue
        val = multi_by_species[sp]["CEN5_active_or_equivalent"]
        color = COL["active"] if val.startswith("yes") else COL["missing"]
        ax.add_patch(Circle((0.12, y), 0.020, facecolor=color, edgecolor="white", lw=0.5))
        ax.text(0.16, y, f"{sp}: {val}", va="center", fontsize=7.2, color=COL["text"])
        y -= 0.052


def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=10, lw=1.1, color=COL["grey"]))


def plot_model(ax):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.95, "TAD-guided working model", fontsize=10, fontweight="bold")
    ax.text(0.04, 0.82, "Ancestral x=8", fontsize=8.0, fontweight="bold")
    ax.add_patch(Rectangle((0.08, 0.58), 0.030, 0.18, facecolor=COL["cen4"], edgecolor="white"))
    ax.add_patch(Rectangle((0.16, 0.58), 0.030, 0.18, facecolor=COL["cen5"], edgecolor="white"))
    ax.text(0.095, 0.52, "Chr4\nCEN4", ha="center", fontsize=6.8, color=COL["cen4"])
    ax.text(0.175, 0.52, "Chr5\nCEN5", ha="center", fontsize=6.8, color=COL["cen5"])

    arrow(ax, 0.24, 0.67, 0.39, 0.67)
    ax.text(0.31, 0.75, "reorganization", ha="center", fontsize=7.2, color=COL["grey"])

    ax.text(0.39, 0.82, "Mpo Chr4", fontsize=8.0, fontweight="bold")
    ax.add_patch(Rectangle((0.40, 0.62), 0.28, 0.055, facecolor="#c7ccd3", edgecolor="none"))
    ax.add_patch(Rectangle((0.44, 0.608), 0.030, 0.079, facecolor=COL["cen5"], edgecolor="none"))
    ax.add_patch(Rectangle((0.55, 0.608), 0.095, 0.079, facecolor=COL["cen4"], edgecolor="none"))
    ax.text(0.455, 0.52, "CEN5\nremnant", ha="center", fontsize=6.8, color=COL["cen5"])
    ax.text(0.598, 0.52, "active\nCEN4-like", ha="center", fontsize=6.8, color=COL["cen4"])

    arrow(ax, 0.72, 0.67, 0.84, 0.67)
    ax.text(0.79, 0.74, "x=7", ha="center", fontsize=7, color=COL["grey"])
    xs = [0.83, 0.87, 0.91, 0.95, 0.85, 0.89, 0.93]
    ys = [0.66, 0.66, 0.66, 0.66, 0.56, 0.56, 0.56]
    labels = ["1", "2", "3", "4", "6", "7", "8"]
    for x, y, lab in zip(xs, ys, labels):
        ax.add_patch(Circle((x, y), 0.017, facecolor=COL["active"], edgecolor="white", lw=0.4))
        ax.text(x, y - 0.045, lab, ha="center", fontsize=6.5)
    ax.add_patch(Circle((0.97, 0.56), 0.017, facecolor=COL["missing"], edgecolor="white", lw=0.4))
    ax.text(0.97, 0.515, "5", ha="center", fontsize=6.5, color=COL["missing"])
    ax.text(0.40, 0.27, "Cautious interpretation: CEN5 loss, inactivation or unresolved repositioning remain alternatives until direct breakpoint validation is finalized.", fontsize=7.1, color=COL["grey"])


def write_tad_summary(context, local_blocks, fate, multi) -> Path:
    ctx = {r["feature"]: r for r in context}
    c5 = ctx["R108_CEN5_top_projection_cluster_on_Mpo_Chr4"]
    active = ctx["active_Mpo_Chr4_CENH3_core"]
    rows = [
        {
            "item": "coordinate_system",
            "value": "current R108-renamed Mpo assembly coordinates",
            "interpretation": "Avoids relying on older chromosome-painting PDFs with potentially outdated chromosome order or IDs.",
        },
        {
            "item": "Mpo_CEN5_candidate_remnant",
            "value": f"{c5['Mpo_chr']}:{c5['start']}-{c5['end']}",
            "interpretation": "Local R108 CEN5-associated sequence cluster on Mpo Chr4.",
        },
        {
            "item": "Mpo_active_CEN4_like_core",
            "value": f"{active['Mpo_chr']}:{active['start']}-{active['end']}",
            "interpretation": "Primary active Mpo Chr4 CENH3 core reciprocally linked to ancestral CEN4.",
        },
        {
            "item": "TAD_evidence_scope",
            "value": "Mpo 40kb raw Hi-C matrix plus 40kb TAD domains and score track",
            "interpretation": "Supports the local chromatin context of the CEN5 remnant and active CEN4-like core; it does not by itself prove the exact fusion breakpoint.",
        },
    ]
    out = OUTDIR / "Mpo_8to7_TAD_based_model_evidence.tsv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return out


def main() -> None:
    if cooler is None:
        raise SystemExit(f"cooler is required for this TAD figure: {COOLER_IMPORT_ERROR}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fate = read_tsv(FATE_TSV)
    context = read_tsv(CONTEXT_TSV)
    local_blocks = read_tsv(LOCAL_BLOCKS_TSV)
    multi = read_tsv(MULTI_TSV)

    mpo_domains = load_bed(HIC_ROOT / "genome_Mpo/genome_Mpo.TAD.40kb_domains.bed", ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "rgb"])
    mpo_score = load_bed(HIC_ROOT / "genome_Mpo/genome_Mpo.TAD.40kb_score.bedgraph", ["chrom", "start", "end", "score"])
    r108_domains = load_bed(HIC_ROOT / "genome_R108/genome_R108.TAD.40kb_domains.bed", ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "rgb"])

    fig = plt.figure(figsize=(15.2, 9.8), facecolor="white")
    outer = fig.add_gridspec(2, 2, height_ratios=[1.72, 1.08], width_ratios=[1.45, 1.05], wspace=0.24, hspace=0.30)

    left = outer[0, :].subgridspec(3, 1, height_ratios=[4.0, 0.85, 1.0], hspace=0.08)
    ax_heat = fig.add_subplot(left[0])
    ax_score = fig.add_subplot(left[1], sharex=ax_heat)
    ax_track = fig.add_subplot(left[2], sharex=ax_heat)

    plot_heatmap(
        ax_heat,
        "genome_Mpo",
        "Chr4",
        16_000_000,
        26_500_000,
        "a  Mpo Chr4 40-kb Hi-C/TAD context around the CEN5-associated remnant and active CEN4-like core",
        mpo_domains,
        [
            (18_316_439, 18_432_188, COL["cen5"], "CEN5-associated remnant"),
            (21_900_000, 24_650_000, COL["cen4"], "active CEN4-like core"),
        ],
    )
    ax_heat.tick_params(labelbottom=False)
    plot_score_track(ax_score, mpo_score, "Chr4", 16_000_000, 26_500_000, "TAD score")
    ax_score.tick_params(labelbottom=False)
    plot_linear_context(ax_track, context, local_blocks)

    rt = outer[1, 0].subgridspec(1, 2, wspace=0.18)
    ax_r4 = fig.add_subplot(rt[0])
    ax_r5 = fig.add_subplot(rt[1])
    plot_heatmap(
        ax_r4,
        "genome_R108",
        "Chr4",
        16_500_000,
        21_500_000,
        "b  R108 CEN4 reference TAD context",
        r108_domains,
        [(18_314_873, 19_993_588, COL["cen4"], "R108 CEN4")],
    )
    plot_heatmap(
        ax_r5,
        "genome_R108",
        "Chr5",
        22_000_000,
        27_500_000,
        "R108 CEN5 reference TAD context",
        r108_domains,
        [(23_929_531, 25_607_691, COL["cen5"], "R108 CEN5")],
    )
    ax_r4.set_xlabel("Chr4 (Mb)", fontsize=7)
    ax_r5.set_xlabel("Chr5 (Mb)", fontsize=7)

    right = outer[1, 1].subgridspec(2, 1, height_ratios=[1.05, 1.0], hspace=0.18)
    ax_fate = fig.add_subplot(right[0])
    ax_model = fig.add_subplot(right[1])
    plot_cen_fate(ax_fate, fate, multi)
    plot_model(ax_model)

    fig.text(
        0.015,
        0.012,
        "Core conclusion: In the current R108-renamed coordinate system, Mpo Chr4 carries a local CEN5-associated sequence remnant near a primary active CEN4-like CENH3 core, "
        "while CEN5 is absent from the seven primary active Mpo centromeres.",
        fontsize=7.5,
        color=COL["grey"],
    )
    outbase = OUTDIR / "Mpo_8to7_TAD_based_karyotype_model"
    fig.savefig(f"{outbase}.png", dpi=350, bbox_inches="tight")
    fig.savefig(f"{outbase}.pdf", bbox_inches="tight")
    fig.savefig(f"{outbase}.svg", bbox_inches="tight")
    fig.savefig(f"{outbase}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)
    summary = write_tad_summary(context, local_blocks, fate, multi)
    print(f"wrote {outbase}.png")
    print(f"wrote {outbase}.pdf")
    print(f"wrote {outbase}.svg")
    print(f"wrote {outbase}.tiff")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
