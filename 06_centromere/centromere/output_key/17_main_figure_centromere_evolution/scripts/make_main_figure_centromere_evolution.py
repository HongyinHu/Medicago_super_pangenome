#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.colors import LinearSegmentedColormap


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


SPECIES = [
    ("genome_R108", "R108", 8),
    ("genome_A17", "A17", 8),
    ("genome_474", "474", 8),
    ("genome_Msa", "Msa", 8),
    ("genome_Mpo", "Mpo", 7),
]

PAL = {
    "chr": "#2f3640",
    "cen": "#151515",
    "cen5": "#2f8f83",
    "relic": "#40a486",
    "inactive": "#d95f5f",
    "active_red": "#d95f5f",
    "cen3": "#6f74b8",
    "cen6": "#d99058",
    "neutral": "#9a9a9a",
    "light": "#f2f2f2",
    "text": "#1f1f1f",
    "blue": "#4d83b5",
    "orange": "#c9895b",
    "gray": "#8a8a8a",
}


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fai_lengths(path):
    lengths = {}
    pat = re.compile(r"^Chr([1-8])$")
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            if pat.match(chrom):
                lengths[chrom] = int(parts[1])
    return lengths


def load_domains(run_root):
    rows = []
    rows.extend(read_tsv(run_root / "16_cen5_strong_validation/results/x8_raw_q20_called_CENH3_domains.tsv"))
    raw_rows = read_tsv(run_root / "15_cen5_relic_raw_rerun/results/raw_q20_called_CENH3_domains.tsv")
    rows.extend([r for r in raw_rows if r["species"] == "genome_Mpo"])
    by_species = defaultdict(dict)
    for r in rows:
        by_species[r["species"]][r["chrom"]] = {
            "start": int(float(r["start"])),
            "end": int(float(r["end"])),
            "max": float(r["max_smooth_log2"]),
            "mean": float(r["mean_log2"]),
            "length": int(float(r["end"])) - int(float(r["start"])),
        }
    return by_species


def load_signal(path, chroms):
    data = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for r in reader:
            if r["chrom"] in chroms:
                data[r["chrom"]].append((int(r["start"]), int(r["end"]), float(r["smooth_log2"])))
    return data


def save_pub(fig, prefix):
    fig.savefig(prefix + ".svg", bbox_inches="tight")
    fig.savefig(prefix + ".pdf", bbox_inches="tight")
    fig.savefig(prefix + ".png", dpi=450, bbox_inches="tight")
    fig.savefig(prefix + ".tiff", dpi=600, bbox_inches="tight")


def panel_a(ax, lengths, domains):
    ax.set_title("A  CENH3-defined active centromeres across complete Medicago assemblies",
                 loc="left", fontweight="bold")
    ax.set_xlim(0, 9.7)
    ax.set_ylim(-0.1, len(SPECIES) + 0.35)
    ax.axis("off")
    for row_idx, (species, label, nchr) in enumerate(SPECIES):
        y_base = len(SPECIES) - row_idx - 0.55
        ax.text(0.05, y_base + 0.28, f"{label}  x={nchr}", ha="left", va="center", fontweight="bold")
        species_max_len = max(lengths[species].values())
        for ci in range(1, nchr + 1):
            chrom = f"Chr{ci}"
            x = 1.35 + (ci - 1) * 0.95
            length = lengths[species].get(chrom, 0)
            h = 0.72 * length / species_max_len
            y0 = y_base - h / 2
            ax.add_patch(patches.FancyBboxPatch(
                (x - 0.10, y0), 0.20, h,
                boxstyle="round,pad=0,rounding_size=0.04",
                fc=PAL["chr"], ec="none", alpha=0.88))
            if chrom in domains[species] and length:
                cen = domains[species][chrom]
                mid = (cen["start"] + cen["end"]) / 2 / length
                cy = y0 + h * mid
                color = PAL["active_red"] if (species == "genome_R108" and chrom == "Chr5") else PAL["cen"]
                ax.add_patch(patches.Ellipse((x, cy), 0.28, 0.06, fc=color, ec="none", zorder=5))
            if row_idx == 0:
                ax.text(x, len(SPECIES) + 0.03, str(ci), ha="center", va="center", fontsize=6)
    ax.text(8.9, 4.55, "black: active CENH3\nred: x=8 CEN5", ha="left", va="top", fontsize=6)


def panel_b(ax, domains):
    species_order = [s for s, _, _ in SPECIES]
    labels = [l for _, l, _ in SPECIES]
    mat = np.full((len(species_order), 8), np.nan)
    for i, species in enumerate(species_order):
        for ci in range(1, 9):
            chrom = f"Chr{ci}"
            if chrom in domains[species]:
                mat[i, ci - 1] = domains[species][chrom]["length"] / 1_000_000
    cmap = LinearSegmentedColormap.from_list("cen_len", ["#eef4f2", "#80b9a8", "#236b5f"])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=np.nanmax(mat))
    ax.set_title("B  Active-centromere domain size", loc="left", fontweight="bold")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xticks(range(8))
    ax.set_xticklabels([f"Chr{i}" for i in range(1, 9)], rotation=45, ha="right")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.1f}", ha="center", va="center", fontsize=5.5, color="#111")
            else:
                ax.text(j, i, "-", ha="center", va="center", fontsize=6, color="#777")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("CENH3 domain (Mb)", fontsize=6)
    cb.ax.tick_params(labelsize=6)


def panel_c(ax, projection_rows):
    species_order = ["genome_R108", "genome_Msa", "genome_474", "genome_A17"]
    labels = ["R108", "Msa", "474", "A17"]
    sums = defaultdict(lambda: {"relic": 0.0, "active": 0.0})
    best_identity = defaultdict(float)
    for row in projection_rows:
        species = row["species"]
        sums[species]["relic"] += float(row["overlap_Mpo_CEN5_relic_bp"])
        sums[species]["active"] += float(row["overlap_Mpo_active_CENH3_bp"])
        best_identity[species] = max(best_identity[species], float(row["max_identity"]))
    x = np.arange(len(species_order))
    relic = np.array([sums[s]["relic"] / 1_000_000 for s in species_order])
    active = np.array([sums[s]["active"] / 1_000_000 for s in species_order])
    ax.bar(x, relic, color=PAL["relic"], width=0.58, label="Mpo CEN5 relic overlap")
    ax.bar(x, active, bottom=relic, color=PAL["inactive"], width=0.58, label="active CENH3 overlap")
    ax.set_title("C  x=8 CEN5 projects to relics, not active Mpo centromeres",
                 loc="left", fontweight="bold")
    ax.set_ylabel("Projected overlap (Mb)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, max(2.4, relic.max() + 0.45))
    for i, s in enumerate(species_order):
        ax.text(i, relic[i] + 0.07, f"active=0\nid {best_identity[s]:.2f}",
                ha="center", va="bottom", fontsize=5.8)
    ax.legend(loc="upper right", fontsize=6)


def draw_segmented_chr(ax, y, name, length_mb, segments, active=None, relics=None):
    x0, x1 = 1.2, 10.9
    ax.plot([x0, x1], [y, y], color="#333333", lw=5, solid_capstyle="round")
    ax.text(0.25, y, name, ha="left", va="center", fontweight="bold")
    for start, end, color, label in segments:
        xs = x0 + (start / length_mb) * (x1 - x0)
        xe = x0 + (end / length_mb) * (x1 - x0)
        ax.plot([xs, xe], [y, y], color=color, lw=7, solid_capstyle="butt")
        if label:
            ax.text((xs + xe) / 2, y + 0.18, label, ha="center", va="bottom", fontsize=6, color=color)
    if relics:
        for pos_mb, label in relics:
            x = x0 + (pos_mb / length_mb) * (x1 - x0)
            ax.add_patch(patches.Rectangle((x - 0.11, y - 0.17), 0.22, 0.34,
                                           fc=PAL["relic"], ec="none", alpha=0.95))
            ax.text(x, y - 0.27, label, ha="center", va="top", fontsize=6, color=PAL["relic"])
    if active:
        for pos_mb in active:
            x = x0 + (pos_mb / length_mb) * (x1 - x0)
            ax.add_patch(patches.Ellipse((x, y), 0.26, 0.22, fc=PAL["cen"], ec="none"))
            ax.text(x, y + 0.25, "active\nCENH3", ha="center", va="bottom", fontsize=5.5)


def panel_d(ax):
    ax.set_title("D  RTA-like Chr3/Chr5/Chr6 rearrangement leaves CEN5-derived relics in Mpo",
                 loc="left", fontweight="bold")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 2.45)
    ax.axis("off")
    draw_segmented_chr(
        ax, 1.55, "Mpo Chr3", 62,
        [(18.0, 22.2, PAL["cen5"], "ancestral Chr5"),
         (22.6, 39.5, PAL["cen6"], "ancestral Chr6"),
         (40.4, 50.6, PAL["cen3"], "ancestral Chr3")],
        active=[46.6],
        relics=[(21.4, "CEN5 relic")]
    )
    draw_segmented_chr(
        ax, 0.65, "Mpo Chr5", 60,
        [(25.3, 25.7, PAL["cen5"], ""),
         (26.8, 55.7, PAL["cen6"], "ancestral Chr6"),
         (56.9, 60.0, PAL["cen3"], "ancestral Chr3")],
        active=[32.0],
        relics=[(25.45, "CEN5 relic")]
    )
    ax.annotate("ancestral CEN5 not retained\nas an active centromere",
                xy=(5.0, 1.08), xytext=(6.9, 2.18),
                arrowprops=dict(arrowstyle="->", color=PAL["inactive"], lw=1.2),
                ha="center", va="center", color=PAL["inactive"], fontsize=7)
    legend_items = [
        patches.Patch(color=PAL["cen5"], label="ancestral Chr5/CEN5 block"),
        patches.Patch(color=PAL["cen6"], label="ancestral Chr6 block"),
        patches.Patch(color=PAL["cen3"], label="ancestral Chr3 block"),
        patches.Patch(color=PAL["relic"], label="CEN5-derived relic"),
        patches.Ellipse((0, 0), 0.3, 0.18, color=PAL["cen"], label="current active CENH3"),
    ]
    ax.legend(handles=legend_items, loc="lower right", fontsize=6, ncol=2, handlelength=1.3)


def panel_e(ax, signal):
    ax.set_title("E  CEN5 relics lack CENH3 enrichment", loc="left", fontweight="bold")
    views = [
        ("Chr3", 18_000_000, 48_500_000, [(20_591_033, 22_176_458, "relic")], [(46_300_000, 46_850_000, "active")]),
        ("Chr5", 24_500_000, 33_500_000, [(25_293_760, 25_645_386, "relic")], [(31_850_000, 32_250_000, "active")]),
    ]
    offsets = [2.5, 0.0]
    for idx, (chrom, start, end, relics, actives) in enumerate(views):
        offset = offsets[idx]
        xs = []
        ys = []
        for s, e, val in signal[chrom]:
            mid = (s + e) / 2
            if start <= mid <= end:
                xs.append((mid - start) / 1_000_000)
                ys.append(val + offset)
        ax.plot(xs, ys, color="#355f53", lw=0.9)
        ax.axhline(offset, color="#cfcfcf", lw=0.6)
        ax.text(-0.45, offset, chrom, ha="right", va="center", fontweight="bold")
        for s, e, label in relics:
            x0 = (s - start) / 1_000_000
            x1 = (e - start) / 1_000_000
            ax.add_patch(patches.Rectangle((x0, offset - 0.25), x1 - x0, 0.5,
                                           fc=PAL["relic"], ec="none", alpha=0.25))
            ax.text((x0 + x1) / 2, offset - 0.55, "CEN5 relic", ha="center", va="top", fontsize=6, color=PAL["relic"])
        for s, e, label in actives:
            x0 = (s - start) / 1_000_000
            x1 = (e - start) / 1_000_000
            ax.add_patch(patches.Rectangle((x0, offset - 0.25), x1 - x0, 0.5,
                                           fc=PAL["cen"], ec="none", alpha=0.18))
            ax.text((x0 + x1) / 2, offset + 1.02, "active CENH3", ha="center", va="bottom", fontsize=6)
    ax.set_xlim(0, 30.5)
    ax.set_ylim(-1.1, 5.3)
    ax.set_xlabel("Local coordinate within shown interval (Mb)")
    ax.text(-1.55, 2.15, "CENH3/Input log2", rotation=90, ha="center", va="center")
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])


def panel_f(ax, repeat_rows):
    keep = ["R108_CEN5_core", "Mpo_Chr3_CEN5_relic", "Mpo_Chr5_CEN5_relic"]
    labels = ["R108\nCEN5", "Mpo Chr3\nrelic", "Mpo Chr5\nrelic"]
    by_id = {r["target_id"]: r for r in repeat_rows}
    trash = [float(by_id[k]["trash_coverage_pct"]) for k in keep]
    trf = [float(by_id[k]["trf_coverage_pct"]) for k in keep]
    x = np.arange(len(keep))
    w = 0.32
    ax.bar(x - w / 2, trash, width=w, color=PAL["blue"], label="CEN5 TRASH hit")
    ax.bar(x + w / 2, trf, width=w, color=PAL["orange"], label="TRF tandem array")
    ax.set_title("F  CEN5 repeat signatures are eroded in relics", loc="left", fontweight="bold")
    ax.set_ylim(0, 105)
    ax.set_ylabel("Coverage (%)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    for i, v in enumerate(trash):
        ax.text(i - w / 2, v + 2, f"{v:.1f}", ha="center", va="bottom", fontsize=6)
    for i, v in enumerate(trf):
        ax.text(i + w / 2, v + 2, f"{v:.1f}", ha="center", va="bottom", fontsize=6)
    ax.legend(loc="upper right", fontsize=6)


def write_caption(out_path):
    text = """# Draft main figure layout

Title: Complete assemblies reveal dynamic centromere evolution in Medicago

Core message: complete assemblies and raw CENH3 re-mapping define active centromeres across Medicago, reveal variation in CENH3 domain size, and show that ancestral CEN5 was not retained as an active Mpo centromere during the x=8 to x=7 karyotype transition. Instead, CEN5-derived relics remain near Chr5-to-Chr6 transition regions on Mpo Chr3 and Chr5, lack CENH3 enrichment, and show eroded centromeric repeat signatures.

Panel guide:
- A: CENH3-defined active centromere atlas across R108, A17, 474, Msa and Mpo.
- B: Active-centromere domain size heatmap.
- C: x=8 CEN5 projections to Mpo; active CENH3 overlap is zero.
- D: Mpo Chr3/Chr5 RTA-like rearrangement model with CEN5 relics.
- E: Mpo Chr3/Chr5 CENH3/Input tracks around relic and active-centromere intervals.
- F: Targeted repeat decay of R108 CEN5 core versus Mpo Chr3/Chr5 relics.

Keep out of the main figure for now: nucleotide-level junction repeat-gap intervals and the Chr4 CEN5-like side signal. These are better as Extended Data/Supplementary controls.
"""
    out_path.write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    run_root = Path(args.run_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lengths = {
        s: read_fai_lengths(run_root / f"data/00_genome/{s}.fa.fai")
        for s, _, _ in SPECIES
    }
    domains = load_domains(run_root)
    projection_rows = read_tsv(run_root / "16_cen5_strong_validation/results/x8_CEN5_pm3Mb_projection_to_Mpo_summary.tsv")
    repeat_rows = read_tsv(run_root / "16_cen5_strong_validation/results/targeted_repeat_decay_summary.tsv")
    signal = load_signal(run_root / "15_cen5_relic_raw_rerun/results/signal/genome_Mpo.q20.50k.log2.tsv",
                         {"Chr3", "Chr5"})

    fig = plt.figure(figsize=(8.2, 10.4), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.045, h_pad=0.045, wspace=0.11, hspace=0.12)
    gs = fig.add_gridspec(4, 2, height_ratios=[1.42, 1.05, 1.25, 1.2], width_ratios=[1.05, 1.0])

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])
    ax_e = fig.add_subplot(gs[3, 0])
    ax_f = fig.add_subplot(gs[3, 1])

    panel_a(ax_a, lengths, domains)
    panel_b(ax_b, domains)
    panel_c(ax_c, projection_rows)
    panel_d(ax_d)
    panel_e(ax_e, signal)
    panel_f(ax_f, repeat_rows)

    prefix = str(outdir / "draft_main_figure_centromere_evolution")
    save_pub(fig, prefix)
    write_caption(outdir / "draft_main_figure_centromere_evolution_notes.md")


if __name__ == "__main__":
    main()
