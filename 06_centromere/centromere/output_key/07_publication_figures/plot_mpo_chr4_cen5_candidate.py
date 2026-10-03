#!/usr/bin/env python3
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "path/to/project/10.centromere_analysis"
RUN_ROOT = os.path.join(ROOT, "output_key")
SIGNAL = os.path.join(ROOT, "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.CENH3_vs_Input.q0.50k.log2.bedgraph")
EDTA = os.path.join(RUN_ROOT, "data/03_repeat/genome_Mpo.EDTA/genome_Mpo.fa.mod.EDTA.TEanno.gff3")
TRASH = os.path.join(RUN_ROOT, "data/03_repeat/genome_Mpo.TRASH/genome_Mpo.TRASH_arrays.sorted.bed")
GENES = os.path.join(ROOT, "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.gene.coord.tsv")
BLOCKS = os.path.join(RUN_ROOT, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN_sequence_block_projection_detail.tsv")
OUT_PNG = os.path.join(RUN_ROOT, "07_publication_figures/Mpo_Chr4_CEN5_candidate_track.png")
OUT_PDF = os.path.join(RUN_ROOT, "07_publication_figures/Mpo_Chr4_CEN5_candidate_track.pdf")

CHROM = "Chr4"
START = 16_500_000
END = 26_000_000
BIN = 100_000
CEN5_CLUSTER = (18_316_439, 18_432_188)
MINOR_CLUSTER = (17_900_000, 18_550_000)
ACTIVE_CEN = (21_900_000, 24_650_000)


def ov(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def read_signal():
    xs, ys = [], []
    with open(SIGNAL) as handle:
        for raw in handle:
            if not raw.strip():
                continue
            p = raw.rstrip("\n").split("\t")
            if p[0] != CHROM:
                continue
            s, e, val = int(float(p[1])), int(float(p[2])), float(p[3])
            if ov(START, END, s, e):
                xs.append((max(s, START) + min(e, END)) / 2 / 1e6)
                ys.append(val)
    return xs, ys


def intervals_from_gff(label_filter=None):
    items = []
    with open(EDTA) as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            p = raw.rstrip("\n").split("\t")
            if p[0] != CHROM:
                continue
            attr = p[8]
            if label_filter and label_filter not in attr:
                continue
            s, e = int(p[3]) - 1, int(p[4])
            if ov(START, END, s, e):
                items.append((s, e))
    return items


def intervals_from_bed(path):
    items = []
    with open(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if p[0] != CHROM:
                continue
            s, e = int(float(p[1])), int(float(p[2]))
            if ov(START, END, s, e):
                items.append((s, e))
    return items


def gene_intervals():
    items = []
    with open(GENES) as handle:
        for raw in handle:
            if not raw.strip():
                continue
            p = raw.rstrip("\n").split("\t")
            if p[1] != CHROM:
                continue
            s, e = int(p[2]), int(p[3])
            if ov(START, END, s, e):
                items.append((s, e))
    return items


def density(items):
    xs, ys = [], []
    for s in range(START, END, BIN):
        e = min(END, s + BIN)
        covered = sum(ov(s, e, a, b) for a, b in items)
        xs.append((s + e) / 2 / 1e6)
        ys.append(covered / (e - s))
    return xs, ys


def gene_density(items):
    xs, ys = [], []
    for s in range(START, END, BIN):
        e = min(END, s + BIN)
        n = sum(1 for a, b in items if ov(s, e, a, b))
        xs.append((s + e) / 2 / 1e6)
        ys.append(n / (BIN / 1e6))
    return xs, ys


def read_cen5_blocks():
    items = []
    with open(BLOCKS, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["ancestral_centromere_id"] != "CEN5" or row["target_chr"] != CHROM:
                continue
            s, e = int(row["target_start"]), int(row["target_end"])
            if ov(START, END, s, e):
                items.append((s, e, int(row["query_len"])))
    return items


def shade(ax, region, color, label="", alpha=0.2, y=0.96):
    ax.axvspan(region[0] / 1e6, region[1] / 1e6, color=color, alpha=alpha, lw=0)
    if label:
        ax.text((region[0] + region[1]) / 2 / 1e6, y, label, ha="center", va="top",
                transform=ax.get_xaxis_transform(), fontsize=8, color=color)


def main():
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    sig_x, sig_y = read_signal()
    repeat_x, repeat_y = density(intervals_from_gff())
    ltr_x, ltr_y = density(intervals_from_gff("classification=LTR"))
    trash_x, trash_y = density(intervals_from_bed(TRASH))
    gene_x, gene_y = gene_density(gene_intervals())
    blocks = read_cen5_blocks()

    fig, axes = plt.subplots(4, 1, figsize=(9, 7.2), sharex=True, gridspec_kw={"height_ratios": [1.3, 0.7, 1, 1]})

    ax = axes[0]
    ax.plot(sig_x, sig_y, color="#355C7D", lw=1.2)
    ax.set_ylabel("CENH3\nlog2")
    shade(ax, CEN5_CLUSTER, "#C44E52", "", 0.18)
    shade(ax, MINOR_CLUSTER, "#DD8452", "minor CENH3 signal", 0.12, y=0.98)
    shade(ax, ACTIVE_CEN, "#4C72B0", "primary active CEN", 0.16, y=0.94)
    ax.axhline(0, color="0.7", lw=0.8)

    ax = axes[1]
    for s, e, qlen in blocks:
        ax.add_patch(plt.Rectangle((s / 1e6, 0.25), max((e - s) / 1e6, 0.002), 0.5, color="#C44E52", alpha=0.75))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel("R108\nCEN5")
    shade(ax, CEN5_CLUSTER, "#C44E52", "CEN5 candidate block cluster", 0.08, y=0.98)

    ax = axes[2]
    ax.plot(repeat_x, repeat_y, color="#6A994E", lw=1.0, label="EDTA repeats")
    ax.plot(ltr_x, ltr_y, color="#9467BD", lw=1.0, label="LTR")
    ax.plot(trash_x, trash_y, color="#8C564B", lw=1.0, label="TRASH")
    ax.set_ylabel("Repeat\nfraction")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, ncol=3, fontsize=8, loc="upper right")

    ax = axes[3]
    ax.plot(gene_x, gene_y, color="#222222", lw=1.0)
    ax.set_ylabel("Genes\nper Mb")
    ax.set_xlabel("Mpo Chr4 position (Mb)")

    for ax in axes:
        ax.set_xlim(START / 1e6, END / 1e6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=8)

    fig.suptitle("Mpo Chr4 CEN5 candidate ancestral centromere interval", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
