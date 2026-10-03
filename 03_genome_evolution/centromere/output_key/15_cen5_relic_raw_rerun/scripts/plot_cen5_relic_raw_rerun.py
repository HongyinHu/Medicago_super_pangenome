from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle


RUN = Path(__file__).resolve().parents[1]
RESULTS = RUN / "results"
FIGURES = RUN / "figures"
FIGURES.mkdir(exist_ok=True)

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
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "savefig.dpi": 600,
    }
)

COL = {
    "Chr3": "#D98B2B",
    "Chr5": "#4C78A8",
    "Chr6": "#59A14F",
    "base": "#D7DCE1",
    "transition": "#F2C7BA",
    "relic": "#B3343A",
    "cenh3": "#111111",
    "side": "#7B6D8D",
}


def mb(v):
    return float(v) / 1_000_000


def save_all(fig, prefix):
    for ext in ("png", "svg", "pdf", "tiff"):
        fig.savefig(FIGURES / f"{prefix}.{ext}", bbox_inches="tight")


def draw_chromosome(ax, chrom, length_mb, y, label):
    ax.add_patch(
        Rectangle(
            (0, y - 0.09),
            length_mb,
            0.18,
            facecolor=COL["base"],
            edgecolor="#9AA3AA",
            linewidth=0.6,
            zorder=1,
        )
    )
    ax.text(-2.0, y, label, va="center", ha="right", fontsize=8, fontweight="bold")
    ax.text(length_mb + 0.8, y, f"{length_mb:.0f} Mb", va="center", ha="left", color="#6B7178")


def draw_span(ax, start, end, y, color, label=None, height=0.26, alpha=1.0, z=3):
    ax.add_patch(
        Rectangle(
            (mb(start), y - height / 2),
            mb(end - start),
            height,
            facecolor=color,
            edgecolor="white",
            linewidth=0.6,
            alpha=alpha,
            zorder=z,
        )
    )
    if label:
        ax.text((mb(start) + mb(end)) / 2, y + height / 2 + 0.05, label, ha="center", va="bottom", fontsize=6)


def draw_relic(ax, start, end, y, label):
    x0, x1 = mb(start), mb(end)
    ax.add_patch(
        Rectangle(
            (x0, y - 0.28),
            x1 - x0,
            0.56,
            facecolor=COL["relic"],
            edgecolor=COL["relic"],
            linewidth=0.8,
            alpha=0.20,
            zorder=4,
        )
    )
    ax.plot([x0, x1], [y + 0.38, y + 0.38], color=COL["relic"], linewidth=1.2, zorder=5)
    ax.plot([x0, x0], [y + 0.30, y + 0.43], color=COL["relic"], linewidth=1.2, zorder=5)
    ax.plot([x1, x1], [y + 0.30, y + 0.43], color=COL["relic"], linewidth=1.2, zorder=5)
    ax.text((x0 + x1) / 2, y + 0.50, label, color=COL["relic"], ha="center", va="bottom", fontsize=6)


def draw_cenh3(ax, start, end, y, label):
    x0, x1 = mb(start), mb(end)
    xm = (x0 + x1) / 2
    ax.plot([xm, xm], [y - 0.33, y + 0.33], color=COL["cenh3"], linewidth=1.6, zorder=6)
    ax.scatter([xm], [y], s=95, color=COL["cenh3"], zorder=7)
    ax.text(xm, y - 0.48, label, ha="center", va="top", fontsize=6, color=COL["cenh3"])


def panel_a(ax):
    ax.set_title("a  Raw alignment places CEN5 relics at Chr5-to-Chr6 ancestry switches", loc="left", fontweight="bold")
    ax.set_xlim(-3, 62)
    ax.set_ylim(0.25, 2.70)
    ax.set_yticks([])
    ax.set_xlabel("Mpo chromosome position (Mb)")
    ax.set_xticks(np.arange(0, 61, 10))
    ax.grid(axis="x", color="#E7EAED", linewidth=0.5, zorder=0)

    draw_chromosome(ax, "Chr3", 54.5, 1.75, "Mpo Chr3")
    draw_chromosome(ax, "Chr5", 60.0, 0.85, "Mpo Chr5")

    spans = [
        (1.75, 18_011_115, 22_075_815, "Chr5", "R108 Chr5"),
        (1.75, 22_569_156, 23_054_983, "Chr6", "R108 Chr6"),
        (1.75, 38_119_173, 39_528_301, "Chr6", "R108 Chr6"),
        (1.75, 40_374_862, 45_959_693, "Chr3", "R108 Chr3"),
        (0.85, 18_501_230, 25_451_336, "Chr5", "R108 Chr5"),
        (0.85, 26_794_443, 28_650_988, "Chr6", "R108 Chr6"),
        (0.85, 54_789_758, 55_703_848, "Chr6", "R108 Chr6"),
        (0.85, 56_883_318, 59_331_231, "Chr3", "R108 Chr3"),
    ]
    for y, start, end, ancestral, label in spans:
        draw_span(ax, start, end, y, COL[ancestral])

    transitions = [
        (22_075_815, 22_569_156, 1.75),
        (39_528_301, 40_374_862, 1.75),
        (25_451_336, 26_794_443, 0.85),
        (55_703_848, 56_883_318, 0.85),
    ]
    for start, end, y in transitions:
        draw_span(ax, start, end, y, COL["transition"], height=0.44, alpha=0.65, z=2)

    draw_relic(ax, 20_591_033, 22_176_458, 1.75, "CEN5 relic")
    draw_relic(ax, 25_293_760, 25_645_386, 0.85, "weak CEN5 relic")
    draw_cenh3(ax, 46_300_000, 46_850_000, 1.75, "active CENH3")
    draw_cenh3(ax, 31_850_000, 32_250_000, 0.85, "active CENH3")

    legend_items = [
        ("R108 Chr5", COL["Chr5"]),
        ("R108 Chr6", COL["Chr6"]),
        ("R108 Chr3", COL["Chr3"]),
        ("transition", COL["transition"]),
        ("CEN5 relic", COL["relic"]),
        ("active CENH3", COL["cenh3"]),
    ]
    x0, y0 = 1.0, 2.48
    for i, (lab, color) in enumerate(legend_items):
        xx = x0 + i * 9.6
        ax.add_patch(Rectangle((xx, y0 - 0.05), 0.8, 0.10, facecolor=color, edgecolor="none", alpha=0.85))
        ax.text(xx + 1.0, y0, lab, va="center", ha="left", fontsize=6)


def panel_b(ax, focus):
    ax.set_title("b  CEN5-derived sequence homology", loc="left", fontweight="bold")
    selected = [
        ("Mpo Chr3\nroute relic", focus[(focus.target_chr == "Chr3") & (focus.start == 20591033)].iloc[0], COL["relic"]),
        ("Mpo Chr5\nweak route relic", focus[(focus.target_chr == "Chr5") & (focus.start == 25293760)].iloc[0], COL["relic"]),
        ("Mpo Chr4\nside signal", focus[(focus.target_chr == "Chr4") & (focus.start == 18314124)].iloc[0], COL["side"]),
    ]
    labels = [x[0] for x in selected]
    values = [x[1].aligned_bp / 1000 for x in selected]
    colors = [x[2] for x in selected]
    ax.barh(labels, values, color=colors, alpha=0.85)
    ax.invert_yaxis()
    ax.set_xlabel("Aligned CEN5-derived sequence (kb)")
    ax.grid(axis="x", color="#E7EAED", linewidth=0.5)
    for i, (_, row, _) in enumerate(selected):
        ax.text(row.aligned_bp / 1000 + 8, i, f"id {row.max_identity:.3f}", va="center", fontsize=6, color="#555555")


def panel_c(ax, focus, domains):
    ax.set_title("c  Relics lack active CENH3 enrichment", loc="left", fontweight="bold")
    active = domains[(domains.species == "genome_Mpo") & (domains.chrom.isin(["Chr3", "Chr5"]))]
    data = [
        ("Chr3\nrelic", float(focus[(focus.target_chr == "Chr3") & (focus.start == 20591033)].iloc[0].max_CENH3_log2), COL["relic"]),
        ("Chr5\nrelic", float(focus[(focus.target_chr == "Chr5") & (focus.start == 25293760)].iloc[0].max_CENH3_log2), COL["relic"]),
        ("Chr4\nside", float(focus[(focus.target_chr == "Chr4") & (focus.start == 18314124)].iloc[0].max_CENH3_log2), COL["side"]),
        ("Chr3\nactive", float(active[active.chrom == "Chr3"].iloc[0].max_smooth_log2), COL["cenh3"]),
        ("Chr5\nactive", float(active[active.chrom == "Chr5"].iloc[0].max_smooth_log2), COL["cenh3"]),
    ]
    labels = [d[0] for d in data]
    vals = [d[1] for d in data]
    colors = [d[2] for d in data]
    ax.bar(labels, vals, color=colors, alpha=0.85)
    ax.axhline(0, color="#555555", linewidth=0.7)
    ax.set_ylabel("Max raw CENH3 log2")
    ax.set_ylim(-0.25, 4.2)
    ax.grid(axis="y", color="#E7EAED", linewidth=0.5)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.08, f"{v:.2f}", ha="center", va="bottom", fontsize=6)


def panel_d(ax, kmer):
    ax.set_title("d  CEN5-specific exact 31-mers are depleted in Mpo", loc="left", fontweight="bold")
    wanted = [
        ("R108\nCEN5 core", "R108_CEN5_core", "Chr5:24000001-25700000", COL["Chr5"]),
        ("Mpo Chr3\nrelic", "Mpo_CEN5_relic_candidates", "Chr3:20591034-22176458", COL["relic"]),
        ("Mpo Chr5\nrelic", "Mpo_CEN5_relic_candidates", "Chr5:25293761-25645386", COL["relic"]),
        ("Mpo Chr3\nactive", "Mpo_active_CENH3_domains", "Chr3:46300001-46850000", COL["cenh3"]),
        ("Mpo Chr5\nactive", "Mpo_active_CENH3_domains", "Chr5:31850001-32250000", COL["cenh3"]),
    ]
    vals = []
    colors = []
    labels = []
    for lab, group, seq, color in wanted:
        row = kmer[(kmer.group == group) & (kmer.sequence == seq)].iloc[0]
        vals.append(float(row.top_CEN5_kmer_hits_per_kb))
        colors.append(color)
        labels.append(lab)
    ax.bar(labels, vals, color=colors, alpha=0.85)
    ax.set_ylabel("Top CEN5 31-mer hits per kb")
    ax.set_ylim(0, max(vals) * 1.18)
    ax.grid(axis="y", color="#E7EAED", linewidth=0.5)
    for i, v in enumerate(vals):
        label = f"{v:.1f}" if v else "0"
        ax.text(i, v + max(vals) * 0.025, label, ha="center", va="bottom", fontsize=6)


def main():
    focus = pd.read_csv(RESULTS / "Mpo_CEN5_relic_focus_candidates.tsv", sep="\t")
    domains = pd.read_csv(RESULTS / "raw_q20_called_CENH3_domains.tsv", sep="\t")
    kmer = pd.read_csv(RESULTS / "CEN5_specific_kmer_key_regions.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 6, height_ratios=[1.65, 1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :3])
    ax_c = fig.add_subplot(gs[1, 3:])
    ax_d = fig.add_subplot(gs[2, :])

    panel_a(ax_a)
    panel_b(ax_b, focus)
    panel_c(ax_c, focus, domains)
    panel_d(ax_d, kmer)

    save_all(fig, "Mpo_CEN5_relic_raw_rerun_evidence")


if __name__ == "__main__":
    main()
