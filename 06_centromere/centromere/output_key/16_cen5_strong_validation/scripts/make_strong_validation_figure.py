#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse


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


PALETTE = {
    "cen5": "#2f8f83",
    "cen6": "#d99058",
    "cen3": "#6f74b8",
    "gap": "#d65f5f",
    "active": "#1f1f1f",
    "trash": "#4d83b5",
    "trf": "#c9895b",
    "edta": "#8a8a8a",
    "relic": "#3f9d7a",
    "inactive": "#d95f5f",
    "neutral": "#b7b7b7",
}


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_pub(fig, prefix):
    fig.savefig(prefix + ".svg", bbox_inches="tight")
    fig.savefig(prefix + ".pdf", bbox_inches="tight")
    fig.savefig(prefix + ".png", dpi=450, bbox_inches="tight")
    fig.savefig(prefix + ".tiff", dpi=600, bbox_inches="tight")


def mb(value):
    return float(value) / 1_000_000.0


def draw_route_panel(ax, junction_rows):
    chrom_lengths = {"Chr3": 62.0, "Chr5": 60.0}
    y_positions = {"Chr3": 1.25, "Chr5": 0.35}
    relics = {
        "Chr3": [(20_591_033, 22_176_458, "CEN5 relic")],
        "Chr5": [(25_293_760, 25_645_386, "CEN5 relic")],
    }
    active = {
        "Chr3": [(46_300_000, 46_850_000)],
        "Chr5": [(31_850_000, 32_250_000)],
    }

    ax.set_xlim(0, 62)
    ax.set_ylim(0, 1.8)
    ax.set_yticks([])
    ax.set_xlabel("Mpo coordinate (Mb)")
    ax.set_title("A  CEN5 relics occur near ancestral-block transitions on Mpo Chr3/Chr5",
                 loc="left", fontweight="bold")

    for chrom, length in chrom_lengths.items():
        y = y_positions[chrom]
        ax.plot([0, length], [y, y], color="#333333", lw=3, solid_capstyle="round")
        ax.text(-1.8, y, chrom, va="center", ha="right", fontweight="bold")

    for row in junction_rows:
        chrom = row["Mpo_chr"]
        if chrom not in y_positions:
            continue
        y = y_positions[chrom]
        start = mb(row["refined_interval_start"])
        end = mb(row["refined_interval_end"])
        ax.add_patch(Rectangle((start, y - 0.17), end - start, 0.34,
                               fill=False, ec=PALETTE["gap"], lw=1.3, ls="--"))
        label = row["transition"].replace("_", " ").replace("CEN", "CEN")
        ax.text((start + end) / 2, y + 0.26, f"{label}\n{end - start:.2f} Mb gap",
                color=PALETTE["gap"], ha="center", va="bottom", fontsize=6)

    for chrom, intervals in relics.items():
        y = y_positions[chrom]
        for start, end, label in intervals:
            ax.add_patch(Rectangle((mb(start), y - 0.10), mb(end - start), 0.20,
                                   color=PALETTE["cen5"], alpha=0.85, lw=0))
            ax.text((mb(start) + mb(end)) / 2, y - 0.28, label, color=PALETTE["cen5"],
                    ha="center", va="top", fontsize=6)

    for chrom, intervals in active.items():
        y = y_positions[chrom]
        for start, end in intervals:
            ax.add_patch(Ellipse(((mb(start) + mb(end)) / 2, y), width=0.85, height=0.22,
                                 color=PALETTE["active"]))
            ax.text((mb(start) + mb(end)) / 2, y + 0.20, "active CENH3",
                    color=PALETTE["active"], ha="center", va="bottom", fontsize=6)

    handles = [
        Rectangle((0, 0), 1, 1, color=PALETTE["cen5"], label="CEN5-derived relic"),
        Rectangle((0, 0), 1, 1, fill=False, ec=PALETTE["gap"], ls="--", label="repeat-gap transition"),
        Ellipse((0, 0), 1, 0.3, color=PALETTE["active"], label="current active CENH3"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=6, handlelength=1.4)


def draw_repeat_panel(ax, repeat_rows):
    order = [
        "R108_CEN5_core",
        "Mpo_Chr3_CEN5_relic",
        "Mpo_Chr5_CEN5_relic",
        "Mpo_Chr4_CEN5_side_signal",
        "Mpo_Chr3_active_CENH3",
        "Mpo_Chr5_active_CENH3",
    ]
    labels = {
        "R108_CEN5_core": "R108\nCEN5",
        "Mpo_Chr3_CEN5_relic": "Chr3\nrelic",
        "Mpo_Chr5_CEN5_relic": "Chr5\nrelic",
        "Mpo_Chr4_CEN5_side_signal": "Chr4\nside",
        "Mpo_Chr3_active_CENH3": "Chr3\nactive",
        "Mpo_Chr5_active_CENH3": "Chr5\nactive",
    }
    by_id = {row["target_id"]: row for row in repeat_rows}
    x = list(range(len(order)))
    width = 0.25
    trash = [float(by_id[k]["trash_coverage_pct"]) for k in order]
    trf = [float(by_id[k]["trf_coverage_pct"]) for k in order]
    edta = [float(by_id[k]["edta_pct"]) for k in order]
    ax.bar([i - width for i in x], trash, width=width, color=PALETTE["trash"], label="CEN5 TRASH hit")
    ax.bar(x, trf, width=width, color=PALETTE["trf"], label="TRF tandem array")
    ax.bar([i + width for i in x], edta, width=width, color=PALETTE["edta"], label="EDTA TE")
    ax.set_ylim(0, 105)
    ax.set_ylabel("Coverage (%)")
    ax.set_xticks(x)
    ax.set_xticklabels([labels[k] for k in order], rotation=0)
    ax.set_title("B  CEN5 repeat decay", loc="left", fontweight="bold")
    ax.legend(loc="upper right", fontsize=6, ncol=1)
    for i, value in enumerate(trash):
        ax.text(i - width, value + 2.2, f"{value:.1f}", ha="center", va="bottom", fontsize=5)


def draw_projection_panel(ax, projection_rows):
    species_order = ["genome_R108", "genome_Msa", "genome_474", "genome_A17"]
    labels = ["R108", "Msa", "474", "A17"]
    sums = defaultdict(lambda: {"relic": 0.0, "active": 0.0})
    best_identity = defaultdict(float)
    best_aligned = defaultdict(float)
    for row in projection_rows:
        species = row["species"]
        sums[species]["relic"] += float(row["overlap_Mpo_CEN5_relic_bp"])
        sums[species]["active"] += float(row["overlap_Mpo_active_CENH3_bp"])
        best_identity[species] = max(best_identity[species], float(row["max_identity"]))
        best_aligned[species] = max(best_aligned[species], float(row["aligned_bp"]))

    x = list(range(len(species_order)))
    relic_mb = [sums[s]["relic"] / 1_000_000.0 for s in species_order]
    active_mb = [sums[s]["active"] / 1_000_000.0 for s in species_order]
    ax.bar(x, relic_mb, color=PALETTE["relic"], width=0.55, label="overlap with Mpo CEN5 relic")
    ax.bar(x, active_mb, bottom=relic_mb, color=PALETTE["inactive"], width=0.55, label="overlap with active CENH3")
    ax.axhline(0, color="#333333", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Projected overlap (Mb)")
    ax.set_ylim(0, max(relic_mb + active_mb) + 0.55)
    ax.set_title("C  x=8 CEN5 avoids active Mpo CENH3", loc="left", fontweight="bold")
    ax.legend(loc="upper right", fontsize=6)
    for i, species in enumerate(species_order):
        ax.text(i, relic_mb[i] + 0.08, f"active=0\nmax id {best_identity[species]:.2f}",
                ha="center", va="bottom", fontsize=5.8, color="#222222")


def draw_junction_panel(ax, junction_rows):
    labels = []
    widths = []
    colors = []
    label_map = {
        "Chr3_CEN5_to_CEN6": "Chr3 CEN5->CEN6",
        "Chr3_CEN6_to_CEN3": "Chr3 CEN6->CEN3",
        "Chr5_CEN5_to_CEN6": "Chr5 CEN5->CEN6",
        "Chr5_CEN6_to_CEN3": "Chr5 CEN6->CEN3",
    }
    for row in junction_rows:
        labels.append(label_map.get(row["transition"], row["transition"]))
        widths.append(float(row["refined_interval_bp"]) / 1_000_000.0)
        colors.append(PALETTE["gap"] if "CEN5" in row["transition"] else "#8d8d8d")
    y = list(range(len(labels)))
    ax.barh(y, widths, color=colors, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Unresolved repeat-gap interval (Mb)")
    ax.set_title("D  Junctions remain repeat-gap intervals", loc="left", fontweight="bold")
    for yi, w in zip(y, widths):
        ax.text(w + 0.03, yi, f"{w:.2f} Mb", va="center", fontsize=6)
    ax.set_xlim(0, max(widths) * 1.25)


def write_report(base, repeat_rows, projection_rows, junction_rows, out_path):
    by_id = {row["target_id"]: row for row in repeat_rows}
    proj_active = sum(float(row["overlap_Mpo_active_CENH3_bp"]) for row in projection_rows)
    relic_by_species = defaultdict(float)
    for row in projection_rows:
        relic_by_species[row["species"]] += float(row["overlap_Mpo_CEN5_relic_bp"])
    cen5_gaps = [row for row in junction_rows if "CEN5" in row["transition"]]

    lines = [
        "# Mpo CEN5 relic 强验证结果小结",
        "",
        "## 当前结论",
        "",
        "这轮从原始数据重新补做的三类证据支持一个较稳妥的表述：Mpo 的 Chr3/Chr5 上存在 ancestral CEN5-derived relic，且这些 relic 位于 Chr5-to-Chr6 祖先片段转换附近；这些区域没有形成当前活性 CENH3 peak，CEN5 相关重复序列也呈明显碎片化/退化。因此，目前更适合写作“ancestral CEN5 was not retained as an active centromere during the RTA-like/complex rearrangement”，暂时不写成“已单碱基证明 CEN5 原位失活”。",
        "",
        "## 关键数值",
        "",
        f"- R108 CEN5 core 的 TRASH/CEN5 repeat 覆盖率为 {float(by_id['R108_CEN5_core']['trash_coverage_pct']):.2f}%，TRF tandem array 覆盖率为 {float(by_id['R108_CEN5_core']['trf_coverage_pct']):.2f}%。",
        f"- Mpo Chr3 CEN5 relic 的 TRASH 覆盖率为 {float(by_id['Mpo_Chr3_CEN5_relic']['trash_coverage_pct']):.2f}%，TRF 覆盖率为 {float(by_id['Mpo_Chr3_CEN5_relic']['trf_coverage_pct']):.2f}%。",
        f"- Mpo Chr5 CEN5 relic 的 TRASH 覆盖率为 {float(by_id['Mpo_Chr5_CEN5_relic']['trash_coverage_pct']):.2f}%，TRF 覆盖率为 {float(by_id['Mpo_Chr5_CEN5_relic']['trf_coverage_pct']):.2f}%。",
        f"- x=8 CEN5 +/-3 Mb 投影到 Mpo 当前活性 CENH3 区域的总重叠为 {proj_active:.0f} bp。",
    ]
    for species in ["genome_R108", "genome_Msa", "genome_474", "genome_A17"]:
        lines.append(f"- {species} CEN5 +/-3 Mb 投影与 Mpo CEN5 relic 的总重叠为 {relic_by_species[species]:.0f} bp。")
    for row in cen5_gaps:
        lines.append(
            f"- {row['transition']} refined interval: {row['Mpo_chr']}:{row['refined_interval_start']}-{row['refined_interval_end']} "
            f"({float(row['refined_interval_bp']) / 1_000_000:.2f} Mb)，状态为 {row['status']}。"
        )
    lines.extend([
        "",
        "## 写作边界",
        "",
        "- 可以作为主文/补充图支持：CEN5 relic、repeat decay、无活性 CENH3 重叠、RTA-like 3/5/6 重排路线。",
        "- 不能强写为：已经找到 nucleotide-level fusion junction 或已经直接证明单个位点的 CEN5 原位失活。",
        "- 若要再升级，需要继续做长读长/局部组装级 junction 验证，或针对 repeat-gap 区域设计 PCR/HiFi-read spanning 证据。",
        "",
        f"结果目录：{base}",
    ])
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, help="16_cen5_strong_validation directory")
    args = parser.parse_args()
    base = Path(args.base)
    result_dir = base / "results"
    fig_dir = base / "figures"
    fig_dir.mkdir(exist_ok=True)

    junction_rows = read_tsv(result_dir / "junction_refined_intervals.tsv")
    repeat_rows = read_tsv(result_dir / "targeted_repeat_decay_summary.tsv")
    projection_rows = read_tsv(result_dir / "x8_CEN5_pm3Mb_projection_to_Mpo_summary.tsv")

    fig = plt.figure(figsize=(7.2, 7.6), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.04, h_pad=0.04, wspace=0.12, hspace=0.10)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.0, 0.72], width_ratios=[1.22, 1.0])
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, :])

    draw_route_panel(ax_a, junction_rows)
    draw_repeat_panel(ax_b, repeat_rows)
    draw_projection_panel(ax_c, projection_rows)
    draw_junction_panel(ax_d, junction_rows)

    prefix = str(fig_dir / "Mpo_CEN5_strong_validation_summary")
    save_pub(fig, prefix)
    write_report(base, repeat_rows, projection_rows, junction_rows,
                 base / "Mpo_CEN5_strong_validation_summary_CN.md")


if __name__ == "__main__":
    main()
