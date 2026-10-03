#!/usr/bin/env python3
"""Summarize local R108->Mpo sequence-block context around Mpo Chr4 CEN5 remnant."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def read_blocks(path: Path, min_len: int) -> list[dict[str, object]]:
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            q_len = int(float(row.get("query_len", 0)))
            if q_len < min_len:
                continue
            rows.append(
                {
                    "query_chr": row["query_chr"],
                    "query_start": int(float(row["query_start"])),
                    "query_end": int(float(row["query_end"])),
                    "target_chr": row["target_chr"],
                    "target_start": int(float(row["target_start"])),
                    "target_end": int(float(row["target_end"])),
                    "query_len": q_len,
                    "target_len": int(float(row.get("target_len", q_len))),
                    "strand": row.get("strand", "."),
                }
            )
    return rows


def ov(a: int, b: int, c: int, d: int) -> int:
    return max(0, min(b, d) - max(a, c))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    parser.add_argument("--chrom", default="Chr4")
    parser.add_argument("--start", type=int, default=16_000_000)
    parser.add_argument("--end", type=int, default=26_500_000)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--min-block-len", type=int, default=1000)
    args = parser.parse_args()

    root = Path(args.project_root)
    run = root / "output_key"
    out_dir = run / "05_fusion_chr_CEN_fate"
    fig_dir = run / "07_publication_figures"
    blocks_path = root / "output_synthen/02_cent_synteny/genome_R108_genome_Mpo/genome_R108_vs_genome_Mpo.rawPAF.blocks.tsv"
    blocks = read_blocks(blocks_path, args.min_block_len)

    local_blocks = [
        b
        for b in blocks
        if b["target_chr"] == args.chrom and ov(args.start, args.end, int(b["target_start"]), int(b["target_end"])) > 0
    ]
    local_blocks.sort(key=lambda x: (int(x["target_start"]), int(x["target_end"])))

    block_rows = []
    for b in local_blocks:
        block_rows.append(
            {
                "Mpo_chr": b["target_chr"],
                "Mpo_start": b["target_start"],
                "Mpo_end": b["target_end"],
                "R108_chr": b["query_chr"],
                "R108_start": b["query_start"],
                "R108_end": b["query_end"],
                "query_len": b["query_len"],
                "target_len": b["target_len"],
                "strand": b["strand"],
                "overlaps_CEN5_candidate": ov(int(b["target_start"]), int(b["target_end"]), 18_316_439, 18_432_188),
                "overlaps_active_Chr4_CEN": ov(int(b["target_start"]), int(b["target_end"]), 21_900_000, 24_650_000),
            }
        )

    bins = []
    for s in range(args.start, args.end, args.bin_size):
        e = min(args.end, s + args.bin_size)
        counts = Counter()
        length_by_chr = Counter()
        for b in local_blocks:
            x = ov(s, e, int(b["target_start"]), int(b["target_end"]))
            if x:
                counts[str(b["query_chr"])] += 1
                length_by_chr[str(b["query_chr"])] += x
        dominant_chr = ""
        dominant_bp = 0
        if length_by_chr:
            dominant_chr, dominant_bp = length_by_chr.most_common(1)[0]
        bins.append(
            {
                "Mpo_chr": args.chrom,
                "bin_start": s,
                "bin_end": e,
                "block_count": sum(counts.values()),
                "dominant_R108_chr_by_overlap_bp": dominant_chr,
                "dominant_R108_overlap_bp": dominant_bp,
                "R108_chr_overlap_bp": ";".join(f"{k}:{v}" for k, v in length_by_chr.most_common()),
                "R108_chr_block_count": ";".join(f"{k}:{v}" for k, v in counts.most_common()),
            }
        )

    summary_rows = []
    for qchr, xs in sorted(defaultdict(list, {k: [b for b in local_blocks if b["query_chr"] == k] for k in {str(b["query_chr"]) for b in local_blocks}}).items()):
        starts = [int(b["target_start"]) for b in xs]
        ends = [int(b["target_end"]) for b in xs]
        summary_rows.append(
            {
                "R108_chr": qchr,
                "Mpo_chr": args.chrom,
                "Mpo_span_start": min(starts),
                "Mpo_span_end": max(ends),
                "block_count": len(xs),
                "sum_query_len": sum(int(b["query_len"]) for b in xs),
                "sum_target_len": sum(int(b["target_len"]) for b in xs),
                "overlap_bp_with_CEN5_candidate": sum(ov(int(b["target_start"]), int(b["target_end"]), 18_316_439, 18_432_188) for b in xs),
                "overlap_bp_with_active_Chr4_CEN": sum(ov(int(b["target_start"]), int(b["target_end"]), 21_900_000, 24_650_000) for b in xs),
            }
        )
    summary_rows.sort(key=lambda x: int(x["sum_query_len"]), reverse=True)

    write_tsv(out_dir / "Mpo_Chr4_16_26Mb_R108_sequence_blocks.tsv", block_rows)
    write_tsv(out_dir / "Mpo_Chr4_16_26Mb_R108_sequence_blocks_100kb_windows.tsv", bins)
    write_tsv(out_dir / "Mpo_Chr4_16_26Mb_R108_sequence_blocks_by_chr_summary.tsv", summary_rows)

    primary_chrs = {f"Chr{i}" for i in range(1, 9)}
    primary_blocks = [b for b in local_blocks if str(b["query_chr"]) in primary_chrs]
    primary_summary = [row for row in summary_rows if row["R108_chr"] in primary_chrs]
    write_tsv(out_dir / "Mpo_Chr4_16_26Mb_R108_primary_chr_sequence_blocks_by_chr_summary.tsv", primary_summary)

    # Draft local synteny block plot.
    color_map = {
        "Chr4": "#2F80ED",
        "Chr5": "#D89B2B",
        "Chr3": "#7B61FF",
        "Chr6": "#27AE60",
        "Chr7": "#EB5757",
        "Chr8": "#9B51E0",
        "Chr1": "#56CCF2",
        "Chr2": "#F2994A",
    }
    y_order = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8"]
    y_pos = {chrom: i for i, chrom in enumerate(y_order)}
    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    for b in primary_blocks:
        qchr = str(b["query_chr"])
        y = y_pos.get(qchr, len(y_order))
        x0 = max(args.start, int(b["target_start"]))
        x1 = min(args.end, int(b["target_end"]))
        if x1 <= x0:
            continue
        ax.add_patch(
            Rectangle(
                ((x0 - args.start) / 1_000_000, y - 0.32),
                (x1 - x0) / 1_000_000,
                0.64,
                facecolor=color_map.get(qchr, "#999999"),
                edgecolor="none",
                alpha=0.85,
            )
        )
    ax.axvspan((18_316_439 - args.start) / 1_000_000, (18_432_188 - args.start) / 1_000_000, color="#D89B2B", alpha=0.25, label="CEN5 candidate")
    ax.axvspan((21_900_000 - args.start) / 1_000_000, (24_650_000 - args.start) / 1_000_000, color="#2F80ED", alpha=0.16, label="active Chr4 CEN")
    ax.set_yticks([y_pos[x] for x in y_order])
    ax.set_yticklabels([f"R108 {x}" for x in y_order])
    ax.set_xlim(0, (args.end - args.start) / 1_000_000)
    ax.set_xlabel(f"Mpo {args.chrom} position from {args.start / 1_000_000:.1f} Mb (Mb)")
    ax.set_title("R108 sequence blocks across the Mpo Chr4 CEN5-candidate region", fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    png = fig_dir / "Mpo_Chr4_16_26Mb_R108_local_synteny_blocks.png"
    pdf = fig_dir / "Mpo_Chr4_16_26Mb_R108_local_synteny_blocks.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")

    summary_md = out_dir / "Mpo_Chr4_16_26Mb_R108_local_synteny_context.summary.md"
    chr5 = next((x for x in summary_rows if x["R108_chr"] == "Chr5"), None)
    chr4 = next((x for x in summary_rows if x["R108_chr"] == "Chr4"), None)
    with summary_md.open("w") as handle:
        handle.write("# Mpo Chr4 Local Synteny Context Around CEN5 Candidate\n\n")
        handle.write(f"- Window: Mpo {args.chrom}:{args.start}-{args.end}\n")
        handle.write(f"- Total local blocks: {len(local_blocks)}\n")
        if chr5:
            handle.write(f"- R108 Chr5 blocks: {chr5['block_count']} blocks, {chr5['sum_query_len']} bp total query length, Mpo span {chr5['Mpo_span_start']}-{chr5['Mpo_span_end']}.\n")
        if chr4:
            handle.write(f"- R108 Chr4 blocks: {chr4['block_count']} blocks, {chr4['sum_query_len']} bp total query length, Mpo span {chr4['Mpo_span_start']}-{chr4['Mpo_span_end']}.\n")
        handle.write("- Interpretation: use this as local sequence-block context only. A concentrated R108 Chr5/CEN5-associated block cluster supports a local remnant, but not a broad genic Chr5 fusion segment unless flanking gene-level or chromosome-painting evidence is added.\n")
        handle.write(f"- Draft plot: `{png}`\n")

    print(out_dir / "Mpo_Chr4_16_26Mb_R108_sequence_blocks.tsv")
    print(out_dir / "Mpo_Chr4_16_26Mb_R108_sequence_blocks_by_chr_summary.tsv")
    print(png)
    print(summary_md)


if __name__ == "__main__":
    main()
