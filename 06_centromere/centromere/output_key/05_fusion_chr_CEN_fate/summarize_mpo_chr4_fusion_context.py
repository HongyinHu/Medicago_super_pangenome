#!/usr/bin/env python3
"""Summarize the Mpo Chr4 CEN4/CEN5 centromere-fate context."""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_chrom_sizes(fai: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with fai.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, size, *_ = line.rstrip("\n").split("\t")
            sizes[chrom] = int(size)
    return sizes


def interval_distance(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def summarize_gene_links(gene_links: list[dict[str, str]], target_chr: str) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for row in gene_links:
        if row["sm_chr"] == target_chr:
            grouped[row["sa_chr"]].append(row)

    out: list[dict[str, object]] = []
    for sa_chr, rows in sorted(grouped.items()):
        sm_starts = [int(x["sm_start"]) for x in rows]
        sm_ends = [int(x["sm_end"]) for x in rows]
        sa_starts = [int(x["sa_start"]) for x in rows]
        sa_ends = [int(x["sa_end"]) for x in rows]
        out.append(
            {
                "R108_chr": sa_chr,
                "Mpo_chr": target_chr,
                "anchor_count": len(rows),
                "Mpo_min": min(sm_starts),
                "Mpo_max": max(sm_ends),
                "Mpo_span_bp": max(sm_ends) - min(sm_starts),
                "R108_min": min(sa_starts),
                "R108_max": max(sa_ends),
                "R108_span_bp": max(sa_ends) - min(sa_starts),
            }
        )
    out.sort(key=lambda x: int(x["anchor_count"]), reverse=True)
    return out


def summarize_windows(gene_links: list[dict[str, str]], target_chr: str, chr_len: int, bin_size: int) -> list[dict[str, object]]:
    bins: list[collections.Counter[str]] = [collections.Counter() for _ in range((chr_len + bin_size - 1) // bin_size)]
    for row in gene_links:
        if row["sm_chr"] != target_chr:
            continue
        midpoint = (int(row["sm_start"]) + int(row["sm_end"])) // 2
        idx = min(midpoint // bin_size, len(bins) - 1)
        bins[idx][row["sa_chr"]] += 1

    out: list[dict[str, object]] = []
    for idx, counts in enumerate(bins):
        start = idx * bin_size
        end = min(chr_len, start + bin_size)
        if counts:
            dominant_chr, dominant_count = counts.most_common(1)[0]
            counts_string = ";".join(f"{chrom}:{count}" for chrom, count in counts.most_common())
        else:
            dominant_chr, dominant_count, counts_string = "NA", 0, ""
        out.append(
            {
                "Mpo_chr": target_chr,
                "bin_start": start,
                "bin_end": end,
                "total_anchor_count": sum(counts.values()),
                "dominant_R108_chr": dominant_chr,
                "dominant_anchor_count": dominant_count,
                "R108_chr_counts": counts_string,
            }
        )
    return out


def summarize_clusters(cluster_rows: list[dict[str, str]], target_chr: str, cids: set[str]) -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    for cid in cids:
        rows = [x for x in cluster_rows if x["ancestral_centromere_id"] == cid and x["Mpo_chr"] == target_chr]
        if not rows:
            continue
        rows.sort(key=lambda x: int(x["sum_query_len"]), reverse=True)
        starts = [int(x["cluster_start"]) for x in rows]
        ends = [int(x["cluster_end"]) for x in rows]
        out[cid] = {
            "cluster_count": len(rows),
            "total_sum_query_len": sum(int(x["sum_query_len"]) for x in rows),
            "total_block_count": sum(int(x["block_count"]) for x in rows),
            "span_start": min(starts),
            "span_end": max(ends),
            "top_start": int(rows[0]["cluster_start"]),
            "top_end": int(rows[0]["cluster_end"]),
            "top_block_count": int(rows[0]["block_count"]),
            "top_sum_query_len": int(rows[0]["sum_query_len"]),
            "top_R108_chr": rows[0]["R108_query_chr"],
            "top_R108_min": int(rows[0]["R108_query_min"]),
            "top_R108_max": int(rows[0]["R108_query_max"]),
        }
    return out


def write_dict_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    parser.add_argument("--target-chr", default="Chr4")
    parser.add_argument("--bin-size", type=int, default=1_000_000)
    args = parser.parse_args()

    project_root = Path(args.project_root)
    out_dir = project_root / "output_key/05_fusion_chr_CEN_fate"
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_links = read_tsv(project_root / "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_R108.genome_Mpo.gene_links.tsv")
    clusters = read_tsv(project_root / "output_key/05_fusion_chr_CEN_fate/R108_CEN_projection_block_clusters.tsv")
    active_back = read_tsv(project_root / "output_key/04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_R108.tsv")
    chrom_sizes = read_chrom_sizes(project_root / "output_all/00_genome/genome_Mpo.fa.fai")

    gene_summary = summarize_gene_links(gene_links, args.target_chr)
    window_summary = summarize_windows(gene_links, args.target_chr, chrom_sizes[args.target_chr], args.bin_size)
    cluster_summary = summarize_clusters(clusters, args.target_chr, {"CEN4", "CEN5"})

    write_dict_rows(out_dir / "Mpo_Chr4_R108_gene_anchor_distribution.tsv", gene_summary)
    write_dict_rows(out_dir / "Mpo_Chr4_R108_gene_anchor_1Mb_windows.tsv", window_summary)

    active_chr4 = next(x for x in active_back if x["Mpo_chr"] == args.target_chr)
    active_start = int(active_chr4["Mpo_start"])
    active_end = int(active_chr4["Mpo_end"])
    cen5 = cluster_summary.get("CEN5")
    cen4 = cluster_summary.get("CEN4")
    if not cen4 or not cen5:
        raise ValueError("Missing CEN4 or CEN5 cluster summary on target chromosome")

    nearest_cen4_cluster_distance = min(
        interval_distance(cen5["top_start"], cen5["top_end"], int(row["cluster_start"]), int(row["cluster_end"]))
        for row in clusters
        if row["ancestral_centromere_id"] == "CEN4" and row["Mpo_chr"] == args.target_chr
    )

    context_rows = [
        {
            "feature": "active_Mpo_Chr4_CENH3_core",
            "Mpo_chr": args.target_chr,
            "start": active_start,
            "end": active_end,
            "nearest_R108_CEN": active_chr4["nearest_R108_CEN"],
            "classification": active_chr4["classification"],
            "note": "Primary active CENH3 core on Mpo Chr4 reciprocally projects to R108 CEN4.",
        },
        {
            "feature": "R108_CEN4_projection_span_on_Mpo_Chr4",
            "Mpo_chr": args.target_chr,
            "start": cen4["span_start"],
            "end": cen4["span_end"],
            "nearest_R108_CEN": "CEN4",
            "classification": "ancestral_CEN4_projection",
            "note": f"{cen4['cluster_count']} clusters; total projected query length {cen4['total_sum_query_len']} bp.",
        },
        {
            "feature": "R108_CEN5_top_projection_cluster_on_Mpo_Chr4",
            "Mpo_chr": args.target_chr,
            "start": cen5["top_start"],
            "end": cen5["top_end"],
            "nearest_R108_CEN": "CEN5",
            "classification": "ancestral_CEN5_candidate_remnant",
            "note": f"{cen5['top_block_count']} blocks; projected R108 interval {cen5['top_R108_chr']}:{cen5['top_R108_min']}-{cen5['top_R108_max']}.",
        },
        {
            "feature": "distance_CEN5_top_cluster_to_active_Chr4_CENH3_core",
            "Mpo_chr": args.target_chr,
            "start": cen5["top_end"],
            "end": active_start,
            "nearest_R108_CEN": "CEN5_to_active_CEN4",
            "classification": "distance_bp",
            "note": str(interval_distance(cen5["top_start"], cen5["top_end"], active_start, active_end)),
        },
        {
            "feature": "distance_CEN5_top_cluster_to_nearest_CEN4_projection_cluster",
            "Mpo_chr": args.target_chr,
            "start": cen5["top_end"],
            "end": cen4["span_start"],
            "nearest_R108_CEN": "CEN5_to_CEN4_projection",
            "classification": "distance_bp",
            "note": str(nearest_cen4_cluster_distance),
        },
    ]
    write_dict_rows(out_dir / "Mpo_Chr4_CEN4_CEN5_fusion_context.tsv", context_rows)

    top_gene_sources = ", ".join(f"{row['R108_chr']} ({row['anchor_count']} anchors)" for row in gene_summary[:5])
    summary_path = out_dir / "Mpo_Chr4_CEN4_CEN5_fusion_context.summary.md"
    with summary_path.open("w") as handle:
        handle.write("# Mpo Chr4 CEN4/CEN5 Fusion Context\n\n")
        handle.write(f"- Target chromosome: Mpo {args.target_chr}\n")
        handle.write(f"- Major R108 gene-anchor sources on Mpo {args.target_chr}: {top_gene_sources}\n")
        handle.write(f"- Active Mpo Chr4 CENH3 core: {args.target_chr}:{active_start}-{active_end}; reciprocal nearest R108 centromere: {active_chr4['nearest_R108_CEN']}\n")
        handle.write(f"- R108 CEN5 top projection cluster: {args.target_chr}:{cen5['top_start']}-{cen5['top_end']}; {cen5['top_block_count']} blocks and {cen5['top_sum_query_len']} bp projected query length\n")
        handle.write(f"- Distance from CEN5 top cluster to active Mpo Chr4 CENH3 core: {interval_distance(cen5['top_start'], cen5['top_end'], active_start, active_end)} bp\n")
        handle.write(f"- Distance from CEN5 top cluster to nearest R108 CEN4 projection cluster: {nearest_cen4_cluster_distance} bp\n")
        handle.write("\nInterpretation: Mpo Chr4 carries an active CEN4-like CENH3 core and a nearby CEN5-associated sequence cluster. This supports a fusion/remnant model, but chromosome painting or breakpoint-level evidence is still needed before claiming centromere inactivation as a proven mechanism.\n")

    print(out_dir / "Mpo_Chr4_CEN4_CEN5_fusion_context.tsv")
    print(out_dir / "Mpo_Chr4_R108_gene_anchor_distribution.tsv")
    print(out_dir / "Mpo_Chr4_R108_gene_anchor_1Mb_windows.tsv")
    print(summary_path)


if __name__ == "__main__":
    main()
