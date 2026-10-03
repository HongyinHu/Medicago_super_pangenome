#!/usr/bin/env python3
"""Measure orthogonal chromatin features around the Mpo CEN5 candidate."""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from pathlib import Path

import pyBigWig


def read_chrom_sizes(fai: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with fai.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, size, *_ = line.rstrip("\n").split("\t")
            sizes[chrom] = int(size)
    return sizes


def clamp(chrom: str, start: int, end: int, sizes: dict[str, int]) -> tuple[int, int]:
    if chrom not in sizes:
        raise ValueError(f"Missing chromosome size for {chrom}")
    start = max(0, start)
    end = min(sizes[chrom], end)
    if end <= start:
        raise ValueError(f"Invalid interval after clamping: {chrom}:{start}-{end}")
    return start, end


def load_top_cen5_cluster(path: Path) -> dict[str, object]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["ancestral_centromere_id"] == "CEN5":
                rows.append(row)
    if not rows:
        raise ValueError(f"No CEN5 clusters in {path}")
    rows.sort(key=lambda x: (int(x["sum_query_len"]), int(x["block_count"])), reverse=True)
    top = rows[0]
    return {
        "label": "R108_CEN5_projection_cluster_top",
        "category": "cen5_candidate",
        "chrom": top["Mpo_chr"],
        "start": int(top["cluster_start"]),
        "end": int(top["cluster_end"]),
        "block_count": top["block_count"],
        "sum_query_len": top["sum_query_len"],
    }


def load_active_cens(path: Path) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            out.append(
                {
                    "label": f"Mpo_active_CENH3_core_{chrom}",
                    "category": "active_cen_core",
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                }
            )
    return out


def expand(interval: dict[str, object], flank: int, sizes: dict[str, int], suffix: str) -> dict[str, object]:
    chrom = str(interval["chrom"])
    start, end = clamp(chrom, int(interval["start"]) - flank, int(interval["end"]) + flank, sizes)
    out = dict(interval)
    out["label"] = f"{interval['label']}_{suffix}"
    out["category"] = f"{interval['category']}_{suffix}"
    out["start"] = start
    out["end"] = end
    return out


def build_intervals(project_root: Path, sizes: dict[str, int]) -> list[dict[str, object]]:
    cen5 = load_top_cen5_cluster(project_root / "output_key/05_fusion_chr_CEN_fate/R108_CEN_projection_block_clusters.tsv")
    active = load_active_cens(project_root / "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.cen_core.bed")
    active_chr4 = next(x for x in active if x["chrom"] == "Chr4")

    intervals: list[dict[str, object]] = [
        cen5,
        expand(cen5, 500_000, sizes, "pm500kb"),
        expand(cen5, 1_000_000, sizes, "pm1Mb"),
        dict(active_chr4, label="Mpo_active_CENH3_core_Chr4_primary_compare"),
        expand(active_chr4, 500_000, sizes, "pm500kb"),
        expand(active_chr4, 1_000_000, sizes, "pm1Mb"),
    ]
    intervals.extend(active)
    for interval in intervals:
        start, end = clamp(str(interval["chrom"]), int(interval["start"]), int(interval["end"]), sizes)
        interval["start"] = start
        interval["end"] = end
        interval["length_bp"] = end - start
    return intervals


def bw_weighted_stats(path: Path, chrom: str, start: int, end: int) -> dict[str, object]:
    total = 0
    weighted_sum = 0.0
    min_val: float | None = None
    max_val: float | None = None
    positive = 0
    negative = 0
    zero = 0
    covered = 0
    with pyBigWig.open(str(path)) as bw:
        intervals = bw.intervals(chrom, start, end) or []
    for s, e, value in intervals:
        ov_start = max(start, int(s))
        ov_end = min(end, int(e))
        if ov_end <= ov_start:
            continue
        width = ov_end - ov_start
        if value is None or math.isnan(float(value)):
            continue
        val = float(value)
        covered += width
        total += width
        weighted_sum += val * width
        min_val = val if min_val is None else min(min_val, val)
        max_val = val if max_val is None else max(max_val, val)
        if val > 0:
            positive += width
        elif val < 0:
            negative += width
        else:
            zero += width
    interval_len = end - start
    mean = weighted_sum / total if total else float("nan")
    return {
        "mean": mean,
        "min": min_val if min_val is not None else "NA",
        "max": max_val if max_val is not None else "NA",
        "covered_fraction": covered / interval_len if interval_len else 0.0,
        "positive_fraction": positive / covered if covered else "NA",
        "negative_fraction": negative / covered if covered else "NA",
        "zero_fraction": zero / covered if covered else "NA",
    }


def load_bedgraph(path: Path) -> dict[str, list[tuple[int, int, float]]]:
    data: dict[str, list[tuple[int, int, float]]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, value, *_ = line.rstrip("\n").split()
            data.setdefault(chrom, []).append((int(start), int(end), float(value)))
    for chrom in data:
        data[chrom].sort()
    return data


def bedgraph_weighted_stats(data: dict[str, list[tuple[int, int, float]]], chrom: str, start: int, end: int) -> dict[str, object]:
    intervals = data.get(chrom, [])
    starts = [x[0] for x in intervals]
    idx = max(0, bisect.bisect_left(starts, start) - 1)
    total = 0
    weighted_sum = 0.0
    min_val: float | None = None
    max_val: float | None = None
    positive = 0
    negative = 0
    covered = 0
    while idx < len(intervals):
        s, e, val = intervals[idx]
        if s >= end:
            break
        ov_start = max(start, s)
        ov_end = min(end, e)
        if ov_end > ov_start:
            width = ov_end - ov_start
            covered += width
            total += width
            weighted_sum += val * width
            min_val = val if min_val is None else min(min_val, val)
            max_val = val if max_val is None else max(max_val, val)
            if val > 0:
                positive += width
            elif val < 0:
                negative += width
        idx += 1
    interval_len = end - start
    mean = weighted_sum / total if total else float("nan")
    return {
        "mean": mean,
        "min": min_val if min_val is not None else "NA",
        "max": max_val if max_val is not None else "NA",
        "covered_fraction": covered / interval_len if interval_len else 0.0,
        "positive_fraction": positive / covered if covered else "NA",
        "negative_fraction": negative / covered if covered else "NA",
    }


def load_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    data: dict[str, list[tuple[int, int]]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            data.setdefault(chrom, []).append((int(start), int(end)))
    for chrom in data:
        data[chrom].sort()
    return data


def bed_overlap_and_distance(data: dict[str, list[tuple[int, int]]], chrom: str, start: int, end: int) -> dict[str, object]:
    intervals = data.get(chrom, [])
    overlap_count = 0
    overlap_bp = 0
    min_distance: int | None = None
    for s, e in intervals:
        if e < start:
            dist = start - e
        elif s > end:
            dist = s - end
        else:
            dist = 0
            overlap_count += 1
            overlap_bp += max(0, min(end, e) - max(start, s))
        min_distance = dist if min_distance is None else min(min_distance, dist)
    return {
        "overlap_count": overlap_count,
        "overlap_bp": overlap_bp,
        "distance_to_nearest_bp": min_distance if min_distance is not None else "NA",
    }


def fmt(value: object) -> object:
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.6f}"
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    args = parser.parse_args()

    project_root = Path(args.project_root)
    out_dir = project_root / "output_key/06_repeat_gene_mappability_controls"
    out_dir.mkdir(parents=True, exist_ok=True)

    sizes = read_chrom_sizes(project_root / "output_all/00_genome/genome_Mpo.fa.fai")
    intervals = build_intervals(project_root, sizes)

    hic_dir = project_root / "output_all/08_hic_AB_TAD/genome_Mpo"
    meth_dir = project_root / "output_all/07_methylation/genome_Mpo.hifi_methylation"

    ab_bw = hic_dir / "genome_Mpo.AB.100kb.cis.bw"
    cpg_bw = meth_dir / "genome_Mpo.CpG_methylation_density.50k.bw"
    cpg_bg = load_bedgraph(meth_dir / "genome_Mpo.CpG_methylation_density.50k.bedGraph")
    tad_40 = load_bedgraph(hic_dir / "genome_Mpo.TAD.40kb_score.bedgraph")
    tad_100 = load_bedgraph(hic_dir / "genome_Mpo.TAD.100kb_score.bedgraph")
    tad_domains_40 = load_bed(hic_dir / "genome_Mpo.TAD.40kb_domains.bed")
    tad_domains_100 = load_bed(hic_dir / "genome_Mpo.TAD.100kb_domains.bed")
    tad_boundaries_40 = load_bed(hic_dir / "genome_Mpo.TAD.40kb_boundaries.bed")
    tad_boundaries_100 = load_bed(hic_dir / "genome_Mpo.TAD.100kb_boundaries.bed")

    rows: list[dict[str, object]] = []
    for interval in intervals:
        chrom = str(interval["chrom"])
        start = int(interval["start"])
        end = int(interval["end"])
        row: dict[str, object] = {
            "label": interval["label"],
            "category": interval["category"],
            "chrom": chrom,
            "start": start,
            "end": end,
            "length_bp": interval["length_bp"],
            "block_count": interval.get("block_count", ""),
            "sum_query_len": interval.get("sum_query_len", ""),
        }
        ab = bw_weighted_stats(ab_bw, chrom, start, end)
        cpg_bw_stats = bw_weighted_stats(cpg_bw, chrom, start, end)
        cpg_bg_stats = bedgraph_weighted_stats(cpg_bg, chrom, start, end)
        tad40 = bedgraph_weighted_stats(tad_40, chrom, start, end)
        tad100 = bedgraph_weighted_stats(tad_100, chrom, start, end)

        for prefix, stats in [
            ("AB_100kb", ab),
            ("CpG_methylation_50kb_bw", cpg_bw_stats),
            ("CpG_methylation_50kb_bedgraph", cpg_bg_stats),
            ("TAD_score_40kb", tad40),
            ("TAD_score_100kb", tad100),
        ]:
            for key, value in stats.items():
                row[f"{prefix}_{key}"] = fmt(value)

        for prefix, bed_data in [
            ("TAD_domains_40kb", tad_domains_40),
            ("TAD_domains_100kb", tad_domains_100),
            ("TAD_boundaries_40kb", tad_boundaries_40),
            ("TAD_boundaries_100kb", tad_boundaries_100),
        ]:
            stats = bed_overlap_and_distance(bed_data, chrom, start, end)
            for key, value in stats.items():
                row[f"{prefix}_{key}"] = fmt(value)
        rows.append(row)

    out_tsv = out_dir / "Mpo_CEN5_candidate_orthogonal_chromatin_features.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    candidate = next(x for x in rows if x["label"] == "R108_CEN5_projection_cluster_top")
    candidate_pm1 = next(x for x in rows if x["label"] == "R108_CEN5_projection_cluster_top_pm1Mb")
    active_chr4 = next(x for x in rows if x["label"] == "Mpo_active_CENH3_core_Chr4_primary_compare")
    summary = out_dir / "Mpo_CEN5_candidate_orthogonal_chromatin_features.summary.md"
    with summary.open("w") as handle:
        handle.write("# Mpo CEN5 Candidate Orthogonal Chromatin Features\n\n")
        handle.write(f"- Candidate interval: {candidate['chrom']}:{candidate['start']}-{candidate['end']}\n")
        handle.write(f"- Active Chr4 CENH3 core: {active_chr4['chrom']}:{active_chr4['start']}-{active_chr4['end']}\n")
        handle.write(f"- Candidate AB mean: {candidate['AB_100kb_mean']}; active Chr4 CEN AB mean: {active_chr4['AB_100kb_mean']}\n")
        handle.write(f"- Candidate CpG methylation mean (50 kb bw): {candidate['CpG_methylation_50kb_bw_mean']}; active Chr4 CEN: {active_chr4['CpG_methylation_50kb_bw_mean']}\n")
        handle.write(f"- Candidate 40 kb TAD score mean: {candidate['TAD_score_40kb_mean']}; active Chr4 CEN: {active_chr4['TAD_score_40kb_mean']}\n")
        handle.write(f"- Candidate +/-1 Mb 40 kb TAD boundary count: {candidate_pm1['TAD_boundaries_40kb_overlap_count']}; active Chr4 CEN boundary count: {active_chr4['TAD_boundaries_40kb_overlap_count']}\n")
        handle.write(f"- Output table: `{out_tsv}`\n")

    print(out_tsv)
    print(summary)


if __name__ == "__main__":
    main()
