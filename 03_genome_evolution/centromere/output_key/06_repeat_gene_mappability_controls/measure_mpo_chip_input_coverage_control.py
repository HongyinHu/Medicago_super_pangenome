#!/usr/bin/env python3
"""Measure CENH3/Input coverage controls for the Mpo CEN5 candidate interval."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path


EXCLUDE_FLAGS = "3332"  # unmapped, secondary, duplicate, supplementary
PSEUDO_RPKM = 0.01


def run(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return proc.stdout


def read_chrom_sizes(fai: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with fai.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            sizes[fields[0]] = int(fields[1])
    return sizes


def clamp_interval(chrom: str, start: int, end: int, sizes: dict[str, int]) -> tuple[int, int]:
    chrom_len = sizes.get(chrom)
    if chrom_len is None:
        raise ValueError(f"Missing chromosome size for {chrom}")
    start = max(0, start)
    end = min(chrom_len, end)
    if end <= start:
        raise ValueError(f"Invalid interval after clamping: {chrom}:{start}-{end}")
    return start, end


def load_active_cens(path: Path) -> list[dict[str, object]]:
    intervals: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            intervals.append(
                {
                    "label": f"Mpo_active_CENH3_core_{chrom}",
                    "category": "active_cen_core",
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                }
            )
    return intervals


def load_top_cen5_cluster(path: Path) -> dict[str, object]:
    candidates: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["ancestral_centromere_id"] != "CEN5":
                continue
            row["_sum_query_len"] = int(row["sum_query_len"])
            row["_block_count"] = int(row["block_count"])
            candidates.append(row)
    if not candidates:
        raise ValueError(f"No CEN5 rows found in {path}")
    top = max(candidates, key=lambda x: (x["_sum_query_len"], x["_block_count"]))
    return {
        "label": "R108_CEN5_projection_cluster_top",
        "category": "cen5_candidate",
        "chrom": top["Mpo_chr"],
        "start": int(top["cluster_start"]),
        "end": int(top["cluster_end"]),
        "block_count": str(top["block_count"]),
        "sum_query_len": str(top["sum_query_len"]),
        "R108_query_chr": top["R108_query_chr"],
        "R108_query_min": top["R108_query_min"],
        "R108_query_max": top["R108_query_max"],
    }


def expand_interval(interval: dict[str, object], flank: int, sizes: dict[str, int], label_suffix: str) -> dict[str, object]:
    chrom = str(interval["chrom"])
    start, end = clamp_interval(chrom, int(interval["start"]) - flank, int(interval["end"]) + flank, sizes)
    out = dict(interval)
    out["label"] = f"{interval['label']}_{label_suffix}"
    out["category"] = f"{interval['category']}_{label_suffix}"
    out["start"] = start
    out["end"] = end
    return out


def total_primary_mapped(bam: Path) -> int:
    out = run(["samtools", "view", "-c", "-F", EXCLUDE_FLAGS, str(bam)])
    return int(out.strip())


def coverage_for_region(bam: Path, chrom: str, start: int, end: int) -> dict[str, float]:
    # samtools uses 1-based closed coordinates in region strings.
    region = f"{chrom}:{start + 1}-{end}"
    out = run(["samtools", "coverage", "--ff", EXCLUDE_FLAGS, "-r", region, str(bam)])
    for line in out.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        return {
            "numreads": float(fields[3]),
            "covbases": float(fields[4]),
            "coverage_pct": float(fields[5]),
            "meandepth": float(fields[6]),
            "meanbaseq": float(fields[7]),
            "meanmapq": float(fields[8]),
        }
    return {
        "numreads": 0.0,
        "covbases": 0.0,
        "coverage_pct": 0.0,
        "meandepth": 0.0,
        "meanbaseq": 0.0,
        "meanmapq": 0.0,
    }


def build_intervals(project_root: Path, sizes: dict[str, int]) -> list[dict[str, object]]:
    active_path = project_root / "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.cen_core.bed"
    cluster_path = project_root / "output_key/05_fusion_chr_CEN_fate/R108_CEN_projection_block_clusters.tsv"

    cen5 = load_top_cen5_cluster(cluster_path)
    active = load_active_cens(active_path)

    intervals: list[dict[str, object]] = [cen5]
    intervals.append(expand_interval(cen5, 500_000, sizes, "pm500kb"))
    intervals.append(expand_interval(cen5, 1_000_000, sizes, "pm1Mb"))

    intervals.extend(active)

    active_chr4 = next((x for x in active if x["chrom"] == "Chr4"), None)
    if active_chr4:
        intervals.insert(3, dict(active_chr4, label="Mpo_active_CENH3_core_Chr4_primary_compare"))

    for interval in intervals:
        start, end = clamp_interval(str(interval["chrom"]), int(interval["start"]), int(interval["end"]), sizes)
        interval["start"] = start
        interval["end"] = end
        interval["length_bp"] = end - start
    return intervals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    parser.add_argument("--out-dir", default=None)
    args = parser.parse_args()

    project_root = Path(args.project_root)
    out_dir = Path(args.out_dir) if args.out_dir else project_root / "output_key/06_repeat_gene_mappability_controls"
    out_dir.mkdir(parents=True, exist_ok=True)

    fai = project_root / "output_all/00_genome/genome_Mpo.fa.fai"
    sizes = read_chrom_sizes(fai)
    intervals = build_intervals(project_root, sizes)

    bams = {
        "CENH3_q20": project_root / "output_all/02_chipseq/genome_Mpo/03_chipseq_mapping/genome_Mpo_CENH3.map.q20.rmdup.bam",
        "Input_q20": project_root / "output_all/02_chipseq/genome_Mpo/03_chipseq_mapping/genome_Mpo_Input.map.q20.rmdup.bam",
        "CENH3_q0": project_root / "output_all/02_chipseq/genome_Mpo/03_chipseq_mapping/genome_Mpo_CENH3.map.q0.rmdup.bam",
        "Input_q0": project_root / "output_all/02_chipseq/genome_Mpo/03_chipseq_mapping/genome_Mpo_Input.map.q0.rmdup.bam",
    }

    missing = [str(path) for path in bams.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing BAM files:\n" + "\n".join(missing))

    totals = {name: total_primary_mapped(path) for name, path in bams.items()}

    bed_path = out_dir / "Mpo_CEN5_candidate_CENH3_Input_coverage_control.intervals.bed"
    with bed_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for interval in intervals:
            writer.writerow([interval["chrom"], interval["start"], interval["end"], interval["label"]])

    rows: list[dict[str, object]] = []
    for interval in intervals:
        row: dict[str, object] = {
            "label": interval["label"],
            "category": interval["category"],
            "chrom": interval["chrom"],
            "start": interval["start"],
            "end": interval["end"],
            "length_bp": interval["length_bp"],
            "block_count": interval.get("block_count", ""),
            "sum_query_len": interval.get("sum_query_len", ""),
            "R108_query_chr": interval.get("R108_query_chr", ""),
            "R108_query_min": interval.get("R108_query_min", ""),
            "R108_query_max": interval.get("R108_query_max", ""),
        }
        length_kb = float(interval["length_bp"]) / 1000.0
        for name, bam in bams.items():
            stats = coverage_for_region(bam, str(interval["chrom"]), int(interval["start"]), int(interval["end"]))
            mapped_million = totals[name] / 1_000_000.0
            rpkm = stats["numreads"] / length_kb / mapped_million if length_kb > 0 and mapped_million > 0 else 0.0
            depth_per_million = stats["meandepth"] / mapped_million if mapped_million > 0 else 0.0
            row[f"{name}_total_primary_mapped"] = totals[name]
            row[f"{name}_numreads"] = f"{stats['numreads']:.0f}"
            row[f"{name}_rpkm"] = f"{rpkm:.6f}"
            row[f"{name}_meandepth"] = f"{stats['meandepth']:.6f}"
            row[f"{name}_meandepth_per_million"] = f"{depth_per_million:.6f}"
            row[f"{name}_coverage_pct"] = f"{stats['coverage_pct']:.6f}"
            row[f"{name}_meanmapq"] = f"{stats['meanmapq']:.6f}"

        cenh3_q20 = float(row["CENH3_q20_rpkm"])
        input_q20 = float(row["Input_q20_rpkm"])
        cenh3_q0 = float(row["CENH3_q0_rpkm"])
        input_q0 = float(row["Input_q0_rpkm"])
        row["log2_CENH3_over_Input_q20_rpkm"] = f"{math.log2((cenh3_q20 + PSEUDO_RPKM) / (input_q20 + PSEUDO_RPKM)):.6f}"
        row["log2_CENH3_over_Input_q0_rpkm"] = f"{math.log2((cenh3_q0 + PSEUDO_RPKM) / (input_q0 + PSEUDO_RPKM)):.6f}"
        row["CENH3_q20_over_q0_rpkm_fraction"] = f"{(cenh3_q20 / cenh3_q0):.6f}" if cenh3_q0 > 0 else "NA"
        row["Input_q20_over_q0_rpkm_fraction"] = f"{(input_q20 / input_q0):.6f}" if input_q0 > 0 else "NA"
        rows.append(row)

    out_path = out_dir / "Mpo_CEN5_candidate_CENH3_Input_coverage_control.tsv"
    fieldnames = list(rows[0].keys())
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary_path = out_dir / "Mpo_CEN5_candidate_CENH3_Input_coverage_control.summary.md"
    candidate = next(row for row in rows if row["label"] == "R108_CEN5_projection_cluster_top")
    chr4_active = next(row for row in rows if row["label"] == "Mpo_active_CENH3_core_Chr4_primary_compare")
    with summary_path.open("w") as handle:
        handle.write("# Mpo CEN5 Candidate CENH3/Input Coverage Control\n\n")
        handle.write(f"- Candidate interval: {candidate['chrom']}:{candidate['start']}-{candidate['end']}\n")
        handle.write(f"- Candidate q20 log2(CENH3/Input): {candidate['log2_CENH3_over_Input_q20_rpkm']}\n")
        handle.write(f"- Active Chr4 CEN q20 log2(CENH3/Input): {chr4_active['log2_CENH3_over_Input_q20_rpkm']}\n")
        handle.write(f"- Candidate CENH3 q20/q0 RPKM fraction: {candidate['CENH3_q20_over_q0_rpkm_fraction']}\n")
        handle.write(f"- Candidate Input q20/q0 RPKM fraction: {candidate['Input_q20_over_q0_rpkm_fraction']}\n")
        handle.write(f"- Output table: `{out_path}`\n")
        handle.write(f"- Interval BED: `{bed_path}`\n")

    print(out_path)
    print(summary_path)


if __name__ == "__main__":
    main()
