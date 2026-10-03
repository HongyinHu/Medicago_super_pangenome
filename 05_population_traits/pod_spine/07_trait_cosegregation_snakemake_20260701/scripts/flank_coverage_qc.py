#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess


def read_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mean_depth(samtools, bam, chrom, start, end):
    start = max(1, int(start))
    end = max(start, int(end))
    region = f"{chrom}:{start}-{end}"
    try:
        proc = subprocess.run(
            [samtools, "depth", "-a", "-r", region, bam],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        return "NA", "samtools_not_found"
    if proc.returncode != 0:
        return "NA", proc.stderr.strip().replace("\t", " ")[:200]
    total = 0
    n = 0
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            total += int(parts[2])
            n += 1
    return (total / n if n else 0.0), "OK"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--states", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--bam-manifest", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--flank-size", type=int, default=1000)
    parser.add_argument("--min-flank-depth", type=float, default=5.0)
    parser.add_argument("--max-event-to-flank-depth-ratio-for-del", type=float, default=0.55)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    targets = {r["target_id"]: r for r in read_table(args.targets)}
    bams = {(r["sample"], r["ref"]): r["bam"] for r in read_table(args.bam_manifest)}
    states = read_table(args.states)

    fields = [
        "target_id", "sample", "ref", "chrom", "start", "end", "bam",
        "left_flank_mean_depth", "event_mean_depth", "right_flank_mean_depth",
        "flank_mean_depth", "event_to_flank_depth_ratio", "coverage_qc", "qc_message",
    ]
    with open(args.out, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for state in states:
            target = targets.get(state["target_id"])
            if not target:
                continue
            ref = target["ref"]
            sample = state["sample"]
            chrom = target["chrom"]
            start = int(target["start"])
            end = int(target["end"])
            bam = bams.get((sample, ref), "")
            if not bam or not os.path.exists(bam):
                writer.writerow({
                    "target_id": state["target_id"], "sample": sample, "ref": ref,
                    "chrom": chrom, "start": start, "end": end, "bam": bam or "NA",
                    "left_flank_mean_depth": "NA", "event_mean_depth": "NA",
                    "right_flank_mean_depth": "NA", "flank_mean_depth": "NA",
                    "event_to_flank_depth_ratio": "NA", "coverage_qc": "not_tested",
                    "qc_message": "bam_missing",
                })
                continue
            left, msg_l = mean_depth(args.samtools, bam, chrom, start - args.flank_size, start - 1)
            event, msg_e = mean_depth(args.samtools, bam, chrom, start, end)
            right, msg_r = mean_depth(args.samtools, bam, chrom, end + 1, end + args.flank_size)
            if "NA" in {left, event, right}:
                qc = "not_tested"
                ratio = "NA"
                flank = "NA"
                msg = ";".join([msg_l, msg_e, msg_r])
            else:
                flank = (left + right) / 2
                ratio = event / flank if flank else 999.0
                if flank < args.min_flank_depth:
                    qc = "low_flank_coverage"
                elif state["state"] == "large_DEL" and ratio <= args.max_event_to_flank_depth_ratio_for_del:
                    qc = "pass"
                elif state["state"] != "large_DEL" and ratio > args.max_event_to_flank_depth_ratio_for_del:
                    qc = "pass"
                else:
                    qc = "coverage_conflict"
                msg = "OK"
            writer.writerow({
                "target_id": state["target_id"], "sample": sample, "ref": ref,
                "chrom": chrom, "start": start, "end": end, "bam": bam,
                "left_flank_mean_depth": f"{left:.3f}" if isinstance(left, float) else left,
                "event_mean_depth": f"{event:.3f}" if isinstance(event, float) else event,
                "right_flank_mean_depth": f"{right:.3f}" if isinstance(right, float) else right,
                "flank_mean_depth": f"{flank:.3f}" if isinstance(flank, float) else flank,
                "event_to_flank_depth_ratio": f"{ratio:.3f}" if isinstance(ratio, float) else ratio,
                "coverage_qc": qc,
                "qc_message": msg,
            })


if __name__ == "__main__":
    main()
