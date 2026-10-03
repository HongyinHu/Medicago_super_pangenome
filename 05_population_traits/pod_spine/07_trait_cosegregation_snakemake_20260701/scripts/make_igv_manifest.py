#!/usr/bin/env python3
import argparse
import csv
import os


def read_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--bam-manifest", required=True)
    parser.add_argument("--padding", type=int, default=2000)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    targets = {r["target_id"]: r for r in read_table(args.targets)}
    bams = read_table(args.bam_manifest)
    candidates = read_table(args.candidates)

    fields = ["candidate_id", "ref", "chrom", "view_start", "view_end", "sample", "bam", "track_label"]
    with open(args.out, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for cand in candidates:
            if cand.get("pass_80pct") != "1":
                continue
            ref = cand.get("ref")
            chrom = cand.get("chrom")
            start = cand.get("start")
            end = cand.get("end")
            if not (ref and chrom and start and end):
                target = targets.get(cand["candidate_id"])
                if not target:
                    continue
                ref = target["ref"]
                chrom = target["chrom"]
                start = target["start"]
                end = target["end"]
            view_start = max(1, int(start) - args.padding)
            view_end = int(end) + args.padding
            for bam in bams:
                if bam["ref"] != ref:
                    continue
                writer.writerow({
                    "candidate_id": cand["candidate_id"],
                    "ref": ref,
                    "chrom": chrom,
                    "view_start": view_start,
                    "view_end": view_end,
                    "sample": bam["sample"],
                    "bam": bam["bam"],
                    "track_label": f"{bam['sample']} {ref}",
                })


if __name__ == "__main__":
    main()
