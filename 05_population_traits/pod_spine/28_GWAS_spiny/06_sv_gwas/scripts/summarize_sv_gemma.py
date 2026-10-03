#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
from pathlib import Path


TARGET_CHROM = "Chr4"
TARGET_START = 88_980_831
TARGET_END = 88_981_052
LOCAL_START = TARGET_START - 1_000_000
LOCAL_END = TARGET_END + 1_000_000


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assoc", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def load_mapping(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return {row["marker_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def reciprocal_overlap(start1, end1, start2, end2):
    overlap = max(0, min(end1, end2) - max(start1, start2) + 1)
    if overlap == 0:
        return 0.0
    return min(overlap / (end1 - start1 + 1), overlap / (end2 - start2 + 1))


def lambda_gc(values):
    from scipy.stats import chi2

    chisq = sorted(float(chi2.isf(p, 1)) for p in values if 0 < p < 1)
    if not chisq:
        return float("nan")
    mid = len(chisq) // 2
    median = chisq[mid] if len(chisq) % 2 else (chisq[mid - 1] + chisq[mid]) / 2
    return median / 0.4549364


def main():
    args = parse_args()
    mapping = load_mapping(args.mapping)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []

    with open(args.assoc, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            marker = row.get("rs", "")
            meta = mapping.get(marker)
            if meta is None:
                raise SystemExit(f"Association marker missing from mapping: {marker}")
            try:
                p = float(row["p_wald"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (0 < p <= 1):
                continue
            merged = {
                "class": "SV",
                "chr": meta["chrom"],
                "pos": int(meta["pos"]),
                "end": int(meta["end"]),
                "id": meta["original_id"],
                "marker_id": marker,
                "svtype": meta["svtype"],
                "svlen": meta["svlen"],
                "p_wald": p,
                "minus_log10_p": -math.log10(p),
                "beta": row.get("beta", "NA"),
                "se": row.get("se", "NA"),
                "af": row.get("af", "NA"),
                "n_miss": row.get("n_miss", "NA"),
            }
            rows.append(merged)

    rows.sort(key=lambda item: item["p_wald"])
    fields = [
        "class", "chr", "pos", "end", "id", "marker_id", "svtype", "svlen",
        "p_wald", "minus_log10_p", "beta", "se", "af", "n_miss"
    ]

    def write_tsv(path, selected, compressed=False):
        opener = gzip.open if compressed else open
        kwargs = {"mode": "wt" if compressed else "w", "encoding": "utf-8", "newline": ""}
        with opener(path, **kwargs) as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(selected)

    write_tsv(out / "SV.primary.assoc.tsv.gz", rows, compressed=True)
    write_tsv(out / "top_hits.tsv", rows[:200])
    local = [
        row for row in rows
        if row["chr"] == TARGET_CHROM and LOCAL_START <= row["pos"] <= LOCAL_END
    ]
    write_tsv(out / "Chr23997_plusminus1Mb.tsv", local)

    exact = []
    for row in local:
        if row["svtype"] != "DEL":
            continue
        bp_match = abs(row["pos"] - TARGET_START) <= 100 and abs(row["end"] - TARGET_END) <= 100
        ro = reciprocal_overlap(row["pos"], row["end"], TARGET_START, TARGET_END)
        if bp_match or ro >= 0.5:
            exact.append(row)
    write_tsv(out / "Chr23997_exact_target_like.tsv", exact)

    with open(out / "association_summary.tsv", "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"tested_SV\t{len(rows)}\n")
        handle.write(f"bonferroni_0.05\t{(0.05 / len(rows) if rows else float('nan')):.12g}\n")
        handle.write(f"suggestive_1_over_n\t{(1 / len(rows) if rows else float('nan')):.12g}\n")
        handle.write(f"lambda_gc\t{lambda_gc([row['p_wald'] for row in rows]):.6g}\n")
        handle.write(f"Chr23997_window_markers\t{len(local)}\n")
        handle.write(f"Chr23997_exact_target_like_markers\t{len(exact)}\n")
        if local:
            best = min(local, key=lambda item: item["p_wald"])
            handle.write(f"Chr23997_window_best_id\t{best['id']}\n")
            handle.write(f"Chr23997_window_best_position\t{best['pos']}\n")
            handle.write(f"Chr23997_window_best_p_wald\t{best['p_wald']:.12g}\n")
        if exact:
            best = min(exact, key=lambda item: item["p_wald"])
            handle.write(f"Chr23997_exact_best_id\t{best['id']}\n")
            handle.write(f"Chr23997_exact_best_interval\t{best['chr']}:{best['pos']}-{best['end']}\n")
            handle.write(f"Chr23997_exact_best_p_wald\t{best['p_wald']:.12g}\n")
        else:
            handle.write("Chr23997_exact_target_genotyped\tno\n")

    print(f"Summarized {len(rows)} SV tests; local={len(local)}; exact_target_like={len(exact)}")


if __name__ == "__main__":
    main()
