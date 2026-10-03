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
    parser.add_argument("--snp-indel", required=True)
    parser.add_argument("--sv", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_rows(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def normalize_chrom(value):
    value = str(value)
    return value if value.startswith("Chr") else f"Chr{value}"


def main():
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []

    for row in read_rows(args.snp_indel):
        try:
            p = float(row["p_wald"])
            pos = int(row["pos"])
        except (KeyError, TypeError, ValueError):
            continue
        if not (0 < p <= 1):
            continue
        rows.append(
            {
                "class": row["class"],
                "chr": normalize_chrom(row["chr"]),
                "pos": pos,
                "end": pos,
                "id": row["id"],
                "svtype": "NA",
                "p_wald": p,
                "minus_log10_p": -math.log10(p),
                "beta": row.get("beta", "NA"),
                "se": row.get("se", "NA"),
                "af": row.get("af", "NA"),
                "n_miss": row.get("n_miss", "NA"),
            }
        )

    for row in read_rows(args.sv):
        try:
            p = float(row["p_wald"])
            pos = int(row["pos"])
            end = int(row["end"])
        except (KeyError, TypeError, ValueError):
            continue
        if not (0 < p <= 1):
            continue
        rows.append(
            {
                "class": "SV",
                "chr": normalize_chrom(row["chr"]),
                "pos": pos,
                "end": end,
                "id": row["id"],
                "svtype": row.get("svtype", "NA"),
                "p_wald": p,
                "minus_log10_p": -math.log10(p),
                "beta": row.get("beta", "NA"),
                "se": row.get("se", "NA"),
                "af": row.get("af", "NA"),
                "n_miss": row.get("n_miss", "NA"),
            }
        )

    chr_rank = {f"Chr{i}": i for i in range(1, 9)}
    rows.sort(key=lambda item: (chr_rank.get(item["chr"], 999), item["pos"], item["class"]))
    fields = [
        "class", "chr", "pos", "end", "id", "svtype", "p_wald",
        "minus_log10_p", "beta", "se", "af", "n_miss"
    ]

    combined = out / "SNP_INDEL_SV.primary.assoc.tsv.gz"
    with gzip.open(combined, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    local = [
        row for row in rows
        if row["chr"] == TARGET_CHROM and LOCAL_START <= row["pos"] <= LOCAL_END
    ]
    with open(out / "Chr23997_plusminus1Mb.all_classes.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(local, key=lambda item: item["p_wald"]))

    counts = {name: sum(row["class"] == name for row in rows) for name in ("SNP", "INDEL", "SV")}
    with open(out / "joint_association_summary.tsv", "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"tested_total\t{len(rows)}\n")
        for name in ("SNP", "INDEL", "SV"):
            handle.write(f"tested_{name}\t{counts[name]}\n")
        handle.write(f"bonferroni_0.05\t{0.05 / len(rows):.12g}\n")
        handle.write(f"suggestive_1_over_n\t{1 / len(rows):.12g}\n")
        handle.write(f"Chr23997_window_markers\t{len(local)}\n")
        for name in ("SNP", "INDEL", "SV"):
            subset = [row for row in local if row["class"] == name]
            handle.write(f"Chr23997_window_{name}_markers\t{len(subset)}\n")
            if subset:
                best = min(subset, key=lambda item: item["p_wald"])
                handle.write(f"Chr23997_window_{name}_best_id\t{best['id']}\n")
                handle.write(f"Chr23997_window_{name}_best_position\t{best['pos']}\n")
                handle.write(f"Chr23997_window_{name}_best_p_wald\t{best['p_wald']:.12g}\n")
        handle.write("Chr23997_exact_221bp_DEL_jointly_genotyped\tno\n")
        handle.write("Chr23997_exact_221bp_DEL_interval\tChr4:88980831-88981052\n")

    print(f"Integrated {len(rows)} tests; Chr23997 window={len(local)}")


if __name__ == "__main__":
    main()
