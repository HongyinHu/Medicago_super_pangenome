#!/usr/bin/env python3
import argparse
import csv
import os
from collections import Counter, defaultdict

META_COLS = {
    "SV_ID", "source_methods", "reference", "CHROM", "POS", "END", "SVTYPE", "SVLEN",
    "read_SV_ID", "svgap_ids", "final_confidence"
}

CATS = ["Intergenic", "2kb Upstream", "2kb downstream", "Intron", "Exon"]


def presence_value(v):
    if v is None:
        return False
    s = str(v).strip()
    if s in ("", ".", "NA", "NaN", "nan"):
        return False
    try:
        return float(s) > 0
    except ValueError:
        return s not in ("0", "./.", "0/0")


def classify_gene_structure(row, flank_bp):
    detail = row.get("gene_body_detail", "")
    if detail == "exon_overlap":
        return "Exon"
    if detail == "intron_only":
        return "Intron"

    signed = row.get("signed_distance_bp", ".")
    try:
        signed = int(float(signed))
    except ValueError:
        return "Intergenic"
    category4 = row.get("category4", "")
    if category4 == "Upstream" and abs(signed) <= flank_bp:
        return "2kb Upstream"
    if category4 == "Downstream" and abs(signed) <= flank_bp:
        return "2kb downstream"
    return "Intergenic"


def main():
    ap = argparse.ArgumentParser(description="Per-species PAV gene-structure distribution excluding genome_Msa.")
    ap.add_argument("--pav", required=True)
    ap.add_argument("--gene-context-events", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--exclude-sample", default="genome_Msa")
    ap.add_argument("--flank-bp", type=int, default=2000)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)

    cat_by_id = {}
    with open(args.gene_context_events, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sv_id = row.get("SV_ID")
            if sv_id:
                cat_by_id[sv_id] = classify_gene_structure(row, args.flank_bp)

    counts = defaultdict(Counter)
    totals = Counter()
    missing_context = 0

    with open(args.pav, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        sample_cols = [c for c in reader.fieldnames if c not in META_COLS and c != args.exclude_sample]
        for row in reader:
            sv_id = row.get("SV_ID")
            cat = cat_by_id.get(sv_id)
            if cat is None:
                # Events absent from remaining 17 species were skipped in the previous gene-context table.
                missing_context += 1
                continue
            for sample in sample_cols:
                if presence_value(row.get(sample)):
                    counts[sample][cat] += 1
                    totals[sample] += 1

    out_long = args.out_prefix + ".per_species.tsv"
    with open(out_long, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "category", "count", "total_PAV", "proportion_percent"])
        for sample in sample_cols:
            total = totals[sample]
            for cat in CATS:
                cnt = counts[sample][cat]
                pct = cnt / total * 100 if total else 0
                w.writerow([sample, cat, cnt, total, f"{pct:.6f}"])

    out_summary = args.out_prefix + ".category_summary.tsv"
    vals = {cat: [] for cat in CATS}
    for sample in sample_cols:
        total = totals[sample]
        for cat in CATS:
            vals[cat].append(counts[sample][cat] / total * 100 if total else 0)

    def quantile(x, q):
        if not x:
            return "."
        x = sorted(x)
        pos = (len(x) - 1) * q
        lo = int(pos)
        hi = min(lo + 1, len(x) - 1)
        frac = pos - lo
        return x[lo] * (1 - frac) + x[hi] * frac

    with open(out_summary, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["category", "n_species", "mean_percent", "median_percent", "q1_percent", "q3_percent", "min_percent", "max_percent"])
        for cat in CATS:
            x = vals[cat]
            w.writerow([
                cat, len(x), f"{sum(x)/len(x):.6f}", f"{quantile(x, 0.5):.6f}",
                f"{quantile(x, 0.25):.6f}", f"{quantile(x, 0.75):.6f}",
                f"{min(x):.6f}", f"{max(x):.6f}",
            ])

    out_run = args.out_prefix + ".run_summary.tsv"
    with open(out_run, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerow(["input_pav", args.pav])
        w.writerow(["gene_context_events", args.gene_context_events])
        w.writerow(["excluded_sample", args.exclude_sample])
        w.writerow(["flank_bp", args.flank_bp])
        w.writerow(["sample_columns", ",".join(sample_cols)])
        w.writerow(["missing_context_rows", missing_context])

    print(f"DONE samples={len(sample_cols)} missing_context_rows={missing_context}")
    print(out_long)
    print(out_summary)


if __name__ == "__main__":
    main()
