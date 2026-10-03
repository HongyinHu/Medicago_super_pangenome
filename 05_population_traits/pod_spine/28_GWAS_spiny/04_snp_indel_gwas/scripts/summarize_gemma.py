#!/usr/bin/env python3
import csv
import gzip
import math
import sys
from pathlib import Path


if len(sys.argv) != 4:
    raise SystemExit("Usage: summarize_gemma.py snp.assoc indel.assoc outdir")

snp_path, indel_path, outdir = sys.argv[1:]
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)


def load(path, variant_class):
    rows = []
    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                p = float(row["p_wald"])
                chrom = str(row["chr"])
                pos = int(float(row["ps"]))
            except (KeyError, TypeError, ValueError):
                continue
            if not (0 < p <= 1):
                continue
            rows.append(
                {
                    "class": variant_class,
                    "chr": chrom,
                    "pos": pos,
                    "id": row.get("rs", f"{chrom}:{pos}"),
                    "p_wald": p,
                    "minus_log10_p": -math.log10(p),
                    "beta": row.get("beta", "NA"),
                    "se": row.get("se", "NA"),
                    "af": row.get("af", "NA"),
                    "n_miss": row.get("n_miss", "NA"),
                }
            )
    return rows


rows = load(snp_path, "SNP") + load(indel_path, "INDEL")
rows.sort(key=lambda x: x["p_wald"])
fields = ["class", "chr", "pos", "id", "p_wald", "minus_log10_p", "beta", "se", "af", "n_miss"]

with gzip.open(out / "SNP_INDEL.primary.assoc.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

with open(out / "top_hits.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows[:200])

local = [
    row
    for row in rows
    if row["chr"].replace("Chr", "") == "4" and 87_980_831 <= row["pos"] <= 89_981_052
]
with open(out / "Chr23997_plusminus1Mb.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(local)


def lambda_gc(values):
    # Wald p-values are upper-tail probabilities from a 1-df chi-square test.
    from scipy.stats import chi2

    chisq = [float(chi2.isf(p, 1)) for p in values if 0 < p < 1]
    if not chisq:
        return float("nan")
    chisq.sort()
    mid = len(chisq) // 2
    median = chisq[mid] if len(chisq) % 2 else (chisq[mid - 1] + chisq[mid]) / 2
    return median / 0.4549364


with open(out / "association_summary.tsv", "w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"tested_total\t{len(rows)}\n")
    handle.write(f"tested_SNP\t{sum(r['class'] == 'SNP' for r in rows)}\n")
    handle.write(f"tested_INDEL\t{sum(r['class'] == 'INDEL' for r in rows)}\n")
    handle.write(f"bonferroni_0.05\t{0.05 / len(rows):.12g}\n")
    handle.write(f"suggestive_1_over_n\t{1 / len(rows):.12g}\n")
    handle.write(f"lambda_gc_all\t{lambda_gc([r['p_wald'] for r in rows]):.6g}\n")
    handle.write(f"Chr23997_window_markers\t{len(local)}\n")
    if local:
        best = min(local, key=lambda x: x["p_wald"])
        handle.write(f"Chr23997_best_class\t{best['class']}\n")
        handle.write(f"Chr23997_best_id\t{best['id']}\n")
        handle.write(f"Chr23997_best_position\t{best['pos']}\n")
        handle.write(f"Chr23997_best_p_wald\t{best['p_wald']:.12g}\n")

print(f"Summarized {len(rows)} SNP/INDEL association tests; local markers={len(local)}")
