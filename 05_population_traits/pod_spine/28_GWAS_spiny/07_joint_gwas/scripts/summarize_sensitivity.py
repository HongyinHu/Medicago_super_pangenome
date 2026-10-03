#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

from scipy.stats import chi2


LOCAL_START = 87_980_831
LOCAL_END = 89_981_052


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-root", required=True)
    parser.add_argument("--sensitivity-dir", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def lambda_gc(values):
    chisq = sorted(float(chi2.isf(p, 1)) for p in values if 0 < p < 1)
    if not chisq:
        return float("nan")
    mid = len(chisq) // 2
    median = chisq[mid] if len(chisq) % 2 else (chisq[mid - 1] + chisq[mid]) / 2
    return median / 0.4549364


def load_assoc(path, variant_class):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                p = float(row["p_wald"])
                chrom = str(row["chr"]).replace("Chr", "")
                pos = int(float(row["ps"]))
            except (KeyError, TypeError, ValueError):
                continue
            if 0 < p <= 1:
                rows.append((p, chrom, pos, row.get("rs", f"{chrom}:{pos}"), variant_class))
    return rows


def main():
    args = parse_args()
    root = Path(args.stage_root)
    sens = Path(args.sensitivity_dir)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    records = []

    for model in ("Konly", "K_PC3"):
        for variant in ("SNP", "INDEL", "SV"):
            rows = load_assoc(sens / f"{variant}.{model}.assoc.txt", variant)
            local = [row for row in rows if row[1] == "4" and LOCAL_START <= row[2] <= LOCAL_END]
            best = min(local) if local else (float("nan"), "NA", -1, "NA", variant)
            records.append(
                {
                    "model": model,
                    "variant_class": variant,
                    "tested": len(rows),
                    "lambda_gc": f"{lambda_gc([row[0] for row in rows]):.6g}",
                    "Chr23997_window_markers": len(local),
                    "Chr23997_best_id": best[3],
                    "Chr23997_best_position": best[2],
                    "Chr23997_best_p_wald": best[0],
                }
            )

    primary = {
        "SNP": root / "04_snp_indel_gwas/results/gemma/output/SNP.pc5.assoc.txt",
        "INDEL": root / "04_snp_indel_gwas/results/gemma/output/INDEL.pc5.assoc.txt",
        "SV": root / "06_sv_gwas/results/gemma/output/SV.pc5.assoc.txt",
    }
    for variant, path in primary.items():
        rows = load_assoc(path, variant)
        local = [row for row in rows if row[1] == "4" and LOCAL_START <= row[2] <= LOCAL_END]
        best = min(local) if local else (float("nan"), "NA", -1, "NA", variant)
        records.append(
            {
                "model": "K_PC5_primary",
                "variant_class": variant,
                "tested": len(rows),
                "lambda_gc": f"{lambda_gc([row[0] for row in rows]):.6g}",
                "Chr23997_window_markers": len(local),
                "Chr23997_best_id": best[3],
                "Chr23997_best_position": best[2],
                "Chr23997_best_p_wald": best[0],
            }
        )

    fields = list(records[0])
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(f"Wrote {len(records)} model-by-variant sensitivity rows")


if __name__ == "__main__":
    main()
