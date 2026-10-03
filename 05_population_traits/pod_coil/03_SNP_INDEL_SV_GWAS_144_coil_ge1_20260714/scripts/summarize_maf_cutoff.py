#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from summarize_gmmat import candidate_summary, draw_manhattan, draw_qq, lambda_gc, load_assoc


def normalise(frame):
    result = frame.copy()
    result["maf"] = np.minimum(pd.to_numeric(result["AF"], errors="coerce"), 1 - pd.to_numeric(result["AF"], errors="coerce"))
    result["effect"] = pd.to_numeric(result.get("beta_score_approx"), errors="coerce")
    result["se"] = pd.to_numeric(result.get("se_score_approx"), errors="coerce")
    return result


def write_top20(label, frame, outdir, prefix):
    columns = ["chrom", "pos", "id", "p_value", "maf", "effect", "se", "AF", "SCORE", "VAR"]
    top = frame.nsmallest(20, "p_value").reindex(columns=columns)
    top = top.rename(columns={"chrom": "CHR", "pos": "POS", "id": "SNP", "p_value": "P", "maf": "MAF", "effect": "Effect_score_approx", "se": "SE_score_approx"})
    top.insert(0, "Variant_class", label)
    top.to_csv(outdir / "top20" / f"{prefix}_{label}_top20.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assoc", action="append", required=True, help="CLASS=path")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--maf", type=float, required=True)
    args = parser.parse_args()

    out = Path(args.outdir)
    (out / "top20").mkdir(parents=True, exist_ok=True)
    frames = {}
    rows = []
    for item in args.assoc:
        label, path = item.split("=", 1)
        frame = normalise(load_assoc(label, path))
        if frame.empty:
            raise RuntimeError(f"No valid tests in {path}")
        if frame.maf.dropna().min() + 1e-9 < args.maf:
            raise RuntimeError(f"{label} contains a marker below MAF {args.maf}")
        frames[label] = frame
        write_top20(label, frame, out, args.prefix)
        best = frame.loc[frame.p_value.idxmin()]
        rows.append({
            "MAF_cutoff": args.maf,
            "Variant_class": label,
            "Variant_number": len(frame),
            "lambdaGC": lambda_gc(frame.p_value),
            "Best_SNP": best.id,
            "Best_CHR": int(best.chrom),
            "Best_POS": int(best.pos),
            "Best_P": float(best.p_value),
            "Best_MAF": float(best.maf),
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(out / f"{args.prefix}.variant_summary.tsv", sep="\t", index=False)
    combined = pd.concat(frames.values(), ignore_index=True).sort_values("p_value")
    combined.to_csv(out / f"{args.prefix}.assoc.tsv.gz", sep="\t", index=False, compression="gzip")
    combined.head(500).to_csv(out / f"{args.prefix}.top500.tsv", sep="\t", index=False)
    candidate_summary(frames, out)
    (out / "candidate_gene_peak_summary.tsv").replace(out / f"{args.prefix}.candidate_gene_peak_summary.tsv")
    draw_qq(frames, out, args.prefix)
    draw_manhattan(frames, out, args.prefix)
    print(f"{args.prefix}: summarized {len(combined):,} tests")


if __name__ == "__main__":
    main()
