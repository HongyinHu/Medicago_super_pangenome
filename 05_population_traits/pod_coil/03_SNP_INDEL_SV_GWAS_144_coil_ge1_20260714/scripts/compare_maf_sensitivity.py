#!/usr/bin/env python3
import argparse
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cutoff", action="append", required=True, help="numeric_cutoff=results_dir")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    sources = {}
    summary_rows = []
    for item in args.cutoff:
        cutoff, location = item.split("=", 1)
        location = Path(location)
        prefix = f"Model_A_MAF{int(round(float(cutoff) * 100)):03d}"
        table = pd.read_csv(location / f"{prefix}.variant_summary.tsv", sep="\t")
        table["MAF_cutoff"] = float(cutoff)
        summary_rows.append(table)
        sources[float(cutoff)] = (location, prefix)
    long = pd.concat(summary_rows, ignore_index=True)
    wide = long.pivot(index="MAF_cutoff", columns="Variant_class", values=["Variant_number", "lambdaGC"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index().sort_values("MAF_cutoff")
    ordered = ["MAF_cutoff", "Variant_number_SNP", "Variant_number_INDEL", "Variant_number_SV", "lambdaGC_SNP", "lambdaGC_INDEL", "lambdaGC_SV"]
    wide.reindex(columns=ordered).to_csv(out / "Model_A_MAF_sensitivity.tsv", sep="\t", index=False)
    long.to_csv(out / "Model_A_MAF_sensitivity_long.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.7), sharex=True)
    for ax, label in zip(axes, ("SNP", "INDEL", "SV")):
        part = long[long.Variant_class == label].sort_values("MAF_cutoff")
        ax.plot(part.MAF_cutoff, part.lambdaGC, marker="o", color="#2C7FB8", lw=1.8)
        ax.axhline(1, color="#D73027", ls="--", lw=0.9)
        ax.set_title(label)
        ax.set_xlabel("MAF cutoff")
        ax.set_ylabel(r"$\lambda_{GC}$")
        ax.set_xticks(sorted(sources))
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "Model_A_MAF_sensitivity_lambda.png", dpi=400)
    fig.savefig(out / "Model_A_MAF_sensitivity_lambda.pdf", dpi=300)
    plt.close(fig)

    comparison_rows = []
    for label in ("SNP", "INDEL", "SV"):
        loaded = {}
        for cutoff, (location, prefix) in sources.items():
            frame = pd.read_csv(location / f"{prefix}.assoc.tsv.gz", sep="\t", compression="gzip")
            loaded[cutoff] = frame[frame["class"] == label][["id", "p_value"]].copy()
        for first, second in combinations(sorted(loaded), 2):
            a = loaded[first].rename(columns={"p_value": "p_a"})
            b = loaded[second].rename(columns={"p_value": "p_b"})
            merged = a.merge(b, on="id", how="inner")
            if len(merged) >= 3:
                correlation = pearsonr(-np.log10(merged.p_a), -np.log10(merged.p_b)).statistic
            else:
                correlation = np.nan
            top_a = set(a.nsmallest(20, "p_a").id)
            top_b = set(b.nsmallest(20, "p_b").id)
            overlap = len(top_a & top_b)
            comparison_rows.append({"Variant_class": label, "MAF_a": first, "MAF_b": second,
                                    "shared_tested_variants": len(merged), "minus_log10P_pearson_r": correlation,
                                    "top20_shared": overlap, "top20_jaccard": overlap / len(top_a | top_b)})
    pd.DataFrame(comparison_rows).to_csv(out / "Model_A_MAF_top20_overlap_and_pvalue_correlation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
