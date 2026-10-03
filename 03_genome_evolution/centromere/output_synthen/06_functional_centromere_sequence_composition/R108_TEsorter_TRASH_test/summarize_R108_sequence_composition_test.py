#!/usr/bin/env python3
import pandas as pd

wide = pd.read_csv("genome_R108.functional_centromere_sequence_composition.broad.wide.tsv", sep="\t")
fine = pd.read_csv("genome_R108.functional_centromere_sequence_composition.TEsorter_fine.tsv", sep="\t")
avg = pd.read_csv("genome_R108.functional_centromere_sequence_composition.broad.average.tsv", sep="\t")

ratio_cols = [c for c in wide.columns if c.endswith("_ratio_pct")]
wide["ratio_sum_pct"] = wide[ratio_cols].sum(axis=1)
print("== broad wide ==")
print(wide.to_string(index=False))
print("\n== ratio sum check ==")
print(wide[["centromere", "chr", "ratio_sum_pct"]].to_string(index=False))
print("\n== broad average ==")
print(avg.to_string(index=False))
print("\n== fine top by total bp ==")
top = (
    fine.groupby(["broad_category", "fine_category"], as_index=False)
    .agg(length_bp=("length_bp", "sum"), mean_ratio_pct=("ratio_pct", "mean"))
    .sort_values("length_bp", ascending=False)
    .head(30)
)
print(top.to_string(index=False))
