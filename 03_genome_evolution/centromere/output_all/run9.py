import pandas as pd

SP="genome_474"

df = pd.read_csv(f"05_centromere_define/{SP}.TRASH_CENH3_array_summary.tsv", sep="\t")

summary = (
    df.groupby("monomer_len", as_index=False)
    .agg(
        total_array_len=("array_len", "sum"),
        total_CENH3_overlap_bp=("CENH3_overlap_bp", "sum"),
        array_count=("TRASH_id", "count"),
        chromosome_count=("TRASH_chr", "nunique"),
        max_array_len=("array_len", "max"),
        max_CENH3_overlap_bp=("CENH3_overlap_bp", "max")
    )
)

summary["estimated_total_copy"] = summary["total_array_len"] / summary["monomer_len"]

summary = summary.sort_values(
    ["total_CENH3_overlap_bp", "total_array_len", "chromosome_count"],
    ascending=False
)

summary.to_csv(f"05_centromere_define/{SP}.TRASH_CENH3_by_monomer_len.tsv", sep="\t", index=False)

print(f"written: {SP}.TRASH_CENH3_by_monomer_len.tsv")
