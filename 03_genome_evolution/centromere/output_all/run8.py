import pandas as pd

SP="genome_474"

infile = f"05_centromere_define/{SP}.TRASH_overlap_CENH3peak.with_header.tsv"
outfile = f"05_centromere_define/{SP}.TRASH_CENH3_array_summary.tsv"

cols = [
    "TRASH_chr", "TRASH_start", "TRASH_end", "TRASH_id", "monomer_len", "array_len",
    "peak_chr", "peak_start", "peak_end", "peak_id", "peak_score", "peak_strand",
    "peak_signal", "peak_pvalue", "peak_qvalue"
]

df = pd.read_csv(infile, sep="\t", header=None, names=cols)

# 转换为数值
for c in ["TRASH_start", "TRASH_end", "monomer_len", "array_len", "peak_start", "peak_end",
          "peak_score", "peak_signal", "peak_pvalue", "peak_qvalue"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# 计算每一行真实 overlap bp
df["overlap_bp"] = (
    df[["TRASH_end", "peak_end"]].min(axis=1)
    - df[["TRASH_start", "peak_start"]].max(axis=1)
)

df = df[df["overlap_bp"] > 0]

summary = (
    df.groupby(
        ["TRASH_chr", "TRASH_start", "TRASH_end", "TRASH_id", "monomer_len", "array_len"],
        as_index=False
    )
    .agg(
        CENH3_overlap_bp=("overlap_bp", "sum"),
        CENH3_peak_count=("peak_id", "count"),
        max_peak_signal=("peak_signal", "max"),
        mean_peak_signal=("peak_signal", "mean"),
        max_peak_qvalue=("peak_qvalue", "max")
    )
)

summary["copy_number_est"] = summary["array_len"] / summary["monomer_len"]
summary["CENH3_overlap_ratio"] = summary["CENH3_overlap_bp"] / summary["array_len"]

summary = summary.sort_values(
    ["CENH3_overlap_bp", "array_len", "copy_number_est"],
    ascending=False
)

summary.to_csv(outfile, sep="\t", index=False)

print(f"written: {outfile}")
