#做成窗口化的Methylation density

SP=genome_474

GENOME=00_genome/${SP}.fa
OUTDIR=07_methylation/${SP}.hifi_methylation

#准备chrom sizes和窗口
cut -f1,2 ${GENOME}.fai > ${OUTDIR}/${SP}.chrom.sizes

bedtools makewindows \
  -g ${OUTDIR}/${SP}.chrom.sizes \
  -w 50000 \
  > ${OUTDIR}/${SP}.50k.windows.bed


#从CpG bed 中提取甲基化位点
zcat ${OUTDIR}/${SP}.CpG.combined.bed.gz | \
awk 'BEGIN{OFS="\t"}
$0 !~ /^#/ && $6 >= 4 {
    print $1,$2,$3,$4,$6
}' \
> ${OUTDIR}/${SP}.CpG.site.methylation.cov4.bed


windows_file=${OUTDIR}/${SP}.50k.windows.bed
meth_file=${OUTDIR}/${SP}.CpG.site.methylation.cov4.bed
out_file=${OUTDIR}/${SP}.CpG_methylation_density.50k.bedGraph

export windows_file
export meth_file
export out_file


#窗口计算 coverage-weighted methylation
cat > ${OUTDIR}/weighted_methylation_window.py <<'PY'
import sys
import os
import pandas as pd
import numpy as np

windows_file = os.environ["windows_file"]
meth_file = os.environ["meth_file"]
out_file = os.environ["out_file"]

# windows: chr start end
windows = pd.read_csv(
    windows_file, sep="\t", header=None,
    names=["chr", "start", "end"]
)

meth = pd.read_csv(
    meth_file, sep="\t", header=None,
    names=["chr", "start", "end", "score", "cov"]
)

meth["bin_start"] = (meth["start"] // 50000) * 50000
meth["weighted"] = meth["score"] * meth["cov"]

summary = (
    meth.groupby(["chr", "bin_start"], as_index=False)
        .agg(weighted_sum=("weighted", "sum"),
             cov_sum=("cov", "sum"))
)

summary["value"] = summary["weighted_sum"] / summary["cov_sum"]

summary = summary.rename(columns={"bin_start": "start"})
summary["end"] = summary["start"] + 50000

result = windows.merge(
    summary[["chr", "start", "value"]],
    on=["chr", "start"],
    how="left"
)

result["value"] = result["value"].fillna(0)

result[["chr", "start", "end", "value"]].to_csv(
    out_file, sep="\t", header=False, index=False,
    float_format="%.6f"
)

print("write file to %s"%out_file)
PY

#执行脚本
python3 ${OUTDIR}/weighted_methylation_window.py

#转成bigWig
bedGraphToBigWig \
  ${OUTDIR}/${SP}.CpG_methylation_density.50k.bedGraph \
  ${OUTDIR}/${SP}.chrom.sizes \
  ${OUTDIR}/${SP}.CpG_methylation_density.50k.bw


