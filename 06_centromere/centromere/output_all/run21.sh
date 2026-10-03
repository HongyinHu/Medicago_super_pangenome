#直接从 bigWig 里提取高信号区间，再合并成 CENH3 functional region

#!/bin/bash
set -euo pipefail

SP=genome_Msa

GENOME=00_genome/${SP}.fa
BW=path/to/project/10.centromere_analysis/output_all/02_chipseq/genome_Msa/05_chipseq_signal/genome_Msa_CENH3_vs_Input.q0.rmdup.50000.smooth250000.log2.bw

OUTDIR=04_cenh3_signal/${SP}.functional_centromere_from_bigwig
mkdir -p ${OUTDIR}

# 根据你的 log2 CENH3/Input 信号调整
SIGNAL_THRESHOLD=0.5

# subpeak 之间小于这个距离就合并
MERGE_GAP=1000000

# 最小功能区长度
MIN_REGION_LEN=50000

# 每条染色体保留几个候选区
TOP_N=4

samtools faidx ${GENOME}
cut -f1,2 ${GENOME}.fai > ${OUTDIR}/${SP}.genome.sizes

echo "[1] Convert bigWig to bedGraph"

bigWigToBedGraph \
  ${BW} \
  ${OUTDIR}/${SP}.CENH3_signal.bedGraph

echo "[2] Select high CENH3 signal bins"

awk -v cutoff=${SIGNAL_THRESHOLD} 'BEGIN{OFS="\t"}
$4>=cutoff {
    print $1,$2,$3,$4
}' ${OUTDIR}/${SP}.CENH3_signal.bedGraph \
| sort -k1,1 -k2,2n \
> ${OUTDIR}/${SP}.CENH3_signal.high_bins.bed

echo "[3] Merge nearby high-signal bins into candidate domains"

bedtools merge \
  -i ${OUTDIR}/${SP}.CENH3_signal.high_bins.bed \
  -d ${MERGE_GAP} \
  -c 4 \
  -o count,max,mean \
| awk 'BEGIN{OFS="\t"}
{
    len=$3-$2;
    bin_count=$4;
    max_signal=$5;
    mean_signal=$6;
    score=len*mean_signal;

    print $1,$2,$3,bin_count,max_signal,mean_signal,len,score;
}' \
> ${OUTDIR}/${SP}.CENH3_signal.cluster.all.bed

echo "[4] Filter candidate domains"

awk -v minlen=${MIN_REGION_LEN} 'BEGIN{OFS="\t"}
$7>=minlen {
    print $0
}' ${OUTDIR}/${SP}.CENH3_signal.cluster.all.bed \
> ${OUTDIR}/${SP}.CENH3_signal.cluster.filtered.bed

echo "[5] Keep top ${TOP_N} domains per chromosome"

sort -k1,1V -k8,8nr ${OUTDIR}/${SP}.CENH3_signal.cluster.filtered.bed \
| awk -v topn=${TOP_N} 'BEGIN{OFS="\t"}
{
    chr=$1;
    rank[chr]++;

    if(rank[chr]<=topn){
        if(rank[chr]==1){
            type="primary_CENH3_core";
        } else {
            type="secondary_CENH3_peak";
        }

        name=chr"_"type"_"rank[chr];

        print $1,$2,$3,name,type,rank[chr],$4,$5,$6,$7,$8;
    }
}' \
| sort -k1,1V -k2,2n \
> ${OUTDIR}/${SP}.CENH3.bigwig.top${TOP_N}.clusters.detail.bed

awk 'BEGIN{OFS="\t"}{
    print $1,$2,$3,$4
}' ${OUTDIR}/${SP}.CENH3.bigwig.top${TOP_N}.clusters.detail.bed \
> ${OUTDIR}/${SP}.CENH3.bigwig.top${TOP_N}.clusters.bed

echo "[6] Make summary table"

awk 'BEGIN{OFS="\t"}
NR==FNR{
    chr_len[$1]=$2;
    next;
}
{
    chr=$1;
    start=$2;
    end=$3;
    name=$4;
    type=$5;
    rank=$6;
    bin_count=$7;
    max_signal=$8;
    mean_signal=$9;
    size_bp=end-start;
    size_mb=size_bp/1000000;
    ratio=size_bp/chr_len[chr]*100;
    score=$11;

    print chr,chr_len[chr],start,end,size_mb,ratio,type,rank,bin_count,max_signal,mean_signal,score;
}' ${OUTDIR}/${SP}.genome.sizes \
   ${OUTDIR}/${SP}.CENH3.bigwig.top${TOP_N}.clusters.detail.bed \
| awk 'BEGIN{
    OFS="\t";
    print "Chr","Chr_length_bp","Start_bp","End_bp","Size_Mb","Ratio_percent","Type","Rank","Bin_count","Max_signal","Mean_signal","Score";
}
{print}' \
> ${OUTDIR}/${SP}.CENH3.bigwig.top${TOP_N}.clusters.summary.tsv

echo "Done."
echo "BED:"
echo "${OUTDIR}/${SP}.CENH3.bigwig.top${TOP_N}.clusters.bed"
echo "Summary:"
echo "${OUTDIR}/${SP}.CENH3.bigwig.top${TOP_N}.clusters.summary.tsv"