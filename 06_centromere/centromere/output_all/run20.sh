#生成CENHH3功能着丝粒坐标

#!/bin/bash
set -euo pipefail

SP=genome_Mpo

GENOME=00_genome/${SP}.fa

# MACS2/MACS3 broadPeak 文件
REP1=path/to/project/10.centromere_analysis/output_all/02_chipseq/genome_Mpo/06_chipseq_peak/genome_Mpo_CENH3.q0.rmdup.strict_peaks.broadPeak

OUTDIR=04_cenh3_signal/${SP}.functional_centromere_topN
mkdir -p ${OUTDIR}

############################################################
# Parameters
############################################################

SIGNAL_CUTOFF=0
QVALUE_CUTOFF=0

# 150 kb 是澳洲稻文章类似的合并尺度
# 如果两个峰距离很近，不想被合并，可以调小
MERGE_D=500000

MIN_CLUSTER_LEN=1
MIN_PEAK_COUNT=1

# 每条染色体保留几个 CENH3 cluster
TOP_N=10

############################################################
# Genome sizes
############################################################

samtools faidx ${GENOME}
cut -f1,2 ${GENOME}.fai > ${OUTDIR}/${SP}.genome.sizes

############################################################
# 1. Filter MACS broadPeak
############################################################

echo "[1] Filter broadPeak"

awk -v sig=${SIGNAL_CUTOFF} -v qv=${QVALUE_CUTOFF} 'BEGIN{OFS="\t"}
$1 !~ /^#/ && $3>$2 && $7>=sig && $9>=qv {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9;
}' ${REP1} \
| sort -k1,1 -k2,2n \
> ${OUTDIR}/${SP}.CENH3.filtered.broadPeak.bed

############################################################
# 2. Merge nearby peaks into CENH3 clusters
############################################################

echo "[2] Merge nearby peaks into clusters"

bedtools merge \
  -i ${OUTDIR}/${SP}.CENH3.filtered.broadPeak.bed \
  -d ${MERGE_D} \
  -c 4,7,9 \
  -o count,max,max \
| awk 'BEGIN{OFS="\t"}
{
    len=$3-$2;
    peak_count=$4;
    max_signal=$5;
    max_qvalue=$6;
    score=peak_count*max_signal;

    print $1,$2,$3,peak_count,max_signal,max_qvalue,len,score;
}' \
> ${OUTDIR}/${SP}.CENH3.cluster.all.bed

############################################################
# 3. Filter clusters
############################################################

echo "[3] Filter clusters"

awk -v minlen=${MIN_CLUSTER_LEN} -v mincount=${MIN_PEAK_COUNT} 'BEGIN{OFS="\t"}
$7>=minlen && $4>=mincount {
    print $0;
}' ${OUTDIR}/${SP}.CENH3.cluster.all.bed \
> ${OUTDIR}/${SP}.CENH3.cluster.filtered.bed

############################################################
# 4. Keep top N clusters per chromosome
############################################################

echo "[4] Keep top ${TOP_N} CENH3 clusters per chromosome"

sort -k1,1V -k8,8nr ${OUTDIR}/${SP}.CENH3.cluster.filtered.bed \
| awk -v topn=${TOP_N} 'BEGIN{OFS="\t"}
{
    chr=$1;
    rank[chr]++;

    if(rank[chr] <= topn){
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
> ${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.detail.bed

############################################################
# 5. Simple BED for pyGenomeTracks
############################################################

awk 'BEGIN{OFS="\t"}{
    print $1,$2,$3,$4;
}' ${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.detail.bed \
> ${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.bed

############################################################
# 6. Primary-only functional centromere BED
############################################################

awk 'BEGIN{OFS="\t"} $5=="primary_CENH3_core" {
    print $1,$2,$3,$4;
}' ${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.detail.bed \
> ${OUTDIR}/${SP}.CENH3.primary_functional_centromere.bed

############################################################
# 7. Summary table
############################################################

echo "[5] Make summary table"

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
    peak_count=$7;
    max_signal=$8;
    max_qvalue=$9;
    size_bp=end-start;
    size_mb=size_bp/1000000;
    ratio=size_bp/chr_len[chr]*100;
    score=$11;

    print chr,chr_len[chr],start,end,size_mb,ratio,type,rank,peak_count,max_signal,max_qvalue,score;
}' ${OUTDIR}/${SP}.genome.sizes \
   ${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.detail.bed \
| awk 'BEGIN{
    OFS="\t";
    print "Chr","Chr_length_bp","Start_bp","End_bp","Size_Mb","Ratio_percent","Type","Rank","Peak_count","Max_signal","Max_qvalue","Score";
}
{print}' \
> ${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.summary.tsv

echo "Done."

echo "All top clusters:"
echo "${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.detail.bed"

echo "Simple BED for plotting:"
echo "${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.bed"

echo "Primary functional centromere only:"
echo "${OUTDIR}/${SP}.CENH3.primary_functional_centromere.bed"

echo "Summary:"
echo "${OUTDIR}/${SP}.CENH3.top${TOP_N}.clusters.summary.tsv"