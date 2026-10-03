#只画主峰附近的trf repeat

#!/bin/bash
set -euo pipefail

SP=genome_474
CHR=Chr4
START=80500000
END=92500000 

GENOME=00_genome/${SP}.fa
OUTDIR=05_centromere_define/${SP}.TRF_peak/${CHR}_${START}_${END}

mkdir -p ${OUTDIR}

echo "[1] Index genome"
samtools faidx ${GENOME}

echo "[2] Extract CENH3 peak region"
samtools faidx ${GENOME} ${CHR}:${START}-${END} > ${OUTDIR}/${SP}.${CHR}_${START}_${END}.fa

echo "[3] Run TRF"
cd ${OUTDIR}

trf ${SP}.${CHR}_${START}_${END}.fa 2 7 7 80 10 50 2000 -d -h

DAT=$(ls *.dat)

echo "[4] Convert TRF dat to genome-coordinate BED"

awk -v chr=${CHR} -v offset=${START} 'BEGIN{OFS="\t"}
$1 ~ /^[0-9]+$/ {
    trf_start=$1;
    trf_end=$2;
    period=$3;
    copy_number=$4;
    consensus_size=$5;
    perc_match=$6;
    perc_indel=$7;
    score=$8;
    entropy=$13;

    bed_start = offset + trf_start - 2;
    bed_end   = offset + trf_end - 1;
    array_len = trf_end - trf_start + 1;

    name = "TRF_" period "bp_len" array_len "_copy" copy_number;

    print chr, bed_start, bed_end, name, score, ".", period, copy_number, consensus_size, perc_match, perc_indel, array_len, entropy;
}' ${DAT} \
| sort -k1,1 -k2,2n \
> ${SP}.${CHR}_${START}_${END}.TRF.all.bed

echo "[5] Filter TRF repeats for plotting"

# 局部图建议保留 array length >=100 bp 且 score >=50
awk 'BEGIN{OFS="\t"} $12>=100 && $5>=50 {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13
}' ${SP}.${CHR}_${START}_${END}.TRF.all.bed \
> ${SP}.${CHR}_${START}_${END}.TRF.len100.score50.bed

# pyGenomeTracks 用简化 BED
awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,$4
}' ${SP}.${CHR}_${START}_${END}.TRF.len100.score50.bed \
> ${SP}.${CHR}_${START}_${END}.TRF.plot.bed


#如果trf repeat 太密，过滤条件严格一点
awk 'BEGIN{OFS="\t"} $12>=300 && $5>=100 {print}' \
  ${SP}.${CHR}_${START}_${END}.TRF.all.bed \
> ${SP}.${CHR}_${START}_${END}.TRF.len300.score100.bed


echo "Done."
echo "Output:"
echo "${OUTDIR}/${SP}.${CHR}_${START}_${END}.TRF.plot.bed"