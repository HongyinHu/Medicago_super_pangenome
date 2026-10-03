#!/bin/bash
set -euo pipefail

SP=genome_R108

GENOME=data/${SP}.fa
PERI_BED=data/${SP}.pericentromere.final.bed

# 原始轨道
CENH3_BW=data/${SP}.CENH3_vs_Input.q20.10k.log2.bw
CENH3_DOMAIN=data/${SP}.CENH3.q20_peaks.broadPeak.bed
#TRF_BED=data/${SP}.TRF_peak/${SP}.TRF.plot.bed

# satellite bed，可以是多个文件
SAT_DIR=data/${SP}.satellite_bed

OUTDIR=data/${SP}.pericentromere_only_tracks
mkdir -p ${OUTDIR}

echo "[1] Genome sizes"
samtools faidx ${GENOME}
cut -f1,2 ${GENOME}.fai > ${OUTDIR}/${SP}.chrom.sizes

echo "[2] Sort pericentromere BED"
sort -k1,1 -k2,2n ${PERI_BED} \
| bedtools merge -i - \
> ${OUTDIR}/${SP}.pericentromere.merge.bed

PERI=${OUTDIR}/${SP}.pericentromere.merge.bed

echo "[3] Clip CENH3 domain to pericentromere"

awk 'BEGIN{OFS="\t"}
{
    name = (NF>=4 ? $4 : $1"_CENH3_domain")
    print $1,$2,$3,name
}' ${CENH3_DOMAIN} \
| sort -k1,1 -k2,2n \
> ${OUTDIR}/${SP}.CENH3.domain.tmp.bed

bedtools intersect \
  -a ${OUTDIR}/${SP}.CENH3.domain.tmp.bed \
  -b ${PERI} \
  -wa -wb \
| awk 'BEGIN{OFS="\t"}
{
    s=($2>$6?$2:$6);
    e=($3<$7?$3:$7);
    if(s<e) print $1,s,e,$4;
}' \
| sort -k1,1 -k2,2n \
> ${OUTDIR}/${SP}.CENH3.domain.pericentromere_only.bed

# echo "[4] Clip TRF repeats to pericentromere"# 

# cut -f1-4 ${TRF_BED} \
# | sort -k1,1 -k2,2n \
# > ${OUTDIR}/${SP}.TRF.tmp.bed# 

# bedtools intersect \
#   -a ${OUTDIR}/${SP}.TRF.tmp.bed \
#   -b ${PERI} \
#   -wa -wb \
# | awk 'BEGIN{OFS="\t"}
# {
#     s=($2>$6?$2:$6);
#     e=($3<$7?$3:$7);
#     if(s<e) print $1,s,e,$4;
# }' \
# | sort -k1,1 -k2,2n \
# > ${OUTDIR}/${SP}.TRF.pericentromere_only.bed

echo "[5] Clip satellite beds to pericentromere"

mkdir -p ${OUTDIR}/satellite

for bed in ${SAT_DIR}/*.bed
do
    name=$(basename ${bed} .bed)

    cut -f1-4 ${bed} \
    | sort -k1,1 -k2,2n \
    > ${OUTDIR}/satellite/${name}.tmp.bed

    bedtools intersect \
      -a ${OUTDIR}/satellite/${name}.tmp.bed \
      -b ${PERI} \
      -wa -wb \
    | awk 'BEGIN{OFS="\t"}
    {
        s=($2>$6?$2:$6);
        e=($3<$7?$3:$7);
        if(s<e) print $1,s,e,$4;
    }' \
    | sort -k1,1 -k2,2n \
    > ${OUTDIR}/satellite/${name}.pericentromere_only.bed
done

echo "[6] Clip CENH3 bigWig to pericentromere"

TMP_BG=${OUTDIR}/${SP}.CENH3.pericentromere_only.tmp.bedGraph
OUT_BG=${OUTDIR}/${SP}.CENH3.pericentromere_only.bedGraph
OUT_BW=${OUTDIR}/${SP}.CENH3.pericentromere_only.bw

rm -f ${TMP_BG}

while read chr start end
do
    echo "Extract ${chr}:${start}-${end}"

    bigWigToBedGraph \
      -chrom=${chr} \
      -start=${start} \
      -end=${end} \
      ${CENH3_BW} \
      stdout >> ${TMP_BG}

done < ${PERI}

sort -k1,1 -k2,2n ${TMP_BG} > ${OUT_BG}

bedGraphToBigWig \
  ${OUT_BG} \
  ${OUTDIR}/${SP}.chrom.sizes \
  ${OUT_BW}

echo "Done."
echo "Output directory:"
echo "${OUTDIR}"