#!/bin/bash
set -euo pipefail

############################################################
# CENH3 ChIP-seq and Input processing pipeline
############################################################

SP=genome_Msa

PREFIX=${SP}
OUTDIR=02_chipseq/${SP}
GENOME=00_genome/${SP}.fa

CENH3_R1=02_chipseq/${SP}_CENH3_R1.fq.gz
CENH3_R2=02_chipseq/${SP}_CENH3_R2.fq.gz
INPUT_R1=02_chipseq/${SP}_Input_R1.fq.gz
INPUT_R2=02_chipseq/${SP}_Input_R2.fq.gz

THREADS=24
MAPQ=20

BIN_SIZE=10000
SMOOTH_LENGTH=50000

CLEAN_DIR=${OUTDIR}/02_chipseq_clean
MAP_DIR=${OUTDIR}/03_chipseq_mapping
QC_DIR=${OUTDIR}/04_chipseq_qc
BW_DIR=${OUTDIR}/05_chipseq_signal
PEAK_DIR=${OUTDIR}/06_chipseq_peak


echo "[12] MACS3 broad peak calling"

GENOME_SIZE=$(awk '{sum+=$2} END{print sum}' ${GENOME}.fai)

for Q in 0 10 20
do
    echo "Processing MAPQ ${Q}"
    macs3 callpeak \
      -t ${MAP_DIR}/${PREFIX}_CENH3.map.q${Q}.rmdup.bam \
      -c ${MAP_DIR}/${PREFIX}_Input.map.q${Q}.rmdup.bam \
      -f BAMPE \
      -g ${GENOME_SIZE} \
      --broad \
      --broad-cutoff 0.1 \
      -q 0.05 \
      --keep-dup all \
      -n ${PREFIX}_CENH3.q${Q}.rmdup.strict \
      --outdir ${PEAK_DIR}

    echo "[13] Filter MACS3 broadPeak"

    awk 'BEGIN{OFS="\t"} $9>=5 && $7>=2 {print}' \
      ${PEAK_DIR}/${PREFIX}_CENH3.q${Q}.rmdup.strict_peaks.broadPeak \
    > ${PEAK_DIR}/${PREFIX}_CENH3.q${Q}.rmdup.strict.q5.signal2.broadPeak.bed
done
