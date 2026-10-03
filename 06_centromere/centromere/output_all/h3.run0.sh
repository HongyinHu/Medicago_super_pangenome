#!/bin/bash
set -euo pipefail

############################################################
# CENH3 ChIP-seq and Input processing pipeline
############################################################

SP=genome_474

PREFIX=${SP}
OUTDIR=02_chipseq/${SP}
GENOME=00_genome/${SP}.fa

CENH3_R1=02_chipseq/${SP}_CENH3_R1.fq.gz
CENH3_R2=02_chipseq/${SP}_CENH3_R2.fq.gz
INPUT_R1=02_chipseq/${SP}_Input_R1.fq.gz
INPUT_R2=02_chipseq/${SP}_Input_R2.fq.gz

THREADS=24
MAPQ=20

BIN_SIZE=100000
SMOOTH_LENGTH=500000

CLEAN_DIR=${OUTDIR}/02_chipseq_clean
MAP_DIR=${OUTDIR}/03_chipseq_mapping
QC_DIR=${OUTDIR}/04_chipseq_qc
BW_DIR=${OUTDIR}/05_chipseq_signal
PEAK_DIR=${OUTDIR}/06_chipseq_peak

############################################################
# Generate CENH3/Input enrichment tracks under different MAPQ
############################################################

echo "[11] Generate CENH3/Input log2 enrichment bigWig for q0/q10/q20/q30"

for Q in 0 10 20 30
do
    echo "Processing MAPQ ${Q}"

    samtools view -@ ${THREADS} -b -q ${Q} -F 3844 \
      ${MAP_DIR}/${PREFIX}_CENH3.map.rmdup.raw.bam \
    > ${MAP_DIR}/${PREFIX}_CENH3.map.q${Q}.rmdup.bam

    samtools index -@ ${THREADS} \
      ${MAP_DIR}/${PREFIX}_CENH3.map.q${Q}.rmdup.bam

    samtools view -@ ${THREADS} -b -q ${Q} -F 3844 \
      ${MAP_DIR}/${PREFIX}_Input.map.rmdup.raw.bam \
    > ${MAP_DIR}/${PREFIX}_Input.map.q${Q}.rmdup.bam

    samtools index -@ ${THREADS} \
      ${MAP_DIR}/${PREFIX}_Input.map.q${Q}.rmdup.bam

    bamCompare \
      -b1 ${MAP_DIR}/${PREFIX}_CENH3.map.q${Q}.rmdup.bam \
      -b2 ${MAP_DIR}/${PREFIX}_Input.map.q${Q}.rmdup.bam \
      --operation log2 \
      --binSize ${BIN_SIZE} \
      --smoothLength ${SMOOTH_LENGTH} \
      --pseudocount 5 \
      --skipZeroOverZero \
      -p ${THREADS} \
      -o ${BW_DIR}/${PREFIX}_CENH3_vs_Input.q${Q}.rmdup.${BIN_SIZE}.smooth${SMOOTH_LENGTH}.log2.bw
done