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

mkdir -p ${CLEAN_DIR} ${MAP_DIR} ${QC_DIR} ${BW_DIR} ${PEAK_DIR}

echo "[0] Check input files"

for f in ${GENOME} ${CENH3_R1} ${CENH3_R2} ${INPUT_R1} ${INPUT_R2}
do
    if [ ! -s ${f} ]; then
        echo "ERROR: missing file: ${f}"
        exit 1
    fi
done

echo "[1] Index genome"

samtools faidx ${GENOME}

if [ ! -f ${GENOME}.0123 ]; then
    bwa-mem2 index ${GENOME}
fi

############################################################
# Function for one sample
############################################################

process_sample () {

    TYPE=$1
    R1=$2
    R2=$3

    echo "============================================================"
    echo "Processing ${TYPE}"
    echo "============================================================"

    echo "[2] fastp clean ${TYPE}"

    fastp \
      -i ${R1} \
      -I ${R2} \
      -o ${CLEAN_DIR}/${PREFIX}_${TYPE}.clean.R1.fq.gz \
      -O ${CLEAN_DIR}/${PREFIX}_${TYPE}.clean.R2.fq.gz \
      -h ${QC_DIR}/${PREFIX}_${TYPE}.fastp.html \
      -j ${QC_DIR}/${PREFIX}_${TYPE}.fastp.json \
      -w ${THREADS}

    echo "[3] Map ${TYPE} reads"

    bwa-mem2.avx2 mem \
      -t ${THREADS} \
      ${GENOME} \
      ${CLEAN_DIR}/${PREFIX}_${TYPE}.clean.R1.fq.gz \
      ${CLEAN_DIR}/${PREFIX}_${TYPE}.clean.R2.fq.gz \
    | samtools sort \
      -@ ${THREADS} \
      -n \
      -o ${MAP_DIR}/${PREFIX}_${TYPE}.map.namesort.bam

    echo "[4] Add mate information"

    samtools fixmate \
      -@ ${THREADS} \
      -m \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.namesort.bam \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.fixmate.bam

    echo "[5] Coordinate sort"

    samtools sort \
      -@ ${THREADS} \
      -o ${MAP_DIR}/${PREFIX}_${TYPE}.map.positionsort.bam \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.fixmate.bam

    samtools index \
      -@ ${THREADS} \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.positionsort.bam

    echo "[6] Remove PCR duplicates"

    samtools markdup \
      -@ ${THREADS} \
      -r \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.positionsort.bam \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.rmdup.raw.bam

    samtools index \
      -@ ${THREADS} \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.rmdup.raw.bam

    echo "[7] Filter by MAPQ"

    samtools view \
      -@ ${THREADS} \
      -b \
      -q ${MAPQ} \
      -F 3844 \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.rmdup.raw.bam \
    > ${MAP_DIR}/${PREFIX}_${TYPE}.map.q${MAPQ}.rmdup.bam

    samtools index \
      -@ ${THREADS} \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.q${MAPQ}.rmdup.bam

    echo "[8] QC statistics for ${TYPE}"

    samtools flagstat \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.positionsort.bam \
    > ${QC_DIR}/${PREFIX}_${TYPE}.map.flagstat.txt

    samtools flagstat \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.rmdup.raw.bam \
    > ${QC_DIR}/${PREFIX}_${TYPE}.map.rmdup.raw.flagstat.txt

    samtools flagstat \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.q${MAPQ}.rmdup.bam \
    > ${QC_DIR}/${PREFIX}_${TYPE}.map.q${MAPQ}.rmdup.flagstat.txt

    samtools idxstats \
      ${MAP_DIR}/${PREFIX}_${TYPE}.map.q${MAPQ}.rmdup.bam \
    > ${QC_DIR}/${PREFIX}_${TYPE}.map.q${MAPQ}.rmdup.idxstats.txt

    total_before_rmdup=$(samtools view -c -F 4 ${MAP_DIR}/${PREFIX}_${TYPE}.map.positionsort.bam)
    total_after_rmdup=$(samtools view -c -F 4 ${MAP_DIR}/${PREFIX}_${TYPE}.map.rmdup.raw.bam)
    total_after_filter=$(samtools view -c -F 4 ${MAP_DIR}/${PREFIX}_${TYPE}.map.q${MAPQ}.rmdup.bam)

    awk -v a=${total_before_rmdup} -v b=${total_after_rmdup} -v c=${total_after_filter} \
    'BEGIN{
        print "sample\ttotal_mapped_before_rmdup\ttotal_mapped_after_rmdup\ttotal_mapped_after_q_filter\tduplication_rate";
        print "'${TYPE}'\t"a"\t"b"\t"c"\t"1-b/a;
    }' > ${QC_DIR}/${PREFIX}_${TYPE}.map.rmdup.summary.tsv

    echo "[9] Generate ${TYPE} coverage bigWig"

    bamCoverage \
      -b ${MAP_DIR}/${PREFIX}_${TYPE}.map.q${MAPQ}.rmdup.bam \
      --binSize ${BIN_SIZE} \
      --smoothLength ${SMOOTH_LENGTH} \
      --normalizeUsing CPM \
      --extendReads \
      -p ${THREADS} \
      -o ${BW_DIR}/${PREFIX}_${TYPE}.map.rmdup.bw

    echo "[10] Remove temporary files for ${TYPE}"

    rm -f ${MAP_DIR}/${PREFIX}_${TYPE}.map.namesort.bam
    rm -f ${MAP_DIR}/${PREFIX}_${TYPE}.map.fixmate.bam
    rm -f ${MAP_DIR}/${PREFIX}_${TYPE}.map.positionsort.bam
    rm -f ${MAP_DIR}/${PREFIX}_${TYPE}.map.positionsort.bam.bai
}

############################################################
# Process CENH3 and Input
############################################################

process_sample CENH3 ${CENH3_R1} ${CENH3_R2}
process_sample Input ${INPUT_R1} ${INPUT_R2}

############################################################
# ChIP/Input enrichment track
############################################################

echo "[11] Generate CENH3/Input log2 enrichment bigWig"

bamCompare \
  -b1 ${MAP_DIR}/${PREFIX}_CENH3.map.q${MAPQ}.rmdup.bam \
  -b2 ${MAP_DIR}/${PREFIX}_Input.map.q${MAPQ}.rmdup.bam \
  --operation log2 \
  --binSize ${BIN_SIZE} \
  --smoothLength ${SMOOTH_LENGTH} \
  --pseudocount 1 \
  --skipZeroOverZero \
  -p ${THREADS} \
  -o ${BW_DIR}/${PREFIX}_CENH3_vs_Input.q${MAPQ}.rmdup.${BIN_SIZE}.smooth${SMOOTH_LENGTH}.log2.bw

############################################################
# MACS3 broad peak calling
############################################################

echo "[12] MACS3 broad peak calling"

GENOME_SIZE=$(awk '{sum+=$2} END{print sum}' ${GENOME}.fai)

macs3 callpeak \
  -t ${MAP_DIR}/${PREFIX}_CENH3.map.q${MAPQ}.rmdup.bam \
  -c ${MAP_DIR}/${PREFIX}_Input.map.q${MAPQ}.rmdup.bam \
  -f BAMPE \
  -g ${GENOME_SIZE} \
  --broad \
  --broad-cutoff 0.01 \
  -q 0.01 \
  --keep-dup all \
  -n ${PREFIX}_CENH3.q${MAPQ}.rmdup.strict \
  --outdir ${PEAK_DIR}

echo "[13] Filter MACS3 broadPeak"

awk 'BEGIN{OFS="\t"} $9>=5 && $7>=2 {print}' \
  ${PEAK_DIR}/${PREFIX}_CENH3.q${MAPQ}.rmdup.strict_peaks.broadPeak \
> ${PEAK_DIR}/${PREFIX}_CENH3.q${MAPQ}.rmdup.strict.q5.signal2.broadPeak.bed

############################################################
# FRiP
############################################################

echo "[14] Calculate FRiP"

PEAK_BED=${PEAK_DIR}/${PREFIX}_CENH3.q${MAPQ}.rmdup.strict.q5.signal2.broadPeak.bed
CHIP_BAM=${MAP_DIR}/${PREFIX}_CENH3.map.q${MAPQ}.rmdup.bam

reads_in_peaks=$(bedtools intersect -a ${CHIP_BAM} -b ${PEAK_BED} -u | samtools view -c -)
total_reads=$(samtools view -c -F 4 ${CHIP_BAM})

awk -v a=${reads_in_peaks} -v b=${total_reads} \
'BEGIN{
    print "reads_in_peaks\t"a;
    print "total_reads\t"b;
    print "FRiP\t"a/b;
}' > ${QC_DIR}/${PREFIX}_CENH3.FRiP.txt

echo "============================================================"
echo "Done."
echo "Final CENH3 BAM:"
echo "${MAP_DIR}/${PREFIX}_CENH3.map.q${MAPQ}.rmdup.bam"
echo "Final Input BAM:"
echo "${MAP_DIR}/${PREFIX}_Input.map.q${MAPQ}.rmdup.bam"
echo "CENH3 only bigWig:"
echo "${BW_DIR}/${PREFIX}_CENH3.map.rmdup.bw"
echo "Input only bigWig:"
echo "${BW_DIR}/${PREFIX}_Input.map.rmdup.bw"
echo "CENH3/Input enrichment bigWig:"
echo "${BW_DIR}/${PREFIX}_CENH3_vs_Input.q${MAPQ}.rmdup.${BIN_SIZE}.smooth${SMOOTH_LENGTH}.log2.bw"
echo "Filtered broadPeak:"
echo "${PEAK_DIR}/${PREFIX}_CENH3.q${MAPQ}.rmdup.strict.q5.signal2.broadPeak.bed"
echo "============================================================"
