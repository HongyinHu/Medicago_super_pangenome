#计算甲基化水平
#切换到env环境pb_methylation

SP=genome_474

GENOME=00_genome/${SP}.fa
HIFI=00_genome/${SP}.hifi_reads.bam

OUTDIR=07_methylation/${SP}.hifi_methylation
mkdir -p ${OUTDIR}

samtools faidx ${GENOME}

pbmm2 index \
  ${GENOME} \
  ${OUTDIR}/${SP}.mmi

pbmm2 align \
  --preset HiFi \
  --sort \
  -j 32 \
  ${OUTDIR}/${SP}.mmi \
  ${HIFI} \
  ${OUTDIR}/${SP}.hifi.MMML.pbmm2.bam

samtools index -@ 8 ${OUTDIR}/${SP}.hifi.MMML.pbmm2.bam


samtools view ${OUTDIR}/${SP}.hifi.MMML.pbmm2.bam | head -n 1000 | \
awk '
/MM:Z:/ {mm++}
/ML:B:C/ {ml++}
END{
  print "MM tags:", mm+0
  print "ML tags:", ml+0
}' 


aligned_bam_to_cpg_scores \
  --bam ${OUTDIR}/${SP}.hifi.MMML.pbmm2.bam \
  --output-prefix ${OUTDIR}/${SP}.CpG \
  --threads 32