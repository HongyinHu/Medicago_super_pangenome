#使用HIC数据分析A/B区室和TAD
#使用需要切换到env环境HIC_AB_TAD

#创建基因组索引
SP=genome_474

GENOME=00_genome/${SP}.fa
GENE_DENSITY=05_centromere_define/${SP}.density_tracks/${SP}.gene_density.100000bp.bedgraph  #run.sh 设置windows
HIC_R1=00_genome/${SP}.HiC_R1.fq.gz
HIC_R2=00_genome/${SP}.HiC_R2.fq.gz
OUTDIR=08_hic_AB_TAD/${SP}

mkdir -p ${OUTDIR}

samtools faidx ${GENOME}
cut -f1,2 ${GENOME}.fai > ${OUTDIR}/${SP}.chrom.sizes

bwa-mem2 index ${GENOME}

#HIC reads 比对并生成valid pairs
bwa-mem2 mem -5SP -t 32 ${GENOME} ${HIC_R1} ${HIC_R2} | \
pairtools parse \
  --chroms-path ${OUTDIR}/${SP}.chrom.sizes \
  --min-mapq 30 \
  --walks-policy 5unique \
  --drop-sam \
  --drop-seq | \
pairtools sort --nproc 16 | \
pairtools dedup \
  --output-stats ${OUTDIR}/${SP}.dedup.stats.txt | \
pairtools select 'pair_type=="UU"' | \
bgzip -c > ${OUTDIR}/${SP}.validPairs.pairs.gz

pairix ${OUTDIR}/${SP}.validPairs.pairs.gz

#构建 100 kb 和 40 kb Hi-C matrix
for RES in 100000 40000
do
  cooler cload pairs \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    ${OUTDIR}/${SP}.chrom.sizes:${RES} \
    ${OUTDIR}/${SP}.validPairs.pairs.gz \
    ${OUTDIR}/${SP}.${RES}.raw.cool

  cp ${OUTDIR}/${SP}.${RES}.raw.cool ${OUTDIR}/${SP}.${RES}.cool

  cooler balance \
    -p 16 \
    ${OUTDIR}/${SP}.${RES}.cool
done

