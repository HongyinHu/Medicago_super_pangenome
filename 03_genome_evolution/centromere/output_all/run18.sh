#使用HIC数据分析A/B区室和TAD
#使用需要切换到env环境hic_ab_tad
#本脚本是TAD分析
#调用 TAD domains

#创建基因组索引
SP=genome_474

GENOME=00_genome/${SP}.fa
GENE_DENSITY=05_centromere_define/${SP}.density_tracks/${SP}.gene_density.100000bp.bedgraph
HIC_R1=00_genome/${SP}.HiC_R1.fq.gz
HIC_R2=00_genome/${SP}.HiC_R2.fq.gz
OUTDIR=08_hic_AB_TAD/${SP}

#TAD分析,调用 TAD domains

#先把cooler转为h5
hicConvertFormat \
  -m ${OUTDIR}/${SP}.100000.cool \
  --inputFormat cool \
  --outputFormat h5 \
  -o ${OUTDIR}/${SP}.100000.h5

#调用 TAD：
hicFindTADs \
  --matrix ${OUTDIR}/${SP}.100000.h5 \
  --outPrefix ${OUTDIR}/${SP}.TAD.100kb \
  --minDepth 300000   \
  --maxDepth 1000000 \
  --step 100000 \
  --thresholdComparisons 0.01 \
  --delta 0.01 \
  --correctForMultipleTesting fdr \
  --numberOfProcessors 16

#ini文件写入

# [TAD domains]
# file = 06_hic_AB_TAD/genome_R108/genome_R108.TAD.40kb_domains.bed
# title = TAD domains
# height = 0.35
# color = #A65E2E
# alpha = 0.4
# border_color = none
# display = collapsed
# labels = false
# file_type = bed
