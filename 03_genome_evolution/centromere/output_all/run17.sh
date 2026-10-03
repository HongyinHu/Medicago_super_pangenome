#使用HIC数据分析A/B区室和TAD
#使用需要切换到env环境HIC_AB_TAD环境
#本脚本是AB区室分析

#创建基因组索引
SP=genome_474

GENOME=00_genome/${SP}.fa
GENE_DENSITY=05_centromere_define/${SP}.density_tracks/${SP}.gene_density.100000bp.bedgraph  #run13.sh
HIC_R1=00_genome/${SP}.HiC_R1.fq.gz
HIC_R2=00_genome/${SP}.HiC_R2.fq.gz
OUTDIR=08_hic_AB_TAD/${SP}

#根据gene density确定A/B正负方向
{
echo -e "chrom\tstart\tend\tgene_density"
cat ${GENE_DENSITY}
} > ${OUTDIR}/${SP}.gene_density.100kb.phasing.tsv

#计算 A/B compartment eigenvector
cooltools eigs-cis \
  --phasing-track ${OUTDIR}/${SP}.gene_density.100kb.phasing.tsv::gene_density \
  --n-eigs 3 \
  --bigwig \
  -o ${OUTDIR}/${SP}.AB.100kb \
  ${OUTDIR}/${SP}.100000.cool


#如果发现正值区域反而是低基因密度区，可以把 E1 乘以 -1：
# awk 'BEGIN{OFS="\t"}
# NR==1 {print; next}
# {
#     $4 = -$4;
#     print
# }' ${OUTDIR}/${SP}.AB.100kb.cis.vecs.tsv \
# > ${OUTDIR}/${SP}.AB.100kb.cis.vecs.flipped.tsv

#然后转成 bedGraph 和 bigWig
# awk 'BEGIN{OFS="\t"}
# NR>1 {
#     print $1,$2,$3,$4
# }' ${OUTDIR}/${SP}.AB.100kb.cis.vecs.flipped.tsv \
# > ${OUTDIR}/${SP}.AB.100kb.E1.flipped.bedGraph

# bedGraphToBigWig \
#   ${OUTDIR}/${SP}.AB.100kb.E1.flipped.bedGraph \
#   ${OUTDIR}/${SP}.chrom.sizes \
#   ${OUTDIR}/${SP}.AB.100kb.E1.flipped.bw


#TAD分析
#“TAD” 这一层本质上是 Hi-C contact map 的三角热图。可以直接把 40 kb .cool 文件放进 pyGenomeTrack
#直接使用文件08_hic_AB_TAD/${SP}/${SP}.40000.cool

#ini文件写入

# [spacer]
# height = 0.15

# [TAD]
# file = 06_hic_AB_TAD/genome_R108/genome_R108.40000.cool
# title = TAD
# height = 4.5
# depth = 30000000
# transform = log1p
# min_value = 1
# max_value = 100
# colormap = OrRd
# show_masked_bins = false
# file_type = hic_matrix