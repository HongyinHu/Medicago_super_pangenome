#make_R108_density_tracks.sh
#Repeat Density 和 Gene Density 最好都做成 100 kb window 的 bigWig

#!/bin/bash
set -euo pipefail

SP=genome_474


GENOME=00_genome/${SP}.fa
GENE_GFF=01_gff/${SP}.gff3  # 改成你的基因注释 GFF3
REPEAT_GFF=03_repeat/${SP}.EDTA/${SP}.fa.mod.EDTA.TEanno.gff3  # 改成你的重复序列注释 GFF3

OUTDIR=05_centromere_define/${SP}.density_tracks
WIN=100000

mkdir -p ${OUTDIR}

echo "[1] Prepare chromosome sizes"
samtools faidx ${GENOME}
cut -f1,2 ${GENOME}.fai > ${OUTDIR}/${SP}.chrom.sizes

echo "[2] Make ${WIN} bp windows"
bedtools makewindows \
  -g ${OUTDIR}/${SP}.chrom.sizes \
  -w ${WIN} \
  > ${OUTDIR}/${SP}.${WIN}bp.windows.bed

echo "[3] Convert repeat GFF to BED"

awk 'BEGIN{OFS="\t"}
$0 !~ /^#/ && NF >= 5 {
    if ($4 < $5) {
        print $1, $4-1, $5
    }
}' ${REPEAT_GFF} \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> ${OUTDIR}/${SP}.repeats.merged.bed

echo "[4] Convert gene GFF to BED"

awk 'BEGIN{OFS="\t"}
$0 !~ /^#/ && $3=="gene" {
    print $1, $4-1, $5
}' ${GENE_GFF} \
| sort -k1,1 -k2,2n \
> ${OUTDIR}/${SP}.genes.bed

echo "[5] Calculate repeat density"

bedtools coverage \
  -a ${OUTDIR}/${SP}.${WIN}bp.windows.bed \
  -b ${OUTDIR}/${SP}.repeats.merged.bed \
| awk 'BEGIN{OFS="\t"}{
    print $1,$2,$3,$7
}' \
> ${OUTDIR}/${SP}.repeat_density.${WIN}bp.bedgraph

echo "[6] Calculate gene density"

bedtools intersect \
  -a ${OUTDIR}/${SP}.${WIN}bp.windows.bed \
  -b ${OUTDIR}/${SP}.genes.bed \
  -c \
| awk -v win=${WIN} 'BEGIN{OFS="\t"}{
    gene_per_Mb = $4 / (win / 1000000);
    print $1,$2,$3,gene_per_Mb
}' \
> ${OUTDIR}/${SP}.gene_density.${WIN}bp.bedgraph

echo "[7] Convert bedGraph to bigWig"

bedGraphToBigWig \
  ${OUTDIR}/${SP}.repeat_density.${WIN}bp.bedgraph \
  ${OUTDIR}/${SP}.chrom.sizes \
  ${OUTDIR}/${SP}.repeat_density.${WIN}bp.bw

bedGraphToBigWig \
  ${OUTDIR}/${SP}.gene_density.${WIN}bp.bedgraph \
  ${OUTDIR}/${SP}.chrom.sizes \
  ${OUTDIR}/${SP}.gene_density.${WIN}bp.bw

echo "Done."
echo "Repeat density:"
echo "${OUTDIR}/${SP}.repeat_density.${WIN}bp.bw"
echo "Gene density:"
echo "${OUTDIR}/${SP}.gene_density.${WIN}bp.bw"