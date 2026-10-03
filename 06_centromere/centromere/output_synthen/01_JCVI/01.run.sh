#批量修复多个物种的gff文件，使得生成的bed id和cds id 对应上
#同时获取jcvi的bed文件

species="A17 R108 Msa Mpo Cbu"

for sp in genome_A17 genome_R108 genome_Msa genome_Mpo genome_474
do

    python -m jcvi.formats.gff bed \
      --type=mRNA \
      --key=ID \
      00_genome/${sp}.gff3 \
      -o 01_JCVI/${sp}.raw.bed

    python -m jcvi.formats.fasta format \
      00_genome/${sp}.cds.fa \
      01_JCVI/${sp}.cds

    (
    cd 01_JCVI
    awk 'BEGIN{OFS="\t"}{
    sub(/%3[Bb]Parent.*/, "", $4);
    print
    }' ${sp}.raw.bed > ${sp}.bed

    rm ${sp}.raw.bed

	)
done