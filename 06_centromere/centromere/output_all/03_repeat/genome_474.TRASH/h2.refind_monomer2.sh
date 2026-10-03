#先定义Chr6 CENH3 核心区
cat > Chr6.CENH3_core.bed <<EOF
Chr6	93000000	105000000	Chr6_CENH3_core
EOF


TRASH_BED=genome_474.TRASH_arrays.bed

bedtools intersect \
  -a ${TRASH_BED} \
  -b Chr6.CENH3_core.bed \
  -wa \
| awk 'BEGIN{OFS="\t"}
{
    len=$5;
    count[len]++;
    total_array_len[len]+=$6;
    total_bp[len]+=$3-$2;
}
END{
    print "monomer_len","array_count","total_array_len","total_interval_bp";
    for(i in count){
        print i,count[i],total_array_len[i],total_bp[i];
    }
}' \
| sort -k3,3nr \
> Chr6.CENH3_core.TRASH_monomer_summary.tsv  #染色体


awk 'BEGIN{OFS="\t"}
$1=="Chr6" && $3>=90000000 && $2<=105000000 {
    print
}' ${TRASH_BED} \
| sort -k6,6nr \
| head -n 50 > Chr6.CENH3_core.TRASH_monomer_summary.txt
