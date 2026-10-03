TRASH_BED=genome_474.TRASH_arrays.bed
START=90000000
END=120000000

awk -v start="$START" -v end="$END" 'BEGIN{OFS="\t"}
$1=="Chr2" && $3>=start && $2<=end {
    print
}' ${TRASH_BED} \
| sort -k6,6nr \
| head -n 50

echo "========================================="

awk -v start="$START" -v end="$END" 'BEGIN{OFS="\t"}
$1=="Chr2" && $3>=start && $2<=end {
    len=$5;
    count[len]++;
    array_bp[len]+=$6;
    overlap_bp[len]+=($3<$e?$3:e)-($2>s?$2:s);
}
END{
    print "monomer_len","array_count","total_array_len";
    for(len in count){
        print len,count[len],array_bp[len];
    }
}' s=100000000 e=115000000 ${TRASH_BED} \
| sort -k3,3nr \
| head -n 30