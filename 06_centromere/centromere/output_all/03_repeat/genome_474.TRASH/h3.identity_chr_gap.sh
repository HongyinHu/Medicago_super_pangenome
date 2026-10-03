#检查Chr6主峰区域是否组装完整

GENOME=path/to/project/6.genome_gap_fill_new/output/genome_474.near_T2T.ctg.final.fa

samtools faidx ${GENOME}

samtools faidx ${GENOME} Chr6:90000000-105000000 | \
grep -v ">" | \
awk '{
    total += length($0);
    n += gsub(/[Nn]/,"N");
}
END{
    print "total_bp =", total;
    print "N_bp =", n;
    print "N_ratio =", n/total;
}'



samtools faidx ${GENOME} Chr6:93000000-100500000 > Chr6_93_100.5Mb.fa

trf Chr6_93_100.5Mb.fa 2 7 7 80 10 50 2000 -d -h

DAT=Chr6_93_100.5Mb.fa.2.7.7.80.10.50.2000.dat

awk 'BEGIN{OFS="\t"}
$1 ~ /^[0-9]+$/ {
    len=$2-$1+1;
    period=$3;
    n[period]++;
    total_len[period]+=len;
    total_copy[period]+=$4;
    if(len>max_len[period]) max_len[period]=len;
    if($8>max_score[period]) max_score[period]=$8;
}
END{
    print "period","array_count","total_array_len","max_array_len","mean_copy","max_score";
    for(p in n){
        print p,n[p],total_len[p],max_len[p],total_copy[p]/n[p],max_score[p];
    }
}' ${DAT} | sort -k3,3nr | column -t | head -n 30