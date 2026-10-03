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

echo "================================="

#专门查找长阵列
awk 'BEGIN{OFS="\t"}
$1 ~ /^[0-9]+$/ {
    len=$2-$1+1;
    if(len>=10000 || $4>=50 || $8>=1000){
        print $1,$2,$3,$4,$5,$6,$7,$8,len,$NF;
    }
}' ${DAT} | column -t | head -n 50



# 将TRF结果转为BED
DAT=Chr6_93_100.5Mb.fa.2.7.7.80.10.50.2000.dat
CHR=Chr6
REGION_START=93000000

awk -v chr=${CHR} -v offset=${REGION_START} 'BEGIN{OFS="\t"}
$1 ~ /^[0-9]+$/ {
    bed_start = offset + $1 - 2;
    bed_end   = offset + $2 - 1;
    len       = $2 - $1 + 1;
    name      = "TRF_"$3"bp_copy"$4"_len"len;
    print chr,bed_start,bed_end,name,$3,len,$4,$6,$7,$8;
}' ${DAT} \
> Chr6_93_100.5Mb.TRF.bed

