FILE=genome_474_CENH3.q0.rmdup.strict.q5.signal2.broadPeak.bed

#先过滤弱峰
awk 'BEGIN{OFS="\t"}
$9>=2 && $7>=2 {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9
}' ${FILE} \
| sort -k1,1 -k2,2n \
> ${FILE}.filtered.broadPeak

#合并相邻peak为候选domain
bedtools merge \
  -i ${FILE}.filtered.broadPeak \
  -d 100000 \
  -c 4,7,9 \
  -o count,max,max \
> ${FILE}.filtered.broadPeak.merge100k.bed

#在筛选较大的peak cluster
awk 'BEGIN{OFS="\t"}
{
    len=$3-$2;
    peak_count=$4;
    max_signal=$5;
    max_q=$6;

    if(len>=100000 && peak_count>=3){
        print $1,$2,$3,"CENH3_cluster_"NR,peak_count,max_signal,max_q,len;
    }
}' ${FILE}.filtered.broadPeak.merge100k.bed \
> ${FILE}.filtered.broadPeak.domain.candidates.bed