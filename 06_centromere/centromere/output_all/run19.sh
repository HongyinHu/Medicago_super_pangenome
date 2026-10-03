#从 TRASH GFF 转成 array BED
#!/bin/bash
set -euo pipefail

SP=genome_R108
#CENTER_LENS=65

# for CENTER_LENS in 179 32 48 317 102 118 29 53 65
# for CENTER_LENS in 185 144 583 374 168 366 181 50 154
# for CENTER_LENS in 346 691 313 549 587 51 288 171 138
# for CENTER_LENS in 183 168 51 321 316
for CENTER_LENS in 168 59 287 132 515 307

do
echo "percess ${CENTER_LENS}======================="
TRASH_GFF=03_repeat/${SP}.TRASH/TRASH_${SP}.fa.gff
OUTDIR=05_centromere_define/${SP}.TRASH_bed
SATNAME=CentM${CENTER_LENS}

TOL=${6:-5}
MERGE_D=${7:-500}
MIN_ARRAY_LEN=${8:-2000}
MIN_COPY=${9:-5}
MIN_DENSITY=${10:-0.20}

mkdir -p ${OUTDIR}

RAW_ALL=${OUTDIR}/${SP}.TRASH.raw_intervals.bed
RAW_SAT=${OUTDIR}/${SP}.${SATNAME}.raw.tol${TOL}.bed

OUT_DETAIL=${OUTDIR}/${SP}.${SATNAME}.TRASH_array.tol${TOL}.merge${MERGE_D}.len${MIN_ARRAY_LEN}.copy${MIN_COPY}.den${MIN_DENSITY}.detail.tsv
OUT_BED=${OUTDIR}/${SP}.${SATNAME}.TRASH_array.tol${TOL}.merge${MERGE_D}.len${MIN_ARRAY_LEN}.copy${MIN_COPY}.den${MIN_DENSITY}.bed

echo "[1] Convert TRASH GFF to raw interval BED"

awk 'BEGIN{OFS="\t"}
$0 !~ /^#/ && NF >= 5 {
    len = $5 - $4 + 1;
    id = $1"_"$4"_"$5;
    print $1, $4-1, $5, id, len, len;
}' ${TRASH_GFF} \
| sort -k1,1 -k2,2n \
> ${RAW_ALL}

echo "[2] Extract target satellite family"
echo "    centers = ${CENTER_LENS}"
echo "    tolerance = +/- ${TOL} bp"

awk -v centers="${CENTER_LENS}" -v tol="${TOL}" 'BEGIN{
    OFS="\t";
    split(centers,a,",");
    for(i in a) center[i]=a[i];
}
{
    monomer_len=$5;
    keep=0;
    for(i in center){
        if(monomer_len >= center[i]-tol && monomer_len <= center[i]+tol){
            keep=1;
        }
    }
    if(keep==1){
        print $1,$2,$3,$4,$5,$6;
    }
}' ${RAW_ALL} \
| sort -k1,1 -k2,2n \
> ${RAW_SAT}

echo "[3] Merge nearby intervals into satellite arrays"
echo "    merge distance = ${MERGE_D}"
echo "    minimum array length = ${MIN_ARRAY_LEN}"
echo "    minimum copy count = ${MIN_COPY}"
echo "    minimum repeat density = ${MIN_DENSITY}"

bedtools merge \
  -i ${RAW_SAT} \
  -d ${MERGE_D} \
  -c 4,5,6 \
  -o count,distinct,sum \
| awk -v name="${SATNAME}" \
      -v minlen="${MIN_ARRAY_LEN}" \
      -v mincopy="${MIN_COPY}" \
      -v mindensity="${MIN_DENSITY}" 'BEGIN{OFS="\t"}
{
    array_len = $3 - $2;
    copy_count = $4;
    monomer_lens = $5;
    repeat_bp = $6;
    density = repeat_bp / array_len;

    if(array_len >= minlen && copy_count >= mincopy && density >= mindensity){
        print $1,$2,$3,name,copy_count,monomer_lens,array_len,repeat_bp,density;
    }
}' \
> ${OUT_DETAIL}

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' ${OUT_DETAIL} \
> ${OUT_BED}

echo "[Done]"
echo "Raw family intervals:"
echo "${RAW_SAT}"
echo "Detailed filtered array file:"
echo "${OUT_DETAIL}"
echo "pyGenomeTracks BED:"
echo "${OUT_BED}"

done