SP=genome_Msa
MONOMER=Sat215_like
GENOME=${SP}.fa

QUERY=${SP}_${MONOMER}.fa
OUTDIR=${SP}.satellite_bed

mkdir -p ${OUTDIR}

makeblastdb \
  -in ${GENOME} \
  -dbtype nucl \
  -out ${OUTDIR}/${SP}.blastdb

blastn \
  -query ${QUERY} \
  -db ${OUTDIR}/${SP}.blastdb \
  -task blastn-short \
  -word_size 7 \
  -dust no \
  -soft_masking false \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qlen sstart send evalue bitscore" \
  -num_threads 16 \
  -out ${OUTDIR}/${SP}_${MONOMER}_vs_genome.blastn.tsv


awk -v monomer="${MONOMER}" 'BEGIN{OFS="\t"}
{
    cov=$4/$5;

    if($3>=80 && cov>=0.70){
        start=($6<$7?$6:$7)-1;
        end=($6>$7?$6:$7);
        if(start<end) print $2,start,end,monomer;
    }
}' ${OUTDIR}/${SP}_${MONOMER}_vs_genome.blastn.tsv \
| sort -k1,1 -k2,2n \
> ${OUTDIR}/${SP}_${MONOMER}.80_70.bed

awk -v monomer="${MONOMER}" 'BEGIN{OFS="\t"}
{
    cov=$4/$5;

    if($3>=70 && cov>=0.60){
        start=($6<$7?$6:$7)-1;
        end=($6>$7?$6:$7);
        if(start<end) print $2,start,end,monomer;
    }
}' ${OUTDIR}/${SP}_${MONOMER}_vs_genome.blastn.tsv \
| sort -k1,1 -k2,2n \
> ${OUTDIR}/${SP}_${MONOMER}.70_60.bed