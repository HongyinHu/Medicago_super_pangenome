#!/bin/bash
set -euo pipefail

THREADS=20
SPECIES_LIST="genome_A17 genome_R108 genome_Msa genome_Mpo genome_474"

mkdir -p 01_jcvi
cd 01_jcvi

# 1. prepare bed and cds
for sp in ${SPECIES_LIST}
do
    echo "[INFO] Preparing ${sp}"

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

# 2. pairwise synteny
pairs=(
"genome_A17 genome_R108"
"genome_R108 genome_Msa"
"genome_Msa genome_Mpo"
"genome_Mpo genome_474"
)


for pair in "${pairs[@]}"
do
    q=$(echo $pair | cut -d' ' -f1)
    s=$(echo $pair | cut -d' ' -f2)

    echo "[INFO] Running JCVI ortholog: ${q} vs ${s}"

    python -m jcvi.compara.catalog ortholog \
        ${q} ${s} \
        --no_strip_names \
        --cpus=${THREADS}

    python -m jcvi.compara.synteny screen \
        --minspan=30 \
        --simple \
        ${q}.${s}.anchors \
        ${q}.${s}.anchors.new
done

echo "[INFO] Done. Now prepare seqids and layout, then run:"
echo "python -m jcvi.graphics.karyotype seqids layout --format pdf"