#!/usr/bin/env bash
set -euo pipefail

RUN=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$RUN/08_Chr23997_targeted_genotyping_20260715"
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
SMOOVE=path/to/home/anaconda3/envs/smoove/bin/smoove
DELLY=path/to/home/anaconda3/envs/delly/bin/delly
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python

export PATH=path/to/home/anaconda3/envs/smoove/bin:path/to/home/anaconda3/envs/delly/bin:path/to/home/anaconda3/envs/panpop/bin:/usr/bin:/bin

mkdir -p "$STAGE/inputs" "$STAGE/pilot" "$STAGE/logs/pilot" "$STAGE/summary"

"$BCFTOOLS" view -Oz -o "$STAGE/inputs/Chr23997.intron2_DEL.svtyper.vcf.gz" \
    "$STAGE/inputs/Chr23997.intron2_DEL.svtyper.vcf"
"$BCFTOOLS" index -f -t "$STAGE/inputs/Chr23997.intron2_DEL.svtyper.vcf.gz"
"$BCFTOOLS" view -Ob -o "$STAGE/inputs/Chr23997.intron2_DEL.delly.bcf" \
    "$STAGE/inputs/Chr23997.intron2_DEL.delly.vcf"
"$BCFTOOLS" index -f "$STAGE/inputs/Chr23997.intron2_DEL.delly.bcf"

for sid in 109 160 203 204; do
    bam=path/to/project/38.medicago_resequence/2.call_SNP_new/04_markdup_all/bam/${sid}.dedup.bam
    outdir="$STAGE/pilot/$sid"
    mkdir -p "$outdir"
    "$SMOOVE" genotype --name "$sid" --outdir "$outdir" --fasta "$REF" \
        --removepr --processes 2 --vcf "$STAGE/inputs/Chr23997.intron2_DEL.svtyper.vcf.gz" "$bam" \
        >"$STAGE/logs/pilot/${sid}.smoove.out" 2>"$STAGE/logs/pilot/${sid}.smoove.err"
    "$DELLY" call -g "$REF" -v "$STAGE/inputs/Chr23997.intron2_DEL.delly.bcf" \
        -o "$outdir/${sid}.delly.bcf" "$bam" \
        >"$STAGE/logs/pilot/${sid}.delly.out" 2>"$STAGE/logs/pilot/${sid}.delly.err"
    "$BCFTOOLS" index -f "$outdir/${sid}.delly.bcf"
    "$PYTHON" "$STAGE/scripts/extract_chr23997_bam_evidence.py" \
        --bam "$bam" --sample "$sid" --out "$outdir/${sid}.evidence.tsv" \
        >"$STAGE/logs/pilot/${sid}.evidence.out" 2>"$STAGE/logs/pilot/${sid}.evidence.err"
done

for sid in 109 160 203 204; do
    outdir="$STAGE/pilot/$sid"
    "$BCFTOOLS" query -f '%CHROM\t%POS\t%END\t[%GT\t%GQ\t%PE\t%SR\t%SU]\n' \
        "$outdir/${sid}-smoove.genotyped.vcf.gz" | sed "s/^/${sid}\tsvtyper\t/"
    "$BCFTOOLS" query -f '%CHROM\t%POS\t%END\t[%GT\t%GQ\t%DR\t%DV\t%RR\t%RV]\n' \
        "$outdir/${sid}.delly.bcf" | sed "s/^/${sid}\tdelly\t/"
    tail -n 1 "$outdir/${sid}.evidence.tsv" | sed "s/^/${sid}\tevidence\t/"
done > "$STAGE/summary/pilot_results.raw.tsv"

date '+completed\t%F %T %Z' > "$STAGE/summary/pilot.done"
