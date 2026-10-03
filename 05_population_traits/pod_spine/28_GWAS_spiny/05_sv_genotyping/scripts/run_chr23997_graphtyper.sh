#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/05_sv_genotyping"
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
SAMTOOLS=path/to/home/anaconda3/envs/panpop/bin/samtools
BGZIP=path/to/home/anaconda3/envs/panpop/bin/bgzip
TABIX=path/to/home/anaconda3/envs/panpop/bin/tabix
GRAPHTYPER=path/to/home/anaconda3/envs/graphtyper/bin/graphtyper
BAM_LIST=${BAM_LIST:-$ROOT/01_input_audit/summary/bam.list}
THREADS=${THREADS:-16}
RUN_LABEL=${RUN_LABEL:-all144}
OUT="$STAGE/target_Chr23997/$RUN_LABEL"
CANDIDATE_RAW="$STAGE/target_Chr23997/Chr23997.intron2.DEL221.candidate.vcf"
CANDIDATE="$CANDIDATE_RAW.gz"

mkdir -p "$STAGE/scripts" "$STAGE/logs" "$STAGE/summary" "$STAGE/target_Chr23997" "$OUT"

ref_sequence=$($SAMTOOLS faidx "$REF" Chr4:88980831-88981052 | tail -n +2 | tr -d '\n')
ref_base=${ref_sequence:0:1}
if [[ ${#ref_sequence} -ne 222 || ${#ref_base} -ne 1 ]]; then
  echo "Unable to retrieve the 222-bp anchored reference allele at Chr4:88980831-88981052" >&2
  exit 1
fi

{
  printf '##fileformat=VCFv4.2\n'
  printf '##source=validated_Chr23997_intron2_candidate\n'
  printf '##contig=<ID=Chr4,length=103542902>\n'
  printf '##ALT=<ID=DEL,Description="Deletion">\n'
  printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">\n'
  printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">\n'
  printf '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'Chr4\t88980831\tChr23997_intron2_DEL221\t%s\t%s\t.\tPASS\tSVTYPE=DEL;END=88981052;SVLEN=-221\n' "$ref_sequence" "$ref_base"
} > "$CANDIDATE_RAW"
"$BGZIP" -f -c "$CANDIDATE_RAW" > "$CANDIDATE"
"$TABIX" -f -p vcf "$CANDIDATE"

rm -f "$OUT/Chr23997_graphtyper.done"
"$GRAPHTYPER" genotype_sv "$REF" "$CANDIDATE" \
  --sams "$BAM_LIST" \
  --region Chr4:88970000-88992000 \
  --threads "$THREADS" \
  --output "$OUT"

touch "$OUT/Chr23997_graphtyper.done"
printf 'run_label\t%s\n' "$RUN_LABEL" > "$STAGE/summary/Chr23997.${RUN_LABEL}.status.tsv"
printf 'bam_count\t%s\n' "$(wc -l < "$BAM_LIST")" >> "$STAGE/summary/Chr23997.${RUN_LABEL}.status.tsv"
printf 'candidate\tChr4:88980831-88981052 DEL 221bp\n' >> "$STAGE/summary/Chr23997.${RUN_LABEL}.status.tsv"
printf 'completed\t%s\n' "$(date '+%F %T %Z')" >> "$STAGE/summary/Chr23997.${RUN_LABEL}.status.tsv"
echo "Chr23997 Graphtyper run completed: $OUT"
