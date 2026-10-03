#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project
WORK="$BASE/N_1.coding_gene_anno/06_mixed_illumina_expression_SV_distance_20260704"
RNADIR="$BASE/N_1.coding_gene_anno/00_inputs/my_species_RNA"
GFF="$BASE/N_4.pod_spiny/00_data/4.reference_anno/genome_Msa.gff"
FA="$BASE/N_3.call_SV/00_data/3.two_ref/genome_Msa.fa"
SVCTX="$BASE/N_3.call_SV/07.read_mapping_flankQC_panSV_Msa/results/Msa/tables/panSV_PAV_gene_context.exclude_genome_Msa.events.tsv"
PAV="$BASE/N_3.call_SV/07.read_mapping_flankQC_panSV_Msa/results/Msa/panSV/panSV.Msa.flankQC.read_primary.PAV.matrix.tsv"

GFFREAD=path/to/home/anaconda3/envs/biosofeware/bin/gffread
SALMON=path/to/home/anaconda3/envs/biosofeware/bin/salmon
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript

THREADS_INDEX=${THREADS_INDEX:-16}
THREADS_QUANT=${THREADS_QUANT:-6}
MAX_QUANT_JOBS=${MAX_QUANT_JOBS:-3}

mkdir -p "$WORK"/{scripts,config,logs,reference,salmon_index,quants,tables,figures}
LOG="$WORK/logs/run_mixed_illumina_expression_SV_distance.log"
exec > >(tee -a "$LOG") 2>&1

date
echo "HOST=$(hostname)"
echo "WORK=$WORK"
echo "THREADS_INDEX=$THREADS_INDEX THREADS_QUANT=$THREADS_QUANT MAX_QUANT_JOBS=$MAX_QUANT_JOBS"

for f in "$GFF" "$FA" "$SVCTX" "$PAV" "$GFFREAD" "$SALMON" "$RSCRIPT"; do
  if [ ! -e "$f" ]; then
    echo "ERROR: missing required input/tool: $f" >&2
    exit 1
  fi
done

MANIFEST="$WORK/config/illumina_pe_samples.tsv"
cat > "$MANIFEST" <<EOF
sample_id	group	r1	r2
Msa_jingA	genome_Msa	$RNADIR/genome_Msa1/illumina-seq/jingA_1.fastq.gz	$RNADIR/genome_Msa1/illumina-seq/jingA_2.fastq.gz
Msa_jingB	genome_Msa	$RNADIR/genome_Msa1/illumina-seq/jingB_1.fastq.gz	$RNADIR/genome_Msa1/illumina-seq/jingB_2.fastq.gz
Msa_jingC	genome_Msa	$RNADIR/genome_Msa1/illumina-seq/jingC_1.fastq.gz	$RNADIR/genome_Msa1/illumina-seq/jingC_2.fastq.gz
Msa_yeA	genome_Msa	$RNADIR/genome_Msa1/illumina-seq/yeA_1.fastq.gz	$RNADIR/genome_Msa1/illumina-seq/yeA_2.fastq.gz
Msa_yeB	genome_Msa	$RNADIR/genome_Msa1/illumina-seq/yeB_1.fastq.gz	$RNADIR/genome_Msa1/illumina-seq/yeB_2.fastq.gz
Msa_yeC	genome_Msa	$RNADIR/genome_Msa1/illumina-seq/yeC_1.fastq.gz	$RNADIR/genome_Msa1/illumina-seq/yeC_2.fastq.gz
genome_410	genome_410	$RNADIR/genome_410/illumina-seq/410_1_clean_R1.fq.gz	$RNADIR/genome_410/illumina-seq/410_1_clean_R2.fq.gz
genome_461	genome_461	$RNADIR/genome_461/illumina-seq/Ms461_2_clean_R1.fq.gz	$RNADIR/genome_461/illumina-seq/Ms461_2_clean_R2.fq.gz
genome_472	genome_472	$RNADIR/genome_472/illumina-seq/472_2_clean_R1.fq.gz	$RNADIR/genome_472/illumina-seq/472_2_clean_R2.fq.gz
genome_474	genome_474	$RNADIR/genome_474a/illumina-seq/474_2_clean_R1.fq.gz	$RNADIR/genome_474a/illumina-seq/474_2_clean_R2.fq.gz
genome_Mar	genome_Mar	$RNADIR/genome_Mar/illumina-seq/Unknown_AJ547-01T0001-02-T00_good_1.fq.gz	$RNADIR/genome_Mar/illumina-seq/Unknown_AJ547-01T0001-02-T00_good_2.fq.gz
Mru_CK_H_1A	genome_Mru	$RNADIR/genome_Mru/illumina-seq/CK-H-1A_1.fq.gz	$RNADIR/genome_Mru/illumina-seq/CK-H-1A_2.fq.gz
Mru_CK_H_2A	genome_Mru	$RNADIR/genome_Mru/illumina-seq/CK-H-2A_1.fq.gz	$RNADIR/genome_Mru/illumina-seq/CK-H-2A_2.fq.gz
Mru_CK_H_3A	genome_Mru	$RNADIR/genome_Mru/illumina-seq/CK-H-3A_1.fq.gz	$RNADIR/genome_Mru/illumina-seq/CK-H-3A_2.fq.gz
Mru_CK_H_4A	genome_Mru	$RNADIR/genome_Mru/illumina-seq/CK-H-4A_1.fq.gz	$RNADIR/genome_Mru/illumina-seq/CK-H-4A_2.fq.gz
genome_R108	genome_R108	$RNADIR/genome_R108/illumina-seq/project_T2T__CRR1954572_r1.fq.gz	$RNADIR/genome_R108/illumina-seq/project_T2T__CRR1954572_r2.fq.gz
ZM4_CRR330267	genome_ZM4	$RNADIR/genome_ZM4/illumina-seq/CRR330267_f1.fq.gz	$RNADIR/genome_ZM4/illumina-seq/CRR330267_r2.fq.gz
ZM4_CRR330268	genome_ZM4	$RNADIR/genome_ZM4/illumina-seq/CRR330268_f1.fq.gz	$RNADIR/genome_ZM4/illumina-seq/CRR330268_r2.fq.gz
ZM4_CRR330269	genome_ZM4	$RNADIR/genome_ZM4/illumina-seq/CRR330269_f1.fq.gz	$RNADIR/genome_ZM4/illumina-seq/CRR330269_r2.fq.gz
EOF

awk 'NR>1 {if(!system("[ -s " $3 " ]") && !system("[ -s " $4 " ]")) ok++; else {print "missing fastq for " $1 > "/dev/stderr"; bad=1}} END{if(bad) exit 1; print "FASTQ pairs OK:", ok}' "$MANIFEST"

TXFA="$WORK/reference/genome_Msa.EVM.transcripts.fa"
TX2GENE="$WORK/reference/tx2gene.tsv"
if [ ! -s "$TXFA" ]; then
  echo "[gffread] extracting transcript fasta"
  "$GFFREAD" -w "$TXFA" -g "$FA" "$GFF"
fi

if [ ! -s "$TX2GENE" ]; then
  echo "[tx2gene] parsing mRNA IDs"
  awk -F'\t' 'BEGIN{OFS="\t"; print "tx_id","gene_id"}
    $0 !~ /^#/ && $3=="mRNA" {
      id=""; parent="";
      n=split($9,a,";");
      for(i=1;i<=n;i++){
        split(a[i],kv,"=");
        if(kv[1]=="ID") id=kv[2];
        if(kv[1]=="Parent") parent=kv[2];
      }
      if(id!="" && parent!="") print id,parent;
    }' "$GFF" > "$TX2GENE"
fi

echo "[reference] transcripts=$(grep -c '^>' "$TXFA") tx2gene=$(($(wc -l < "$TX2GENE")-1))"

INDEX="$WORK/salmon_index/genome_Msa_EVM"
if [ ! -s "$INDEX/ref_k31_fixed.fa" ] && [ ! -s "$INDEX/complete_ref_lens.bin" ]; then
  echo "[salmon index] building index"
  "$SALMON" index -t "$TXFA" -i "$INDEX" -k 31 -p "$THREADS_INDEX"
fi

echo "[salmon quant] start"
declare -a pids=()
run_quant() {
  local sample="$1" r1="$2" r2="$3"
  local out="$WORK/quants/$sample"
  local slog="$WORK/logs/salmon_${sample}.log"
  if [ -s "$out/quant.sf" ]; then
    echo "[salmon quant] skip existing $sample"
    return 0
  fi
  rm -rf "$out"
  mkdir -p "$out"
  echo "[salmon quant] $sample"
  "$SALMON" quant -i "$INDEX" -l A -1 "$r1" -2 "$r2" -p "$THREADS_QUANT" --validateMappings -o "$out" > "$slog" 2>&1
}

while IFS=$'\t' read -r sample group r1 r2; do
  run_quant "$sample" "$r1" "$r2" &
  pids+=("$!")
  while [ "$(jobs -pr | wc -l)" -ge "$MAX_QUANT_JOBS" ]; do
    sleep 15
  done
done < <(tail -n +2 "$MANIFEST")

fail=0
for pid in "${pids[@]:-}"; do
  wait "$pid" || fail=1
done
if [ "$fail" -ne 0 ]; then
  echo "ERROR: one or more salmon quant jobs failed; inspect $WORK/logs/salmon_*.log" >&2
  exit 1
fi

echo "[salmon quant] completed quant.sf files:"
find "$WORK/quants" -maxdepth 2 -name quant.sf -printf '%P\n' | sort

echo "[R] aggregate TPM and plot"
"$RSCRIPT" "$WORK/scripts/plot_mixed_illumina_expression_sv_distance.R" \
  "$WORK" "$MANIFEST" "$TX2GENE" "$SVCTX" "$PAV"

touch "$WORK/summary.done"
date
echo "DONE $WORK"

