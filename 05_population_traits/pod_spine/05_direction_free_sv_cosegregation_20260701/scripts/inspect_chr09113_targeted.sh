#!/usr/bin/env bash
set -euo pipefail

N4="path/to/project/N_4.pod_spiny"
RUN="$N4/05_direction_free_sv_cosegregation_20260701"
ANN="path/to/project/16.T2T_ref_function_anno/output"
OUT="$RUN/summary/Chr09113_function_targeted_20260701.txt"

mkdir -p "$RUN/summary"
: > "$OUT"

add_section() {
  local title="$1"
  shift
  {
    echo
    echo "## $title"
    "$@" || true
  } >> "$OUT"
}

{
  echo "date	$(date '+%F %T')"
  echo "gene	Chr09113"
} >> "$OUT"

add_section "strict_pass_candidate" awk -F $'\t' 'NR==1 || $0 ~ /Chr09113/' "$RUN/results/candidate_strict_pass.tsv"
add_section "Msa_GFF_records" grep -n "Chr09113" "$N4/00_data/4.reference_anno/genome_Msa.gff"
add_section "TAIR_blast_annotation" grep -n "Chr09113" "$ANN/blast_tair/T2T_tait/ref_Msa.T2T_ctg.pep.tab.anno"
add_section "SwissProt_plants_annotation" grep -n "Chr09113" "$ANN/uniprot_sprot_plants/ref_Msa.T2T_ctg.pep_blastp.tab.anno.txt"
add_section "Swiss_annotation_table" grep -n "Chr09113" "$ANN/uniprot_sprot_plants/swiss_annotation.tsv"
add_section "TrEMBL_plants_annotation" grep -n "Chr09113" "$ANN/uniprot_trembl/process_files/ref_Msa.T2T_ctg.pep_blastp.tab.anno.new.txt"
add_section "InterProScan_out" grep -n "Chr09113" "$ANN/interproscan/ref_Msa.T2T_ctg.pep.interproscan.out"
add_section "InterProScan_domains" grep -n "Chr09113" "$ANN/interproscan/ref_Msa.T2T_ctg.pep.domains.txt"
add_section "Pfam_out" grep -n "Chr09113" "$ANN/pfam/ref_Msa.T2T_ctg.pep.pfam.out"
add_section "KOG_blast" grep -n "Chr09113" "$ANN/KOG/ref_Msa.T2T_ctg.pep_blastp.tab"
add_section "candidate_gene_dir_files" find "$N4/00_data/6.candicated_gene" -maxdepth 3 -type f -print
add_section "candidate_gene_dir_hits" grep -RIn "Chr09113\|SPL\|SQUAMOSA\|SBP\|AT1G27370\|AT2G33810\|AT3G15270\|AT5G43270" "$N4/00_data/6.candicated_gene"

cat "$OUT"
