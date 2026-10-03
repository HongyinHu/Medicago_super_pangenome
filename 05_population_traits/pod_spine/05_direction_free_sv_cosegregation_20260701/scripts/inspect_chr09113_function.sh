#!/usr/bin/env bash
set -euo pipefail

N4="path/to/project/N_4.pod_spiny"
RUN="$N4/05_direction_free_sv_cosegregation_20260701"
ANN="path/to/project/16.T2T_ref_function_anno/output"
OUT="$RUN/summary/Chr09113_function_inspection_20260701.txt"

mkdir -p "$RUN/summary"

{
  echo "date	$(date '+%F %T')"
  echo "gene	Chr09113"
  echo
  echo "## strict_pass_candidate"
  awk -F '\t' 'NR==1 || $0 ~ /Chr09113/' "$RUN/results/candidate_strict_pass.tsv"
  echo
  echo "## Msa_GFF_records"
  grep -n "Chr09113" "$N4/00_data/4.reference_anno/genome_Msa.gff" | head -80 || true
  echo
  echo "## files_containing_Chr09113"
  if command -v rg >/dev/null 2>&1; then
    rg -l "Chr09113" "$ANN" "$N4/00_data" 2>/dev/null | head -80
  else
    grep -RIl "Chr09113" "$ANN" "$N4/00_data" 2>/dev/null | head -80
  fi
  echo
  echo "## annotation_context"
  if command -v rg >/dev/null 2>&1; then
    rg -n "Chr09113|AT[1-5CM]G[0-9]{5}|SPL|SQUAMOSA|SBP|squamosa|transcription factor" "$ANN" "$N4/00_data/6.candicated_gene" 2>/dev/null | head -240
  else
    grep -RInE "Chr09113|AT[1-5CM]G[0-9]{5}|SPL|SQUAMOSA|SBP|squamosa|transcription factor" "$ANN" "$N4/00_data/6.candicated_gene" 2>/dev/null | head -240
  fi
  echo
  echo "## candidate_gene_dir"
  find "$N4/00_data/6.candicated_gene" -maxdepth 3 -type f -print 2>/dev/null | sort | head -120
} > "$OUT"

cat "$OUT"
