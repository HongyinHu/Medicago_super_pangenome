#!/usr/bin/env bash
set -euo pipefail

base=path/to/project/N_1.coding_gene_anno
out="$base/09_summary_tables"
summary="$out/current_evm_protein_annotation_summary.tsv"
status="$out/current_annotation_stage_status.tsv"
readme="$out/README_current_evm_protein_annotation_stats.txt"

printf 'OUTPUT_FILES\n'
ls -lh "$summary" "$status" "$readme"
printf 'SUMMARY_LINE_COUNT\n'
wc -l "$summary"
printf 'SUMMARY_HEAD\n'
sed -n '1,8p' "$summary"
printf 'SUMMARY_TAIL\n'
tail -n 6 "$summary"
printf 'STAGE_STATUS\n'
cat "$status"
printf 'README\n'
cat "$readme"
printf 'INDEPENDENT_COUNT_CHECK\n'
printf 'assembly\tsummary_genes\tgff_genes\tsummary_proteins\tfasta_proteins\tdone_marker\n'
awk -F '\t' 'NR > 1 {print $1 "\t" $4 "\t" $8 "\t" $32}' "$summary" |
while IFS=$'\t' read -r sample summary_genes summary_proteins done_marker; do
  final="$base/05_EVM_integration/work/$sample/02_final"
  gff_genes=$(awk -F '\t' '!/^#/ && $3 == "gene" {n++} END {print n+0}' "$final/$sample.evm.gff3")
  fasta_proteins=$(grep -c '^>' "$final/$sample.evm.pep.fa")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$summary_genes" "$gff_genes" "$summary_proteins" "$fasta_proteins" "$done_marker"
  test "$summary_genes" -eq "$gff_genes"
  test "$summary_proteins" -eq "$fasta_proteins"
  test "$done_marker" = yes
done
printf 'ALL_INDEPENDENT_COUNTS_MATCH\n'
