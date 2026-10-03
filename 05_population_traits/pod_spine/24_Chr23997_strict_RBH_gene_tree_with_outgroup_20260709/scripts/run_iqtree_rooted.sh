#!/usr/bin/env bash
set -euo pipefail
OUT=path/to/project/N_4.pod_spiny/24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709
BIN=path/to/home/anaconda3/envs/biosofeware/bin
OG=outgroup_Msa_Chr23998
mkdir -p "$OUT/logs" "$OUT/trees" "$OUT/figures" "$OUT/summary"
: > "$OUT/logs/pad_alignments.log"
python3 "$OUT/scripts/pad_fasta_alignment.py" \
  "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.cds.codon.aln.fa" \
  "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.cds.codon.aln.fa" >> "$OUT/logs/pad_alignments.log" 2>&1
python3 "$OUT/scripts/make_codon12_alignment.py" \
  "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.cds.codon.aln.fa" \
  "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.codon12.aln.fa" > "$OUT/logs/make_codon12.log" 2>&1
python3 "$OUT/scripts/pad_fasta_alignment.py" \
  "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.codon12.aln.fa" \
  "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.codon12.aln.fa" >> "$OUT/logs/pad_alignments.log" 2>&1
if [ ! -s "$OUT/trees/Chr23997_pep_ML_rooted.treefile" ]; then
  "$BIN/iqtree2" -s "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.pep.aln.fa" -st AA -m MFP -bb 1000 -alrt 1000 -nt 4 -o "$OG" -pre "$OUT/trees/Chr23997_pep_ML_rooted" > "$OUT/logs/iqtree_pep.log" 2>&1
fi
if [ ! -s "$OUT/trees/Chr23997_cds_ML_rooted.treefile" ]; then
  "$BIN/iqtree2" -s "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.cds.codon.aln.fa" -st DNA -m MFP -bb 1000 -alrt 1000 -nt 4 -o "$OG" -pre "$OUT/trees/Chr23997_cds_ML_rooted" -redo > "$OUT/logs/iqtree_cds.log" 2>&1
fi
if [ ! -s "$OUT/trees/Chr23997_codon12_ML_rooted.treefile" ]; then
  "$BIN/iqtree2" -s "$OUT/alignments/Chr23997_strict_RBH_plus_outgroup.codon12.aln.fa" -st DNA -m MFP -bb 1000 -alrt 1000 -nt 4 -o "$OG" -pre "$OUT/trees/Chr23997_codon12_ML_rooted" -redo > "$OUT/logs/iqtree_codon12.log" 2>&1
fi
Rscript "$OUT/scripts/plot_newick_baseR.R" "$OUT" > "$OUT/logs/plot_newick_baseR.log" 2>&1
