#!/usr/bin/env bash
set -euo pipefail

source ~/.bashrc 2>/dev/null || true
conda activate cent_synteny

cd path/to/project/10.centromere_analysis/output_synthen/02_cent_synteny/genome_A17_genome_R108

plot_one() {
  local method=$1
  local pair_file=$2
  local suffix=$3

  python plot_sequence_synteny_CENH3.py \
    --blocks "genome_A17_vs_genome_R108.${method}.blocks.tsv" \
    --chr_pairs "${pair_file}" \
    --query_signal genome_A17.CENH3_vs_Input.q20.10k.log2.bedgraph \
    --target_signal genome_R108.CENH3_vs_Input.q20.10k.log2.bedgraph \
    --query_cen genome_A17.pericentromere_6Mb.bed \
    --target_cen genome_R108.pericentromere_6Mb.bed \
    --out "genome_A17_vs_genome_R108_${suffix}_sequence_${method}_CENH3.min1k.pdf" \
    --min_len 1000 \
    --max_links 0 \
    --query_label A17 \
    --target_label R108 \
    --dpi 300

  python plot_sequence_synteny_CENH3.py \
    --blocks "genome_A17_vs_genome_R108.${method}.blocks.tsv" \
    --chr_pairs "${pair_file}" \
    --query_signal genome_A17.CENH3_vs_Input.q20.10k.log2.bedgraph \
    --target_signal genome_R108.CENH3_vs_Input.q20.10k.log2.bedgraph \
    --query_cen genome_A17.pericentromere_6Mb.bed \
    --target_cen genome_R108.pericentromere_6Mb.bed \
    --out "genome_A17_vs_genome_R108_${suffix}_sequence_${method}_CENH3.min1k.png" \
    --min_len 1000 \
    --max_links 0 \
    --query_label A17 \
    --target_label R108 \
    --dpi 300
}

plot_one rawPAF chr_pairs1-4.tsv chr1-4
plot_one rawPAF chr_pairs5-8.tsv chr5-8
plot_one prenetChain chr_pairs1-4.tsv chr1-4
plot_one prenetChain chr_pairs5-8.tsv chr5-8

ls -lh \
  genome_A17_vs_genome_R108_chr*_sequence_rawPAF_CENH3.min1k.* \
  genome_A17_vs_genome_R108_chr*_sequence_prenetChain_CENH3.min1k.* \
  genome_A17_vs_genome_R108_chr*_sequence_netSyntenic_CENH3.min1k.*
