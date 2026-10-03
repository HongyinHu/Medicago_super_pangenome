#!/usr/bin/env bash
set -euo pipefail

OUT=path/to/project/N_4.pod_spiny/29_RcMYB106_promoter_SV_R108_M22_20260714
ROOT=path/to/project/N_4.pod_spiny

python3 "$OUT/scripts/analyze_rcmyb106_promoter.py" \
  --outdir "$OUT" \
  --r108-genome "$ROOT/00_data/1.reference_genome/genome_R108.fa" \
  --r108-gff "$ROOT/00_data/4.reference_anno/genome_R108.gff" \
  --r108-pep path/to/project/9.T2T_gene_anno_new/output/4.evm_combind/genome_R108/test2_finally/R108_genome.anno.pep \
  --m22-genome "$ROOT/00_data/1.reference_genome/genome_M22.fa" \
  --m22-gff "$ROOT/00_data/4.reference_anno/genome_M22.gff" \
  --m22-pep path/to/projects_all/pan_genome/X_genome_M22_0/9.genome_annotation_gene_predict_new/output/5.evidencemodeler_combind/finally.genome.anno.pep \
  --blastp path/to/project/biosofeware/ncbi-blast-2.13.0+/bin/blastp \
  --minimap2 path/to/project/biosofeware/NextPolish/bin/minimap2
