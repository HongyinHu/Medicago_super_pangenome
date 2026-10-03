#!/usr/bin/env bash
set -euo pipefail

source ~/.bashrc 2>/dev/null || true

BASE=path/to/project
OUTDIR=${BASE}/10.centromere_analysis/output_synthen/05_functional_centromere_intact_LTR_similarity
cd "${OUTDIR}"

cp path/to/home/plot_functional_centromere_intact_LTR_similarity.py .
chmod +x plot_functional_centromere_intact_LTR_similarity.py

if [[ ! -s functional_centromere.intact_LTRRT.fa ]]; then
  echo "Missing functional_centromere.intact_LTRRT.fa; run run_05_functional_centromere_intact_LTR_similarity.sh first" >&2
  exit 1
fi

if [[ ! -s functional_centromere.intact_LTRRT.rexdb-plant.cls.tsv ]]; then
  conda activate EDTA_env
  TEsorter functional_centromere.intact_LTRRT.fa \
    -db rexdb-plant \
    -pre functional_centromere.intact_LTRRT.rexdb-plant \
    -p 16
fi

conda activate EDTA_env
if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v blastn >/dev/null 2>&1; then
  conda activate cent_synteny
fi
if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v blastn >/dev/null 2>&1; then
  echo "makeblastdb/blastn are not available in EDTA_env or cent_synteny" >&2
  exit 1
fi

if [[ ! -s functional_centromere.intact_LTRRT.blastdb.nsq ]]; then
  makeblastdb \
    -in functional_centromere.intact_LTRRT.fa \
    -dbtype nucl \
    -out functional_centromere.intact_LTRRT.blastdb
fi

if [[ ! -s functional_centromere.intact_LTRRT.self.blastn.tsv ]]; then
  blastn \
    -query functional_centromere.intact_LTRRT.fa \
    -db functional_centromere.intact_LTRRT.blastdb \
    -task blastn \
    -dust no \
    -evalue 1e-5 \
    -num_threads 32 \
    -max_target_seqs 100000 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -out functional_centromere.intact_LTRRT.self.blastn.tsv
fi

conda activate EDTA_env
python plot_functional_centromere_intact_LTR_similarity.py \
  --metadata functional_centromere.intact_LTRRT.metadata.tsv \
  --fasta functional_centromere.intact_LTRRT.fa \
  --blastn functional_centromere.intact_LTRRT.self.blastn.tsv \
  --tesorter_cls functional_centromere.intact_LTRRT.rexdb-plant.cls.tsv \
  --out_prefix functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident \
  --blast_min_cov 0 \
  --blast_min_aln_len 0 \
  --collapse_copia \
  --gypsy_top_n 5

python plot_functional_centromere_intact_LTR_similarity.py \
  --metadata functional_centromere.intact_LTRRT.metadata.tsv \
  --fasta functional_centromere.intact_LTRRT.fa \
  --blastn functional_centromere.intact_LTRRT.self.blastn.tsv \
  --tesorter_cls functional_centromere.intact_LTRRT.rexdb-plant.cls.tsv \
  --out_prefix functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.minCov20 \
  --blast_min_cov 0.2 \
  --blast_min_aln_len 100 \
  --collapse_copia \
  --gypsy_top_n 5

cat >> README_05_functional_centromere_intact_LTR_similarity.md <<EOF

## BLASTN/TBtools-style Similarity Trial

Added at: $(date '+%F %T')

This trial follows the paper-style idea more closely than the earlier minimap2 whole-element PAF figure:
all functional-centromere intact LTR-RT nucleotide sequences were compared all-vs-all by BLASTN, and each pair was scored by the maximum local alignment identity.

This is intentionally looser than the minimap2 \`--min_cov 0.5\` matrix. It tests whether the weak signal in the first heatmap was caused by requiring large whole-element alignments rather than local LTR-RT similarity.

Commands:

\`\`\`bash
bash run_05_functional_centromere_intact_LTR_similarity_blastn.sh

makeblastdb -in functional_centromere.intact_LTRRT.fa -dbtype nucl -out functional_centromere.intact_LTRRT.blastdb

blastn -query functional_centromere.intact_LTRRT.fa \\
  -db functional_centromere.intact_LTRRT.blastdb \\
  -task blastn -dust no -evalue 1e-5 -num_threads 32 -max_target_seqs 100000 \\
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \\
  -out functional_centromere.intact_LTRRT.self.blastn.tsv
\`\`\`

Outputs:

\`\`\`text
functional_centromere.intact_LTRRT.self.blastn.tsv
functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.pdf
functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.png
functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.similarity_matrix.tsv
functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.minCov20.pdf
functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.minCov20.png
functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.minCov20.similarity_matrix.tsv
\`\`\`
EOF

ls -lh \
  functional_centromere.intact_LTRRT.self.blastn.tsv \
  functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.pdf \
  functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.png \
  functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.minCov20.pdf \
  functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.minCov20.png \
  README_05_functional_centromere_intact_LTR_similarity.md
