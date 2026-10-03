#!/usr/bin/env bash
set -euo pipefail

source ~/.bashrc 2>/dev/null || true

BASE=path/to/project
OUTDIR=${BASE}/10.centromere_analysis/output_synthen/05_functional_centromere_intact_LTR_similarity
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

cp path/to/home/build_functional_centromere_intact_LTR_similarity.py .
cp path/to/home/plot_functional_centromere_intact_LTR_similarity.py .
chmod +x build_functional_centromere_intact_LTR_similarity.py plot_functional_centromere_intact_LTR_similarity.py

cat > functional_centromere_intact_LTRRT.config.tsv <<EOF
genome	label	pass_list	intact_fasta	centromere_bed
genome_A17	A17	${BASE}/8.T2T_TEs_anno_new/output/genome_A17/EDTA/genome_A17.T2T.ctg.final.fa.mod.EDTA.raw/LTR/genome_A17.T2T.ctg.final.fa.mod.pass.list	${BASE}/8.T2T_TEs_anno_new/output/genome_A17/EDTA/genome_A17.T2T.ctg.final.fa.mod.EDTA.intact.fa	${BASE}/10.centromere_analysis/output_all/04_cenh3_signal/genome_A17.functional_centromere/genome_A17.CENH3.functional_centromere.final.bed
genome_R108	R108	${BASE}/8.T2T_TEs_anno_new/output/genome_R108/EDTA/genome_R108.T2T.ctg.final.fa.mod.EDTA.raw/LTR/genome_R108.T2T.ctg.final.fa.mod.pass.list	${BASE}/8.T2T_TEs_anno_new/output/genome_R108/EDTA/genome_R108.T2T.ctg.final.fa.mod.EDTA.intact.fa	${BASE}/10.centromere_analysis/output_all/04_cenh3_signal/genome_R108.functional_centromere/genome_R108.CENH3.functional_centromere.final.bed
genome_Mpo	Mpo	${BASE}/8.T2T_TEs_anno_new/output/genome_Mpo/EDTA/genome_Mpo.T2T.ctg.final.fa.mod.EDTA.raw/LTR/genome_Mpo.T2T.ctg.final.fa.mod.pass.list	${BASE}/8.T2T_TEs_anno_new/output/genome_Mpo/EDTA/genome_Mpo.T2T.ctg.final.fa.mod.EDTA.intact.fa	${BASE}/10.centromere_analysis/output_all/04_cenh3_signal/genome_Mpo.functional_centromere/genome_Mpo.CENH3.functional_centromere.final.bed
genome_Msa	Msa	${BASE}/8.T2T_TEs_anno_new/output/genome_Msa/EDTA/genome_Msa.T2T.ctg.final.fa.mod.EDTA.raw/LTR/genome_Msa.T2T.ctg.final.fa.mod.pass.list	${BASE}/8.T2T_TEs_anno_new/output/genome_Msa/EDTA/genome_Msa.T2T.ctg.final.fa.mod.EDTA.intact.fa	${BASE}/10.centromere_analysis/output_all/04_cenh3_signal/genome_Msa.functional_centromere/genome_Msa.CENH3.functional_centromere.final.bed
genome_474	474	${BASE}/8.T2T_TEs_anno_new/output/genome_474/EDTA/genome_474.near_T2T.ctg.final.fa.mod.EDTA.raw/LTR/genome_474.near_T2T.ctg.final.fa.mod.pass.list	${BASE}/8.T2T_TEs_anno_new/output/genome_474/EDTA/genome_474.near_T2T.ctg.final.fa.mod.EDTA.intact.fa	${BASE}/10.centromere_analysis/output_all/04_cenh3_signal/genome_474.functional_centromere_topN/genome_474.CENH3.primary_functional_centromere.bed
EOF

python build_functional_centromere_intact_LTR_similarity.py \
  --config functional_centromere_intact_LTRRT.config.tsv \
  --outdir "${OUTDIR}"

if [[ ! -s functional_centromere.intact_LTRRT.fa ]]; then
  echo "No functional-centromere intact LTR-RT sequences were extracted" >&2
  exit 1
fi

conda activate EDTA_env
if [[ ! -s functional_centromere.intact_LTRRT.rexdb-plant.cls.tsv ]]; then
  TEsorter functional_centromere.intact_LTRRT.fa \
    -db rexdb-plant \
    -pre functional_centromere.intact_LTRRT.rexdb-plant \
    -p 16
fi

conda activate cent_synteny
if [[ ! -s functional_centromere.intact_LTRRT.self.paf ]]; then
  minimap2 -x asm20 -c --secondary=yes -N 1000 -t 32 \
    functional_centromere.intact_LTRRT.fa \
    functional_centromere.intact_LTRRT.fa \
    > functional_centromere.intact_LTRRT.self.paf
fi

conda activate EDTA_env
python plot_functional_centromere_intact_LTR_similarity.py \
  --metadata functional_centromere.intact_LTRRT.metadata.tsv \
  --fasta functional_centromere.intact_LTRRT.fa \
  --paf functional_centromere.intact_LTRRT.self.paf \
  --tesorter_cls functional_centromere.intact_LTRRT.rexdb-plant.cls.tsv \
  --out_prefix functional_centromere.intact_LTRRT.similarity_heatmap \
  --min_cov 0.5

cat > README_05_functional_centromere_intact_LTR_similarity.md <<EOF
# Functional centromere intact LTR-RT similarity heatmap

Created at: $(date '+%F %T')

This analysis draws a species-level intact LTR-RT similarity heatmap for the five genomes,
using only intact LTR-RTs overlapping CENH3-defined functional centromere regions.

## Scope

- Region: functional centromere only.
- Chromosome filter: only \`^Chr[0-9]+\`; scaffold records excluded.
- LTR source: EDTA intact LTR-RT \`*.EDTA.intact.fa\` and \`*.pass.list\`.
- Insertion time: \`Insertion_Time\` from EDTA pass.list, converted to Mya.
- Clade annotation: TEsorter on extracted functional-centromere intact LTR-RTs with \`-db rexdb-plant\`.
- Similarity: minimap2 all-vs-all PAF, using max identity per pair with alignment coverage >= 0.5 of the shorter sequence.

## Inputs

Input paths are recorded in:

\`\`\`text
functional_centromere_intact_LTRRT.config.tsv
\`\`\`

## Main commands

\`\`\`bash
python build_functional_centromere_intact_LTR_similarity.py --config functional_centromere_intact_LTRRT.config.tsv --outdir ${OUTDIR}
conda activate EDTA_env
TEsorter functional_centromere.intact_LTRRT.fa -db rexdb-plant -pre functional_centromere.intact_LTRRT.rexdb-plant -p 16
conda activate cent_synteny
minimap2 -x asm20 -c --secondary=yes -N 1000 -t 32 functional_centromere.intact_LTRRT.fa functional_centromere.intact_LTRRT.fa > functional_centromere.intact_LTRRT.self.paf
conda activate EDTA_env
python plot_functional_centromere_intact_LTR_similarity.py --metadata functional_centromere.intact_LTRRT.metadata.tsv --fasta functional_centromere.intact_LTRRT.fa --paf functional_centromere.intact_LTRRT.self.paf --tesorter_cls functional_centromere.intact_LTRRT.rexdb-plant.cls.tsv --out_prefix functional_centromere.intact_LTRRT.similarity_heatmap --min_cov 0.5
\`\`\`

## Outputs

\`\`\`text
functional_centromere.intact_LTRRT.fa
functional_centromere.intact_LTRRT.metadata.tsv
functional_centromere.intact_LTRRT.counts.tsv
functional_centromere.intact_LTRRT.rexdb-plant.cls.tsv
functional_centromere.intact_LTRRT.self.paf
functional_centromere.intact_LTRRT.similarity_heatmap.similarity_matrix.tsv
functional_centromere.intact_LTRRT.similarity_heatmap.ordered_metadata.tsv
functional_centromere.intact_LTRRT.similarity_heatmap.pdf
functional_centromere.intact_LTRRT.similarity_heatmap.png
\`\`\`
EOF

ls -lh \
  functional_centromere.intact_LTRRT.counts.tsv \
  functional_centromere.intact_LTRRT.metadata.tsv \
  functional_centromere.intact_LTRRT.rexdb-plant.cls.tsv \
  functional_centromere.intact_LTRRT.self.paf \
  functional_centromere.intact_LTRRT.similarity_heatmap.pdf \
  functional_centromere.intact_LTRRT.similarity_heatmap.png \
  README_05_functional_centromere_intact_LTR_similarity.md
