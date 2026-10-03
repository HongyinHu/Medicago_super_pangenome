#!/usr/bin/env bash
set -euo pipefail

source ~/.bashrc 2>/dev/null || true

BASE=path/to/project
OUT=${BASE}/10.centromere_analysis/output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/TE_first
mkdir -p "${OUT}"
cp path/to/home/calculate_R108_functional_centromere_composition_tesorter_parallel.py "${OUT}/"

conda activate EDTA_env
python "${OUT}/calculate_R108_functional_centromere_composition_tesorter_parallel.py" \
  --centromere_bed "${BASE}/10.centromere_analysis/output_all/04_cenh3_signal/genome_R108.functional_centromere/genome_R108.CENH3.functional_centromere.final.bed" \
  --teanno_gff "${BASE}/8.T2T_TEs_anno_new/output/genome_R108/EDTA/genome_R108.T2T.ctg.final.fa.mod.EDTA.TEanno.gff3" \
  --tesorter_cls "${BASE}/10.centromere_analysis/output_synthen/04_cent_TE_type/TEsorter_R108_test/genome_R108.TElib.rexdb-plant.cls.tsv" \
  --trash_bed "${BASE}/10.centromere_analysis/output_all/03_repeat/genome_R108.TRASH/genome_R108.TRASH_arrays.sorted.bed" \
  --outdir "${OUT}" \
  --label R108 \
  --threads 8 \
  --priority_mode te_first

cp path/to/home/run_06_R108_te_first_patch.sh "${OUT}/run_06_R108_te_first_patch.sh"
ls -lh "${OUT}"
