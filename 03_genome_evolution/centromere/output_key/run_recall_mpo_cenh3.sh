#!/usr/bin/env bash
set -euo pipefail

PROJECT="path/to/project/10.centromere_analysis"
RUN_ROOT="${PROJECT}/output_key"
OUT="${RUN_ROOT}/12_recall_Mpo_CENH3"
SIG="${PROJECT}/output_all/04_cenh3_signal"
SAMTOOLS="path/to/home/anaconda3/bin/samtools"
CONDA="path/to/home/anaconda3/bin/conda"
PYTHON="path/to/home/anaconda3/bin/python"

mkdir -p "${OUT}/logs" "${OUT}/bedgraph" "${OUT}/results" "${OUT}/scripts"

echo "[INFO] project=${PROJECT}"
echo "[INFO] out=${OUT}"
date

"${SAMTOOLS}" idxstats "${SIG}/genome_Mpo.CENH3.raw.bam" \
  | awk '$1!="*" {print $1"\t"$2}' \
  > "${OUT}/results/genome_Mpo.genome.sizes"

echo "[INFO] BAM idxstats"
: > "${OUT}/results/Mpo_CENH3_Input_bam_idxstats_summary.tsv"
for mode in raw relaxed q20; do
  "${SAMTOOLS}" idxstats "${SIG}/genome_Mpo.CENH3.${mode}.bam" \
    | awk -v mode="${mode}" '$1!="*" {mapped+=$3; unmapped+=$4} END{print mode"\tCENH3\t"mapped"\t"unmapped}' \
    >> "${OUT}/results/Mpo_CENH3_Input_bam_idxstats_summary.tsv"
  "${SAMTOOLS}" idxstats "${SIG}/genome_Mpo.Input.${mode}.bam" \
    | awk -v mode="${mode}" '$1!="*" {mapped+=$3; unmapped+=$4} END{print mode"\tInput\t"mapped"\t"unmapped}' \
    >> "${OUT}/results/Mpo_CENH3_Input_bam_idxstats_summary.tsv"
done

echo "[INFO] bamCompare CENH3/Input"
for mode in raw relaxed q20; do
  out_bg="${OUT}/bedgraph/genome_Mpo.CENH3_vs_Input.${mode}.10k.smooth50k.CPM.log2.bedGraph"
  if [[ ! -s "${out_bg}" ]]; then
    echo "[INFO] running bamCompare for ${mode}"
    "${CONDA}" run -n CENH3_env bamCompare \
      -b1 "${SIG}/genome_Mpo.CENH3.${mode}.bam" \
      -b2 "${SIG}/genome_Mpo.Input.${mode}.bam" \
      --operation log2 \
      --binSize 10000 \
      --smoothLength 50000 \
      --normalizeUsing CPM \
      --scaleFactorsMethod None \
      --pseudocount 1 1 \
      --numberOfProcessors 8 \
      --outFileFormat bedgraph \
      -o "${out_bg}" \
      > "${OUT}/logs/bamCompare.${mode}.log" \
      2> "${OUT}/logs/bamCompare.${mode}.err"
  else
    echo "[INFO] using existing ${out_bg}"
  fi
done

echo "[INFO] calling broad CENH3 domains"
"${PYTHON}" "${OUT}/scripts/call_mpo_cenh3_domains.py" \
  --outdir "${OUT}/results" \
  --genome-sizes "${OUT}/results/genome_Mpo.genome.sizes" \
  --old-core-bed "${PROJECT}/output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.cen_core.bed" \
  --threshold 0.5 \
  --cluster-gap 50000 \
  --domain-gap 1200000 \
  --min-seed-score 50000 \
  raw="${OUT}/bedgraph/genome_Mpo.CENH3_vs_Input.raw.10k.smooth50k.CPM.log2.bedGraph" \
  relaxed="${OUT}/bedgraph/genome_Mpo.CENH3_vs_Input.relaxed.10k.smooth50k.CPM.log2.bedGraph" \
  q20="${OUT}/bedgraph/genome_Mpo.CENH3_vs_Input.q20.10k.smooth50k.CPM.log2.bedGraph" \
  > "${OUT}/logs/call_mpo_cenh3_domains.log" \
  2> "${OUT}/logs/call_mpo_cenh3_domains.err"

echo "[INFO] consensus"
cat "${OUT}/results/Mpo_CENH3_recall_consensus.tsv"
date
