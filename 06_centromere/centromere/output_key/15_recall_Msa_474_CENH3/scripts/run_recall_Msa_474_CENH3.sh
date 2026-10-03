#!/usr/bin/env bash
set -euo pipefail

PROJECT="path/to/project/10.centromere_analysis"
RUN_ROOT="${PROJECT}/output_key/15_recall_Msa_474_CENH3"
SIG="${PROJECT}/output_all/04_cenh3_signal"
SAMTOOLS="path/to/home/anaconda3/bin/samtools"
CONDA="path/to/home/anaconda3/bin/conda"
PYTHON="path/to/home/anaconda3/bin/python"

mkdir -p "${RUN_ROOT}/scripts" "${RUN_ROOT}/logs"

echo "[INFO] project=${PROJECT}"
echo "[INFO] run_root=${RUN_ROOT}"
date

old_core_for_species() {
  case "$1" in
    genome_Msa)
      echo "${PROJECT}/output_synthen/01_JCVI/R108_Msa_synteny_lz10_run/genome_Msa.cen_core.bed"
      ;;
    genome_474)
      echo "${PROJECT}/output_synthen/01_JCVI/genome_R108_genome_474_synteny_lz10_run/genome_474.cen_core.bed"
      ;;
    *)
      echo ""
      ;;
  esac
}

run_one_species() {
  local species="$1"
  local out="${RUN_ROOT}/${species}"
  local old_core
  old_core="$(old_core_for_species "${species}")"

  mkdir -p "${out}/logs" "${out}/bedgraph" "${out}/results" "${out}/figures"
  echo "[INFO] ===== ${species} ====="
  date

  "${SAMTOOLS}" idxstats "${SIG}/${species}.CENH3.raw.bam" \
    | awk '$1!="*" {print $1"\t"$2}' \
    > "${out}/results/${species}.genome.sizes"

  {
    echo -e "mode\tsample\tmapped\tunmapped"
    for mode in raw relaxed q20; do
      "${SAMTOOLS}" idxstats "${SIG}/${species}.CENH3.${mode}.bam" \
        | awk -v mode="${mode}" '$1!="*" {mapped+=$3; unmapped+=$4} END{print mode"\tCENH3\t"mapped"\t"unmapped}'
      "${SAMTOOLS}" idxstats "${SIG}/${species}.Input.${mode}.bam" \
        | awk -v mode="${mode}" '$1!="*" {mapped+=$3; unmapped+=$4} END{print mode"\tInput\t"mapped"\t"unmapped}'
    done
  } > "${out}/results/${species}_CENH3_Input_bam_idxstats_summary.tsv"

  for mode in raw relaxed q20; do
    local out_bg="${out}/bedgraph/${species}.CENH3_vs_Input.${mode}.10k.smooth50k.CPM.log2.bedGraph"
    if [[ ! -s "${out_bg}" ]]; then
      echo "[INFO] ${species} bamCompare ${mode}"
      "${CONDA}" run -n CENH3_env bamCompare \
        -b1 "${SIG}/${species}.CENH3.${mode}.bam" \
        -b2 "${SIG}/${species}.Input.${mode}.bam" \
        --operation log2 \
        --binSize 10000 \
        --smoothLength 50000 \
        --normalizeUsing CPM \
        --scaleFactorsMethod None \
        --pseudocount 1 1 \
        --numberOfProcessors 8 \
        --outFileFormat bedgraph \
        -o "${out_bg}" \
        > "${out}/logs/bamCompare.${mode}.log" \
        2> "${out}/logs/bamCompare.${mode}.err"
    else
      echo "[INFO] using existing ${out_bg}"
    fi
  done

  echo "[INFO] ${species} calling domains"
  "${PYTHON}" "${RUN_ROOT}/scripts/call_cenh3_domains_generic.py" \
    --species "${species}" \
    --outdir "${out}/results" \
    --genome-sizes "${out}/results/${species}.genome.sizes" \
    --old-core-bed "${old_core}" \
    --threshold 0.5 \
    --cluster-gap 50000 \
    --domain-gap 1200000 \
    --min-seed-score 50000 \
    raw="${out}/bedgraph/${species}.CENH3_vs_Input.raw.10k.smooth50k.CPM.log2.bedGraph" \
    relaxed="${out}/bedgraph/${species}.CENH3_vs_Input.relaxed.10k.smooth50k.CPM.log2.bedGraph" \
    q20="${out}/bedgraph/${species}.CENH3_vs_Input.q20.10k.smooth50k.CPM.log2.bedGraph" \
    > "${out}/logs/call_cenh3_domains.log" \
    2> "${out}/logs/call_cenh3_domains.err"

  echo "[INFO] ${species} plotting"
  "${PYTHON}" "${RUN_ROOT}/scripts/plot_cenh3_recall_signal_generic.py" \
    --species "${species}" \
    --root "${RUN_ROOT}" \
    > "${out}/logs/plot_cenh3_recall_signal.log" \
    2> "${out}/logs/plot_cenh3_recall_signal.err"

  cat "${out}/results/${species}_CENH3_recall_consensus.tsv"
}

run_one_species genome_Msa
run_one_species genome_474

echo "[INFO] all done"
date
