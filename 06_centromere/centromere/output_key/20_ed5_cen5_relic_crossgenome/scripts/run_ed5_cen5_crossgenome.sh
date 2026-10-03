#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/10.centromere_analysis/output_key
RUN=${ROOT}/20_ed5_cen5_relic_crossgenome
PREV_RAW=${ROOT}/15_cen5_relic_raw_rerun
PREV_STRONG=${ROOT}/16_cen5_strong_validation
CONDA=path/to/home/anaconda3/bin/conda
THREADS=${THREADS:-48}

mkdir -p "${RUN}"/{fasta,index,paf,results,logs,status,scripts}

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${RUN}/logs/run.log"
}

require_file() {
  if [[ ! -s "$1" ]]; then
    printf 'Required file missing or empty: %s\n' "$1" >&2
    exit 1
  fi
}

extract_region() {
  local fasta=$1
  local chrom=$2
  local start0=$3
  local end=$4
  local name=$5
  local start1=$((start0 + 1))
  "${CONDA}" run -n biosofeware samtools faidx \
    "${fasta}" "${chrom}:${start1}-${end}" \
    | awk -v name="${name}" '/^>/{print ">" name; next} {print}'
}

index_genome() {
  local species=$1
  local fasta="${ROOT}path/to/data/${species}.fa"
  local index="${RUN}/index/${species}.asm20.mmi"
  require_file "${fasta}"
  if [[ ! -s "${fasta}.fai" ]]; then
    "${CONDA}" run -n biosofeware samtools faidx "${fasta}"
  fi
  if [[ ! -s "${index}" ]]; then
    log "Build minimap2 index: ${species}"
    "${CONDA}" run -n minimap2_env minimap2 -I 16G -d "${index}" "${fasta}" \
      2> "${RUN}/logs/${species}.index.stderr"
  fi
}

align_query() {
  local target_species=$1
  local query_fasta=$2
  local output_paf=$3
  local max_secondary=$4
  if [[ ! -s "${output_paf}" ]]; then
    log "Align $(basename "${query_fasta}") to ${target_species}"
    "${CONDA}" run -n minimap2_env minimap2 \
      -x asm20 -c --cs=short --secondary=yes \
      -N "${max_secondary}" -p 0.01 -t "${THREADS}" \
      "${RUN}/index/${target_species}.asm20.mmi" "${query_fasta}" \
      > "${output_paf}" \
      2> "${output_paf%.paf}.stderr"
  fi
}

log "Validate source data"
require_file "${PREV_STRONG}/results/x8_raw_q20_called_CENH3_domains.tsv"
require_file "${PREV_RAW}/results/Mpo_CEN5_relic_focus_candidates.tsv"
require_file "${PREV_RAW}/results/Mpo_raw_q20_active_CENH3_domains.bed"
require_file "${RUN}/scripts/summarize_ed5_cen5_crossgenome.py"

for species in genome_Mpo genome_R108 genome_A17 genome_Msa genome_474; do
  index_genome "${species}"
done

MANIFEST="${RUN}/results/input_intervals.tsv"
printf 'feature_id\tspecies\tdisplay_species\tfeature_type\tchromosome\tstart_bp\tend_bp\tsource\n' > "${MANIFEST}"

declare -A DISPLAY
DISPLAY[genome_R108]=Mtru_R108
DISPLAY[genome_A17]=Mtru_A17
DISPLAY[genome_Msa]=Msat_cae
DISPLAY[genome_474]=Mcar

X8_CORES="${RUN}/fasta/x8_CEN5_cores.fa"
: > "${X8_CORES}"
for species in genome_R108 genome_A17 genome_Msa genome_474; do
  read -r chrom start end < <(
    awk -v sp="${species}" 'BEGIN{FS=OFS="\t"} NR>1 && $1==sp && $3=="Chr5"{print $3,$4,$5; exit}' \
      "${PREV_STRONG}/results/x8_raw_q20_called_CENH3_domains.tsv"
  )
  if [[ -z "${chrom:-}" || -z "${start:-}" || -z "${end:-}" ]]; then
    printf 'Unable to find active CEN5 for %s\n' "${species}" >&2
    exit 1
  fi
  feature="${species}_CEN5_core"
  printf '%s\t%s\t%s\tactive_CEN5\t%s\t%s\t%s\traw_q20_CENH3_domain\n' \
    "${feature}" "${species}" "${DISPLAY[$species]}" "${chrom}" "${start}" "${end}" \
    >> "${MANIFEST}"
  extract_region "${ROOT}path/to/data/${species}.fa" \
    "${chrom}" "${start}" "${end}" "${feature}" >> "${X8_CORES}"
done

MPO_RELICS="${RUN}/fasta/Mpo_CEN5_relics.fa"
MPO_FLANKS="${RUN}/fasta/Mpo_CEN5_relics_pm2Mb.fa"
: > "${MPO_RELICS}"
: > "${MPO_FLANKS}"
for spec in "Mpo_Chr3_relic Chr3 20591033" "Mpo_Chr5_relic Chr5 25293760"; do
  read -r feature chrom expected_start <<< "${spec}"
  read -r start end < <(
    awk -v chr="${chrom}" -v s="${expected_start}" \
      'BEGIN{FS=OFS="\t"} NR>1 && $1==chr && $2==s{print $2,$3; exit}' \
      "${PREV_RAW}/results/Mpo_CEN5_relic_focus_candidates.tsv"
  )
  if [[ -z "${start:-}" || -z "${end:-}" ]]; then
    printf 'Unable to find Mpol relic %s %s\n' "${chrom}" "${expected_start}" >&2
    exit 1
  fi
  printf '%s\tgenome_Mpo\tMpol\tcandidate_CEN5_relic\t%s\t%s\t%s\traw_R108_CEN5_projection_candidate\n' \
    "${feature}" "${chrom}" "${start}" "${end}" >> "${MANIFEST}"
  extract_region "${ROOT}path/to/data" \
    "${chrom}" "${start}" "${end}" "${feature}" >> "${MPO_RELICS}"

  chr_len=$(awk -v chr="${chrom}" '$1==chr{print $2; exit}' \
    "${ROOT}path/to/data")
  flank_start=$((start - 2000000))
  [[ "${flank_start}" -lt 0 ]] && flank_start=0
  flank_end=$((end + 2000000))
  [[ "${flank_end}" -gt "${chr_len}" ]] && flank_end="${chr_len}"
  flank_feature="${feature}_pm2Mb"
  printf '%s\tgenome_Mpo\tMpol\trelic_context_pm2Mb\t%s\t%s\t%s\tderived_from_candidate_relic\n' \
    "${flank_feature}" "${chrom}" "${flank_start}" "${flank_end}" >> "${MANIFEST}"
  extract_region "${ROOT}path/to/data" \
    "${chrom}" "${flank_start}" "${flank_end}" "${flank_feature}" >> "${MPO_FLANKS}"
done

while IFS=$'\t' read -r chrom start end name; do
  [[ -z "${chrom}" ]] && continue
  printf '%s\tgenome_Mpo\tMpol\tactive_CENH3\t%s\t%s\t%s\traw_q20_CENH3_domain\n' \
    "${name}" "${chrom}" "${start}" "${end}" >> "${MANIFEST}"
done < "${PREV_RAW}/results/Mpo_raw_q20_active_CENH3_domains.bed"

log "Forward projection: each x=8 active CEN5 core to Mpol"
align_query genome_Mpo "${X8_CORES}" \
  "${RUN}/paf/x8_CEN5_cores_to_Mpo.asm20.paf" 5000

log "Reverse projection: both Mpol relic candidates to each x=8 genome"
for species in genome_R108 genome_A17 genome_Msa genome_474; do
  align_query "${species}" "${MPO_RELICS}" \
    "${RUN}/paf/Mpo_relics_to_${species}.asm20.paf" 5000
  align_query "${species}" "${MPO_FLANKS}" \
    "${RUN}/paf/Mpo_relic_context_pm2Mb_to_${species}.asm20.paf" 1000
done

log "Summarize bidirectional alignments, non-overlapping coverage and threshold sensitivity"
"${CONDA}" run -n biosofeware python \
  "${RUN}/scripts/summarize_ed5_cen5_crossgenome.py" \
  --run "${RUN}"

{
  printf 'completed_at\t%s\n' "$(date -Iseconds)"
  printf 'node\t%s\n' "$(hostname)"
  printf 'threads\t%s\n' "${THREADS}"
  printf 'minimap2_version\t'
  "${CONDA}" run -n minimap2_env minimap2 --version
  printf 'samtools_version\t'
  "${CONDA}" run -n biosofeware samtools --version | sed -n '1p'
} > "${RUN}/status/RUN_METADATA.tsv"

touch "${RUN}/status/DONE"
log "ED5 CEN5 cross-genome refinement finished"
