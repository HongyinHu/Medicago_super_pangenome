#!/usr/bin/env bash
set -eo pipefail

ROOT="path/to/project/10.centromere_analysis"
RUN="${ROOT}/output_key/15_cen5_relic_raw_rerun"
SCRIPT="${RUN}/scripts/analyze_cen5_raw_rerun.py"
THREADS="${THREADS:-32}"
PYTHON_BIN="${PYTHON_BIN:-path/to/home/anaconda3/bin/python}"

R108_FA="${ROOT}/output_key/data/00_genome/genome_R108.fa"
MPO_FA="${ROOT}/output_key/data/00_genome/genome_Mpo.fa"

R108_CENH3_R1="path/to/project/0.raw_data/data_R108/data_CENH3_public/CRR1954573_r1.fq.gz"
R108_CENH3_R2="path/to/project/0.raw_data/data_R108/data_CENH3_public/CRR1954573_r2.fq.gz"
R108_INPUT_R1="path/to/project/0.raw_data/data_R108/data_CENH3_public/R108_chip_input/CRR1954574_r1.fq.gz"
R108_INPUT_R2="path/to/project/0.raw_data/data_R108/data_CENH3_public/R108_chip_input/CRR1954574_r2.fq.gz"

MPO_CENH3_R1="${ROOT}/output_key/data/02_chipseq/genome_Mpo/02_chipseq_clean/genome_Mpo_CENH3.clean.R1.fq.gz"
MPO_CENH3_R2="${ROOT}/output_key/data/02_chipseq/genome_Mpo/02_chipseq_clean/genome_Mpo_CENH3.clean.R2.fq.gz"
MPO_INPUT_R1="${ROOT}/output_key/data/02_chipseq/genome_Mpo/02_chipseq_clean/genome_Mpo_Input.clean.R1.fq.gz"
MPO_INPUT_R2="${ROOT}/output_key/data/02_chipseq/genome_Mpo/02_chipseq_clean/genome_Mpo_Input.clean.R2.fq.gz"

mkdir -p "${RUN}"/{logs,index,mapping,windows,coverage,chrom_sizes,results/signal,fasta,paf}

log() {
  echo "[$(date '+%F %T')] $*" | tee -a "${RUN}/logs/run.log"
}

require_file() {
  if [[ ! -s "$1" ]]; then
    echo "Required file missing or empty: $1" >&2
    exit 1
  fi
}

bed_to_fasta() {
  local bed="$1"
  local fa="$2"
  local out="$3"
  : > "${out}"
  while IFS=$'\t' read -r chrom start end name rest; do
    [[ -z "${chrom}" ]] && continue
    local one_based=$((start + 1))
    samtools faidx "${fa}" "${chrom}:${one_based}-${end}" >> "${out}"
  done < "${bed}"
}

map_chip_bed() {
  local species="$1"
  local sample="$2"
  local ref="$3"
  local index="$4"
  local r1="$5"
  local r2="$6"
  local out="${RUN}/mapping/${species}.${sample}.q20.bed"
  if [[ -s "${out}" ]]; then
    log "Reuse existing ${out}"
    return
  fi
  log "Map ${species} ${sample} from fastq with chromap MAPQ20"
  conda activate chromap
  chromap --preset chip \
    -t "${THREADS}" \
    -x "${index}" \
    -r "${ref}" \
    -1 "${r1}" \
    -2 "${r2}" \
    -q 20 \
    --remove-pcr-duplicates \
    --BED \
    -o "${out}.tmp" \
    --summary "${RUN}/logs/${species}.${sample}.q20.chromap.summary.txt" \
    2> "${RUN}/logs/${species}.${sample}.q20.chromap.stderr"
  conda deactivate
  conda activate biosofeware
  sort -k1,1 -k2,2n "${out}.tmp" > "${out}"
  rm -f "${out}.tmp"
  conda deactivate
}

make_coverage() {
  local species="$1"
  local sample="$2"
  local win="${RUN}/windows/${species}.50k.bed"
  local bed="${RUN}/mapping/${species}.${sample}.q20.bed"
  local out="${RUN}/coverage/${species}.${sample}.q20.50k.counts.tsv"
  if [[ -s "${out}" ]]; then
    log "Reuse existing ${out}"
    return
  fi
  log "Count ${species} ${sample} q20 reads in 50 kb windows"
  conda activate biosofeware
  bedtools coverage -a "${win}" -b "${bed}" -counts > "${out}"
  conda deactivate
}

for f in "${R108_FA}" "${MPO_FA}" "${R108_CENH3_R1}" "${R108_CENH3_R2}" "${R108_INPUT_R1}" "${R108_INPUT_R2}" "${MPO_CENH3_R1}" "${MPO_CENH3_R2}" "${MPO_INPUT_R1}" "${MPO_INPUT_R2}" "${SCRIPT}"; do
  require_file "${f}"
done

source ~/.bashrc 2>/dev/null || true

log "Prepare chromosome sizes and FASTA indexes"
conda activate biosofeware
samtools faidx "${R108_FA}"
samtools faidx "${MPO_FA}"
cut -f1,2 "${R108_FA}.fai" > "${RUN}/chrom_sizes/genome_R108.sizes"
cut -f1,2 "${MPO_FA}.fai" > "${RUN}/chrom_sizes/genome_Mpo.sizes"
bedtools makewindows -g "${RUN}/chrom_sizes/genome_R108.sizes" -w 50000 > "${RUN}/windows/genome_R108.50k.bed"
bedtools makewindows -g "${RUN}/chrom_sizes/genome_Mpo.sizes" -w 50000 > "${RUN}/windows/genome_Mpo.50k.bed"
conda deactivate

log "Build chromap indexes from raw FASTA when needed"
conda activate chromap
if [[ ! -s "${RUN}/index/genome_R108.chromap.index" ]]; then
  chromap -i -r "${R108_FA}" -o "${RUN}/index/genome_R108.chromap.index" 2> "${RUN}/logs/genome_R108.chromap_index.stderr"
fi
if [[ ! -s "${RUN}/index/genome_Mpo.chromap.index" ]]; then
  chromap -i -r "${MPO_FA}" -o "${RUN}/index/genome_Mpo.chromap.index" 2> "${RUN}/logs/genome_Mpo.chromap_index.stderr"
fi
conda deactivate

map_chip_bed "genome_R108" "CENH3" "${R108_FA}" "${RUN}/index/genome_R108.chromap.index" "${R108_CENH3_R1}" "${R108_CENH3_R2}"
map_chip_bed "genome_R108" "Input" "${R108_FA}" "${RUN}/index/genome_R108.chromap.index" "${R108_INPUT_R1}" "${R108_INPUT_R2}"
map_chip_bed "genome_Mpo" "CENH3" "${MPO_FA}" "${RUN}/index/genome_Mpo.chromap.index" "${MPO_CENH3_R1}" "${MPO_CENH3_R2}"
map_chip_bed "genome_Mpo" "Input" "${MPO_FA}" "${RUN}/index/genome_Mpo.chromap.index" "${MPO_INPUT_R1}" "${MPO_INPUT_R2}"

make_coverage "genome_R108" "CENH3"
make_coverage "genome_R108" "Input"
make_coverage "genome_Mpo" "CENH3"
make_coverage "genome_Mpo" "Input"

log "Call raw q20 CENH3 signal domains and R108 CEN5 region"
"${PYTHON_BIN}" "${SCRIPT}" call-signal --run "${RUN}" > "${RUN}/logs/call_signal.stdout" 2> "${RUN}/logs/call_signal.stderr"

log "Extract raw-called R108 CEN5 sequence"
conda activate biosofeware
bed_to_fasta "${RUN}/results/R108_CEN5_raw_q20_core.bed" "${R108_FA}" "${RUN}/fasta/R108_CEN5_raw_q20_core.fa"
bed_to_fasta "${RUN}/results/R108_CEN5_raw_q20_pm1Mb.bed" "${R108_FA}" "${RUN}/fasta/R108_CEN5_raw_q20_pm1Mb.fa"
bed_to_fasta "${RUN}/results/R108_CEN5_raw_q20_pm3Mb.bed" "${R108_FA}" "${RUN}/fasta/R108_CEN5_raw_q20_pm3Mb.fa"
conda deactivate

log "Fresh whole-genome minimap2 alignment: Mpo query to R108 reference"
conda activate minimap2_env
if [[ ! -s "${RUN}/paf/Mpo_to_R108.raw_asm5.paf" ]]; then
  minimap2 -x asm5 -c --cs -t "${THREADS}" "${R108_FA}" "${MPO_FA}" > "${RUN}/paf/Mpo_to_R108.raw_asm5.paf" 2> "${RUN}/logs/Mpo_to_R108.raw_asm5.stderr"
fi
log "Fresh CEN5 homology search: raw-called R108 CEN5 +/-3 Mb to Mpo"
if [[ ! -s "${RUN}/paf/R108_CEN5_pm3Mb_to_Mpo.raw_asm20.paf" ]]; then
  minimap2 -x asm20 -c --cs --secondary=yes -N 500 -t "${THREADS}" "${MPO_FA}" "${RUN}/fasta/R108_CEN5_raw_q20_pm3Mb.fa" > "${RUN}/paf/R108_CEN5_pm3Mb_to_Mpo.raw_asm20.paf" 2> "${RUN}/logs/R108_CEN5_pm3Mb_to_Mpo.raw_asm20.stderr"
fi
conda deactivate

log "Summarize raw minimap2 breakpoint and CEN5 relic evidence"
"${PYTHON_BIN}" "${SCRIPT}" summarize-paf --run "${RUN}" > "${RUN}/logs/summarize_paf.stdout" 2> "${RUN}/logs/summarize_paf.stderr"

log "Extract Mpo CEN5 relic candidates and active CENH3 domains for k-mer decay test"
conda activate biosofeware
bed_to_fasta "${RUN}/results/Mpo_CEN5_relic_candidate_merged_intervals.bed" "${MPO_FA}" "${RUN}/fasta/Mpo_CEN5_relic_candidates.fa"
awk '$6=="genome_Mpo"{print $1"\t"$2"\t"$3"\t"$4}' "${RUN}/results/raw_q20_called_CENH3_domains.bed" > "${RUN}/results/Mpo_raw_q20_active_CENH3_domains.bed"
bed_to_fasta "${RUN}/results/Mpo_raw_q20_active_CENH3_domains.bed" "${MPO_FA}" "${RUN}/fasta/Mpo_raw_q20_active_CENH3_domains.fa"
conda deactivate

log "Run CEN5-specific k-mer relic/decay test"
"${PYTHON_BIN}" "${SCRIPT}" kmer-decay --run "${RUN}" > "${RUN}/logs/kmer_decay.stdout" 2> "${RUN}/logs/kmer_decay.stderr"

log "Raw rerun finished"
