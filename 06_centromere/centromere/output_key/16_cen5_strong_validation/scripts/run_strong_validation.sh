#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/10.centromere_analysis/output_key
RUN=${ROOT}/16_cen5_strong_validation
PREV=${ROOT}/15_cen5_relic_raw_rerun
CONDA=path/to/home/anaconda3/bin/conda
SCRIPT=${RUN}/scripts/analyze_strong_validation.py
THREADS=${THREADS:-32}

mkdir -p "${RUN}"/{01_junction,02_repeat_decay,03_x8_projection/{index,mapping,windows,coverage,signal,fasta,paf},logs,results,figures,tmp}

log() {
  echo "[$(date '+%F %T')] $*" | tee -a "${RUN}/logs/run_strong_validation.log"
}

require_file() {
  if [[ ! -s "$1" ]]; then
    echo "Required file missing or empty: $1" >&2
    exit 1
  fi
}

extract_region() {
  local fa="$1"
  local chrom="$2"
  local start0="$3"
  local end="$4"
  local name="$5"
  local out="$6"
  local start1=$((start0 + 1))
  "${CONDA}" run -n biosofeware samtools faidx "${fa}" "${chrom}:${start1}-${end}" \
    | awk -v n="${name}" 'BEGIN{done=0} /^>/{if(done==0){print ">"n; done=1}; next} {print}'
}

make_windows() {
  local species="$1"
  local fa="${ROOT}path/to/data/${species}.fa"
  local sizes="${RUN}/03_x8_projection/windows/${species}.sizes"
  local win="${RUN}/03_x8_projection/windows/${species}.50k.bed"
  if [[ ! -s "${win}" ]]; then
    "${CONDA}" run -n biosofeware samtools faidx "${fa}"
    cut -f1,2 "${fa}.fai" > "${sizes}"
    "${CONDA}" run -n biosofeware bedtools makewindows -g "${sizes}" -w 50000 > "${win}"
  fi
}

map_chip_bed() {
  local species="$1"
  local sample="$2"
  local fa="${ROOT}path/to/data/${species}.fa"
  local idx="${RUN}/03_x8_projection/index/${species}.chromap.index"
  local r1="${ROOT}path/to/data/${species}_${sample}_R1.fq.gz"
  local r2="${ROOT}path/to/data/${species}_${sample}_R2.fq.gz"
  local out="${RUN}/03_x8_projection/mapping/${species}.${sample}.q20.bed"
  require_file "${fa}"
  require_file "${r1}"
  require_file "${r2}"
  if [[ ! -s "${idx}" ]]; then
    log "Build chromap index for ${species}"
    "${CONDA}" run -n chromap chromap -i -r "${fa}" -o "${idx}" 2> "${RUN}/logs/${species}.chromap_index.stderr"
  fi
  if [[ ! -s "${out}" ]]; then
    log "Map ${species} ${sample} from reads with chromap MAPQ20"
    "${CONDA}" run -n chromap chromap --preset chip \
      -t "${THREADS}" \
      -x "${idx}" \
      -r "${fa}" \
      -1 "${r1}" \
      -2 "${r2}" \
      -q 20 \
      --remove-pcr-duplicates \
      --BED \
      -o "${out}.tmp" \
      --summary "${RUN}/logs/${species}.${sample}.q20.chromap.summary.txt" \
      2> "${RUN}/logs/${species}.${sample}.q20.chromap.stderr"
    "${CONDA}" run -n biosofeware sort -k1,1 -k2,2n "${out}.tmp" > "${out}"
    rm -f "${out}.tmp"
  fi
}

count_windows() {
  local species="$1"
  local sample="$2"
  local win="${RUN}/03_x8_projection/windows/${species}.50k.bed"
  local bed="${RUN}/03_x8_projection/mapping/${species}.${sample}.q20.bed"
  local out="${RUN}/03_x8_projection/coverage/${species}.${sample}.q20.50k.counts.tsv"
  if [[ ! -s "${out}" ]]; then
    log "Count ${species} ${sample} reads in 50 kb windows"
    "${CONDA}" run -n biosofeware bedtools coverage -a "${win}" -b "${bed}" -counts > "${out}"
  fi
}

log "Check required scripts and previous raw rerun outputs"
require_file "${SCRIPT}"
require_file "${PREV}/results/Mpo_CEN5_relic_focus_candidates.tsv"
require_file "${PREV}/results/raw_q20_called_CENH3_domains.tsv"
require_file "${PREV}/results/Mpo_raw_q20_active_CENH3_domains.bed"
require_file "${PREV}/results/R108_CEN5_raw_q20_core.bed"

log "01 junction refinement: extract raw FASTA local windows"
R108_FA="${ROOT}path/to/data"
MPO_FA="${ROOT}path/to/data"
require_file "${R108_FA}"
require_file "${MPO_FA}"
"${CONDA}" run -n biosofeware samtools faidx "${R108_FA}" Chr3 Chr5 Chr6 > "${RUN}/01_junction/R108_Chr3_Chr5_Chr6.fa"
: > "${RUN}/01_junction/Mpo_transition_windows.fa"
"${CONDA}" run -n biosofeware samtools faidx "${MPO_FA}" Chr3:21075816-23569156 >> "${RUN}/01_junction/Mpo_transition_windows.fa"
"${CONDA}" run -n biosofeware samtools faidx "${MPO_FA}" Chr3:38528302-41374862 >> "${RUN}/01_junction/Mpo_transition_windows.fa"
"${CONDA}" run -n biosofeware samtools faidx "${MPO_FA}" Chr5:24451337-27794443 >> "${RUN}/01_junction/Mpo_transition_windows.fa"
"${CONDA}" run -n biosofeware samtools faidx "${MPO_FA}" Chr5:54703849-57883318 >> "${RUN}/01_junction/Mpo_transition_windows.fa"

log "01 junction refinement: high-sensitivity minimap2 local alignment"
if [[ ! -s "${RUN}/01_junction/Mpo_transition_windows_vs_R108_Chr3_5_6.asm20.paf" ]]; then
  "${CONDA}" run -n minimap2_env minimap2 -x asm20 -c --cs --secondary=yes -N 1000 -p 0.05 -t "${THREADS}" \
    "${RUN}/01_junction/R108_Chr3_Chr5_Chr6.fa" \
    "${RUN}/01_junction/Mpo_transition_windows.fa" \
    > "${RUN}/01_junction/Mpo_transition_windows_vs_R108_Chr3_5_6.asm20.paf" \
    2> "${RUN}/logs/junction_minimap2_asm20.stderr"
fi
"${CONDA}" run -n biosofeware python "${SCRIPT}" summarize-junction --run "${RUN}"

log "02 targeted repeat decay: extract target regions from raw FASTA"
TARGET_TSV="${RUN}/02_repeat_decay/target_regions.tsv"
TARGET_FA="${RUN}/02_repeat_decay/target_regions.fa"
cat > "${TARGET_TSV}" <<'EOF'
target_id	species	class	chrom	start	end
R108_CEN5_core	genome_R108	active_reference_CEN5	Chr5	24000000	25700000
Mpo_Chr3_CEN5_relic	genome_Mpo	route_CEN5_relic	Chr3	20591033	22176458
Mpo_Chr5_CEN5_relic	genome_Mpo	weak_route_CEN5_relic	Chr5	25293760	25645386
Mpo_Chr4_CEN5_side_signal	genome_Mpo	side_CEN5_like_signal	Chr4	18314124	18432421
Mpo_Chr3_active_CENH3	genome_Mpo	active_CENH3	Chr3	46300000	46850000
Mpo_Chr5_active_CENH3	genome_Mpo	active_CENH3	Chr5	31850000	32250000
EOF
: > "${TARGET_FA}"
tail -n +2 "${TARGET_TSV}" | while IFS=$'\t' read -r target_id species class chrom start end; do
  extract_region "${ROOT}path/to/data/${species}.fa" "${chrom}" "${start}" "${end}" "${target_id}" "${TARGET_FA}" >> "${TARGET_FA}"
done

TRASH_LIB="${ROOT}path/to/data"
require_file "${TRASH_LIB}"
log "02 targeted repeat decay: BLAST R108 TRASH consensus to target regions"
if [[ ! -s "${RUN}/02_repeat_decay/blast_R108_TRASH_consensus_vs_targets.tsv" ]]; then
  "${CONDA}" run -n blast makeblastdb -in "${TARGET_FA}" -dbtype nucl -out "${RUN}/02_repeat_decay/target_regions.db" > "${RUN}/logs/repeat_makeblastdb.stdout" 2> "${RUN}/logs/repeat_makeblastdb.stderr"
  "${CONDA}" run -n blast blastn -task blastn -dust no -evalue 1e-5 -num_threads "${THREADS}" \
    -query "${TRASH_LIB}" \
    -db "${RUN}/02_repeat_decay/target_regions.db" \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
    > "${RUN}/02_repeat_decay/blast_R108_TRASH_consensus_vs_targets.tsv" \
    2> "${RUN}/logs/repeat_blast.stderr"
fi

log "02 targeted repeat decay: TRF tandem repeat scan"
if ! ls "${RUN}/02_repeat_decay"/target_regions.fa.*.dat >/dev/null 2>&1; then
  (
    cd "${RUN}/02_repeat_decay"
    timeout 600 "${CONDA}" run -n genome_repeat trf target_regions.fa 2 7 7 80 10 50 2000 -d -h > "${RUN}/logs/repeat_trf.stdout" 2> "${RUN}/logs/repeat_trf.stderr" || true
  )
else
  log "Reuse existing TRF dat file"
fi
"${CONDA}" run -n biosofeware python "${SCRIPT}" summarize-repeat --run "${RUN}"

log "03 x8 CENH3 projection: remap CENH3/Input reads and call domains"
X8_SPECIES="genome_A17,genome_R108,genome_474,genome_Msa"
for species in genome_A17 genome_R108 genome_474 genome_Msa; do
  make_windows "${species}"
  map_chip_bed "${species}" CENH3
  map_chip_bed "${species}" Input
  count_windows "${species}" CENH3
  count_windows "${species}" Input
done
"${CONDA}" run -n biosofeware python "${SCRIPT}" call-x8-domains --run "${RUN}" --species "${X8_SPECIES}"

log "03 x8 CEN5 projection: extract x=8 CEN5 +/-3 Mb and align to Mpo"
cp "${PREV}/results/Mpo_raw_q20_active_CENH3_domains.bed" "${RUN}/03_x8_projection/Mpo_active_CENH3_domains.bed"
awk 'BEGIN{OFS="\t"} NR>1 && (($1=="Chr3" && $2==20591033) || ($1=="Chr5" && $2==25293760) || ($1=="Chr4" && $2==18314124)){print $1,$2,$3,$1"_"$2"_"$3}' "${PREV}/results/Mpo_CEN5_relic_focus_candidates.tsv" > "${RUN}/03_x8_projection/Mpo_CEN5_relic_focus.bed"
echo -e "species\tchrom\tCEN5_start\tCEN5_end\tpm3Mb_start\tpm3Mb_end" > "${RUN}/03_x8_projection/x8_CEN5_pm3Mb_regions.tsv"
for species in genome_A17 genome_R108 genome_474 genome_Msa; do
  fa="${ROOT}path/to/data/${species}.fa"
  chr_len=$(awk '$1=="Chr5"{print $2}' "${fa}.fai")
  read start end < <(awk -v sp="${species}" 'BEGIN{FS=OFS="\t"} NR>1 && $1==sp && $3=="Chr5"{print $4,$5; exit}' "${RUN}/results/x8_raw_q20_called_CENH3_domains.tsv")
  if [[ -z "${start:-}" || -z "${end:-}" ]]; then
    echo "Could not call ${species} Chr5 CENH3" >&2
    exit 1
  fi
  pm_start=$((start - 3000000)); [[ "${pm_start}" -lt 0 ]] && pm_start=0
  pm_end=$((end + 3000000)); [[ "${pm_end}" -gt "${chr_len}" ]] && pm_end="${chr_len}"
  echo -e "${species}\tChr5\t${start}\t${end}\t${pm_start}\t${pm_end}" >> "${RUN}/03_x8_projection/x8_CEN5_pm3Mb_regions.tsv"
  extract_region "${fa}" Chr5 "${pm_start}" "${pm_end}" "${species}_CEN5_pm3Mb" "${RUN}/03_x8_projection/fasta/${species}_CEN5_pm3Mb.fa" > "${RUN}/03_x8_projection/fasta/${species}_CEN5_pm3Mb.fa"
  "${CONDA}" run -n minimap2_env minimap2 -x asm20 -c --cs --secondary=yes -N 500 -t "${THREADS}" \
    "${MPO_FA}" "${RUN}/03_x8_projection/fasta/${species}_CEN5_pm3Mb.fa" \
    > "${RUN}/03_x8_projection/paf/${species}_CEN5_pm3Mb_to_Mpo.asm20.paf" \
    2> "${RUN}/logs/${species}_CEN5_pm3Mb_to_Mpo.stderr"
done
"${CONDA}" run -n biosofeware python "${SCRIPT}" summarize-x8-projection --run "${RUN}"

log "Strong validation pipeline finished"
