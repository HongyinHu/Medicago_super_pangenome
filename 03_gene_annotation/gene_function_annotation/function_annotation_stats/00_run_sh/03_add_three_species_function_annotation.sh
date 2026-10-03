#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="${WORK_DIR:-path/to/project/7.function_anno_genes/function_annotation_stats}"
SPECIES_PEP_DIR="${SPECIES_PEP_DIR:-path/to/project/7.function_anno_genes/data/species_pep}"
OUTPUT_DIR="${WORK_DIR}/output"
LOG_DIR="${WORK_DIR}/01_logs"
CMD_DIR="${WORK_DIR}/02_cmds"
LOCK_DIR="${WORK_DIR}/04_locks"
DONE_DIR="${WORK_DIR}/05_done"
RUN_WITH_LOCK="${WORK_DIR}/run_function_task_with_lock.sh"

SAMPLES=(genome_Mar genome_Mru genome_zm4)

CONDA="path/to/home/anaconda3/bin/conda"
PFAM_SCAN="path/to/home/anaconda3/envs/pfam/bin/pfam_scan.pl"
INTERPROSCAN="path/to/home/sofeware/my_interproscan/interproscan-5.30-69.0/interproscan.sh"
EGGNOG_DB_DIR="path/to/home/db/eggnog"

PFAM_DB="path/to/data/database/pfam"
NR_DB="path/to/project/7.function_anno_genes/data/nr.gz.dmnd"
SWISS_DB="path/to/project/7.function_anno_genes/data/uniprot_sprot.fasta.gz.dmnd"
TREMBL_DB="path/to/project/7.function_anno_genes/data/uniprot_trembl.fasta.dmnd"
KOG_DB="path/to/data/database/KOG_database/KOG_dataset.fa.dmnd"

mkdir -p "${LOG_DIR}" "${CMD_DIR}" "${LOCK_DIR}" "${DONE_DIR}" "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"/{pfam_anno,InterproScan_anno,NR_anno,SwissPort_anno,KOG_anno,eggNOG_anno,TrEMBL_anno,KEGG_anno,GO_anno}

normalize_fasta_headers() {
  local input_fasta="$1"
  local output_fasta="$2"
  awk '
    /^>/ {
      header=$1
      sub(/^>/, "", header)
      print ">" header
      next
    }
    {
      seq=toupper($0)
      gsub(/[^ACDEFGHIKLMNPQRSTVWYBXZJUO*]/, "X", seq)
      print seq
    }
  ' "${input_fasta}" > "${output_fasta}"
}

clean_previous_failed_outputs() {
  for sample in "${SAMPLES[@]}"; do
    find "${LOCK_DIR}" -type f -name "${sample}.lock" -delete 2>/dev/null || true
    find "${OUTPUT_DIR}" -type f \
      \( -path "*${sample}*" \) \
      \( -name "*.pfam.out" -o -name "*.interproscan.out" -o -name "*_blastp.tab" -o -name "*.emapper.annotations" -o -name "*.emapper.hits" -o -name "*.emapper.seed_orthologs" \) \
      -delete 2>/dev/null || true
  done
}

prepare_inputs() {
  : > "${LOG_DIR}/add_species.prepare_inputs.log"
  for sample in "${SAMPLES[@]}"; do
    pep="${SPECIES_PEP_DIR}/${sample}.pep"
    if [[ ! -s "${pep}" ]]; then
      echo "[ERROR] missing protein FASTA: ${pep}" | tee -a "${LOG_DIR}/add_species.prepare_inputs.log" >&2
      exit 1
    fi
    for folder in pfam_anno InterproScan_anno NR_anno SwissPort_anno KOG_anno eggNOG_anno TrEMBL_anno; do
      mkdir -p "${OUTPUT_DIR}/${folder}/${sample}"
      normalize_fasta_headers "${pep}" "${OUTPUT_DIR}/${folder}/${sample}/${sample}.pep.format"
    done
    count="$(grep -c '^>' "${OUTPUT_DIR}/InterproScan_anno/${sample}/${sample}.pep.format")"
    printf '[OK] %s proteins=%s source=%s\n' "${sample}" "${count}" "${pep}" | tee -a "${LOG_DIR}/add_species.prepare_inputs.log"
  done
}

write_commands() {
  : > "${CMD_DIR}/add_species.pfam.cmds"
  : > "${CMD_DIR}/add_species.interproscan.cmds"
  : > "${CMD_DIR}/add_species.diamond.cmds"
  : > "${CMD_DIR}/add_species.eggnog.cmds"

  for sample in "${SAMPLES[@]}"; do
    echo "${RUN_WITH_LOCK} ${DONE_DIR}/pfam/${sample}.done ${LOCK_DIR}/pfam/${sample}.lock bash -lc 'cd ${OUTPUT_DIR}/pfam_anno && rm -f ${sample}/${sample}.pep.format.pfam.out && ${CONDA} run -n pfam ${PFAM_SCAN} -fasta ${sample}/${sample}.pep.format -dir ${PFAM_DB} -outfile ${sample}/${sample}.pep.format.pfam.out -cpu 30'" >> "${CMD_DIR}/add_species.pfam.cmds"
    echo "${RUN_WITH_LOCK} ${DONE_DIR}/interproscan/${sample}.done ${LOCK_DIR}/interproscan/${sample}.lock bash -lc 'cd ${OUTPUT_DIR}/InterproScan_anno && rm -f ${sample}/${sample}.pep.format.interproscan.out && ${CONDA} run -n java ${INTERPROSCAN} -f tsv -i ${sample}/${sample}.pep.format -cpu 30 --highmem -o ${sample}/${sample}.pep.format.interproscan.out -iprlookup -goterms -pa -td ${sample}/temp'" >> "${CMD_DIR}/add_species.interproscan.cmds"

    echo "${RUN_WITH_LOCK} ${DONE_DIR}/NR/${sample}.done ${LOCK_DIR}/NR/${sample}.lock bash -lc 'cd ${OUTPUT_DIR}/NR_anno && rm -f ${sample}/${sample}.pep.format_blastp.tab && ${CONDA} run -n biosofeware diamond blastp --threads 30 --db ${NR_DB} --query ${sample}/${sample}.pep.format --out ${sample}/${sample}.pep.format_blastp.tab --outfmt 6 --sensitive --evalue 1e-5 --quiet'" >> "${CMD_DIR}/add_species.diamond.cmds"
    echo "${RUN_WITH_LOCK} ${DONE_DIR}/SwissPort/${sample}.done ${LOCK_DIR}/SwissPort/${sample}.lock bash -lc 'cd ${OUTPUT_DIR}/SwissPort_anno && rm -f ${sample}/${sample}.pep.format_blastp.tab && ${CONDA} run -n biosofeware diamond blastp --threads 30 --db ${SWISS_DB} --query ${sample}/${sample}.pep.format --out ${sample}/${sample}.pep.format_blastp.tab --outfmt 6 --sensitive --evalue 1e-5 --quiet'" >> "${CMD_DIR}/add_species.diamond.cmds"
    echo "${RUN_WITH_LOCK} ${DONE_DIR}/KOG/${sample}.done ${LOCK_DIR}/KOG/${sample}.lock bash -lc 'cd ${OUTPUT_DIR}/KOG_anno && rm -f ${sample}/${sample}.pep.format_blastp.tab && ${CONDA} run -n biosofeware diamond blastp --threads 30 --db ${KOG_DB} --query ${sample}/${sample}.pep.format --out ${sample}/${sample}.pep.format_blastp.tab --outfmt 6 --sensitive --evalue 1e-5 --quiet'" >> "${CMD_DIR}/add_species.diamond.cmds"
    echo "${RUN_WITH_LOCK} ${DONE_DIR}/TrEMBL/${sample}.done ${LOCK_DIR}/TrEMBL/${sample}.lock bash -lc 'cd ${OUTPUT_DIR}/TrEMBL_anno && rm -f ${sample}/${sample}.pep.format_blastp.tab && ${CONDA} run -n biosofeware diamond blastp --threads 30 --db ${TREMBL_DB} --query ${sample}/${sample}.pep.format --out ${sample}/${sample}.pep.format_blastp.tab --outfmt 6 --sensitive --evalue 1e-5 --quiet'" >> "${CMD_DIR}/add_species.diamond.cmds"

    echo "${RUN_WITH_LOCK} ${DONE_DIR}/eggNOG/${sample}.done ${LOCK_DIR}/eggNOG/${sample}.lock bash -lc 'cd ${OUTPUT_DIR}/eggNOG_anno && rm -f ${sample}/${sample}.emapper.annotations && ${CONDA} run -n eggnog_env emapper.py -i ${sample}/${sample}.pep.format --itype proteins -m diamond --data_dir ${EGGNOG_DB_DIR} --cpu 30 --override --output ${sample} --output_dir ${sample}'" >> "${CMD_DIR}/add_species.eggnog.cmds"
  done
}

write_review_scripts() {
  mkdir -p "${WORK_DIR}/00_run_sh"
  cp -f "$0" "${WORK_DIR}/00_run_sh/03_add_three_species_function_annotation.sh"
  for f in "${CMD_DIR}"/add_species.*.cmds; do
    b="$(basename "${f}" .cmds)"
    {
      echo "#!/usr/bin/env bash"
      echo "set -euo pipefail"
      echo "# Command list captured for review: ${f}"
      cat "${f}"
    } > "${WORK_DIR}/00_run_sh/${b}.sh"
    chmod +x "${WORK_DIR}/00_run_sh/${b}.sh"
  done
}

start_queue() {
  local node="$1"
  local name="$2"
  local jobs="$3"
  local cmd_file="$4"
  ssh "${node}" "cd ${WORK_DIR} && nohup bash run_command_queue_no_parallel.sh ${jobs} ${cmd_file} ${LOG_DIR}/${name}.${node}.stdout ${LOG_DIR}/${name}.${node}.stderr >/dev/null 2>&1 & echo \$! > ${name}.${node}.pid"
  echo "[STARTED] ${name} on ${node} pid=$(cat "${WORK_DIR}/${name}.${node}.pid") jobs=${jobs}"
}

clean_previous_failed_outputs
prepare_inputs
write_commands
write_review_scripts

start_queue lz10 add_species_pfam 2 "${CMD_DIR}/add_species.pfam.cmds"
start_queue lz10 add_species_interproscan 1 "${CMD_DIR}/add_species.interproscan.cmds"
start_queue lz31 add_species_diamond 3 "${CMD_DIR}/add_species.diamond.cmds"
start_queue lz31 add_species_eggnog 1 "${CMD_DIR}/add_species.eggnog.cmds"
