#!/bin/bash
set -euo pipefail

# =========================
# User settings
# =========================

THREADS=10        # 每个物种 EDTA 使用线程数
PARALLEL_JOBS=12 # 同时跑 10 个物种

BASE_DIR=$(pwd)
GENOME_DIR=${BASE_DIR}/00_genome
SPECIES_LIST=${BASE_DIR}/species.list

EDTA_DIR=${BASE_DIR}/02_EDTA_single
LOG_DIR=${BASE_DIR}/logs_singleEDTA
CMD_FILE=${BASE_DIR}/run_all_singleEDTA_cmd.sh

# 如果有 CDS，填绝对路径；没有就留空
# CDS=path/to/project/related_species.cds.fa
CDS=""

mkdir -p "${EDTA_DIR}" "${LOG_DIR}"

rm -f "${CMD_FILE}"
touch "${CMD_FILE}"

# =========================
# Generate one run script per species
# =========================

while read SP
do
    [ -z "${SP}" ] && continue
    [[ "${SP}" =~ ^# ]] && continue

    GENOME=${GENOME_DIR}/${SP}.fa
    WORKDIR=${EDTA_DIR}/${SP}
    RUNSH=${WORKDIR}/run_${SP}.sh

    if [ ! -s "${GENOME}" ]; then
        echo "WARNING: genome not found, skip ${SP}: ${GENOME}"
        continue
    fi

    mkdir -p "${WORKDIR}"

    ln -sf "${GENOME}" "${WORKDIR}/${SP}.fa"

    cat > "${RUNSH}" <<EOF
#!/bin/bash
set -euo pipefail

source ~/.bashrc
conda activate EDTA_env

export PYTHONNOUSERSITE=1

SP="${SP}"
THREADS="${THREADS}"
WORKDIR="${WORKDIR}"
LOG_DIR="${LOG_DIR}"
CDS="${CDS}"

cd "\${WORKDIR}"

echo "========================================"
echo "Start species: \${SP}"
echo "Start time: \$(date)"
echo "Workdir: \${WORKDIR}"
echo "========================================"

# =========================
# 1. Run single EDTA
# =========================

if [ -s "\${SP}.fa.mod.EDTA.anno/\${SP}.fa.mod.out" ] && \\
   [ -s "\${SP}.fa.mod.EDTA.raw/LTR/\${SP}.fa.mod.pass.list" ]; then

    echo "EDTA result exists for \${SP}, skip EDTA."

else

    echo "Running single EDTA for \${SP}"

    if [ -n "\${CDS}" ]; then

        EDTA.pl \\
            --genome \${SP}.fa \\
            --species others \\
            --step all \\
            --sensitive 1 \\
            --anno 1 \\
            --evaluate 0 \\
            --overwrite 1 \\
            --cds \${CDS} \\
            --threads \${THREADS} \\
            > "\${LOG_DIR}/\${SP}.EDTA.log" 2>&1

    else

        EDTA.pl \\
            --genome \${SP}.fa \\
            --species others \\
            --step all \\
            --sensitive 1 \\
            --anno 1 \\
            --evaluate 0 \\
            --overwrite 1 \\
            --threads \${THREADS} \\
            > "\${LOG_DIR}/\${SP}.EDTA.log" 2>&1

    fi
fi

# =========================
# 2. Check EDTA output
# =========================

MOD_GENOME="\${WORKDIR}/\${SP}.fa.mod"
INTACT="\${WORKDIR}/\${SP}.fa.mod.EDTA.raw/LTR/\${SP}.fa.mod.pass.list"
ALL_LTR="\${WORKDIR}/\${SP}.fa.mod.EDTA.anno/\${SP}.fa.mod.out"

if [ ! -s "\${MOD_GENOME}" ]; then
    echo "ERROR: missing \${MOD_GENOME}"
    exit 1
fi

if [ ! -s "\${INTACT}" ]; then
    echo "ERROR: missing \${INTACT}"
    exit 1
fi

if [ ! -s "\${ALL_LTR}" ]; then
    echo "ERROR: missing \${ALL_LTR}"
    exit 1
fi

# =========================
# 3. Calculate LAI from single EDTA
# =========================

LAI_DIR="\${WORKDIR}/LAI"
mkdir -p "\${LAI_DIR}"

cd "\${LAI_DIR}"

ln -sf "\${MOD_GENOME}" "\${SP}.fa.mod"
ln -sf "\${INTACT}" "\${SP}.fa.mod.pass.list"
ln -sf "\${ALL_LTR}" "\${SP}.fa.mod.out"

if [ -s "\${SP}.fa.mod.out.LAI" ]; then
    echo "LAI result exists for \${SP}, skip LAI."
else
    echo "Calculating LAI for \${SP}"

    LAI \\
        -genome "\${SP}.fa.mod" \\
        -intact "\${SP}.fa.mod.pass.list" \\
        -all "\${SP}.fa.mod.out" \\
        -q \\
        -t "\${THREADS}" \\
        > "\${SP}.LAI.log" 2>&1
fi

echo "========================================"
echo "Finished species: \${SP}"
echo "End time: \$(date)"
echo "========================================"
EOF

    chmod +x "${RUNSH}"

    echo "bash ${RUNSH} > ${LOG_DIR}/${SP}.pipeline.log 2>&1" >> "${CMD_FILE}"

done < "${SPECIES_LIST}"

echo "========================================"
echo "All per-species scripts generated."
echo "Command file: ${CMD_FILE}"
echo "Run with:"
echo "parallel -j ${PARALLEL_JOBS} --joblog parallel.singleEDTA.joblog < ${CMD_FILE}"
echo "========================================"