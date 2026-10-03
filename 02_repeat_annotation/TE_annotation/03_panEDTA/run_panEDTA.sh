#!/usr/bin/env bash
set -eo pipefail
source ~/.bashrc
conda activate EDTA_env
export PYTHONNOUSERSITE=1
cd path/to/project/N_1.EDTA_single/03_panEDTA
THREADS=${THREADS:-20}
echo "Start panEDTA: $(date)"
echo "Host: $(hostname)"
echo "Workdir: $(pwd)"
echo "Threads: ${THREADS}"
echo "Genome count: $(wc -l < genome.list)"
which panEDTA.sh
bash path/to/project/N_1.EDTA_single/03_panEDTA/EDTA_runner/panEDTA.sh -g genome.list -f 3 -a 1 -o 0 -t "${THREADS}"
echo "Finish panEDTA: $(date)"
