#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/10.centromere_analysis
SRC=${ROOT}/output_synthen/01_JCVI/genome_A17_genome_R108
OUT=${ROOT}/output_key/11_jcvi_A17_R108_only

mkdir -p "${OUT}"

cp "${SRC}/genome_A17.bed" "${OUT}/"
cp "${SRC}/genome_R108.bed" "${OUT}/"
cp "${SRC}/genome_A17.genome_R108.anchors.simple" "${OUT}/"

cat > "${OUT}/seqids" <<'EOF'
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
EOF

cat > "${OUT}/layout" <<'EOF'
# y, xstart, xend, rotation, color, label, va, bed
.64, .08, .92, 0, #72B7B2, A17, top, genome_A17.bed
.40, .08, .92, 0, #4C78A8, R108, top, genome_R108.bed
# edges
e, 0, 1, genome_A17.genome_R108.anchors.simple
EOF

(
    cd "${OUT}"
    for fmt in pdf png svg; do
        python -m jcvi.graphics.karyotype seqids layout \
            --basepair \
            --shadestyle=curve \
            --figsize=12x5 \
            --dpi=450 \
            --format="${fmt}" \
            --notex
        mv "karyotype.${fmt}" "A17_R108_only_JCVI_karyotype_curve.${fmt}"
    done
)

cat > "${OUT}/README.md" <<'EOF'
# A17-R108 only JCVI karyotype

This directory contains a two-track JCVI karyotype plot using only A17 and R108.
No other species or custom CENH3 tracks are included.

Inputs:
- genome_A17.bed
- genome_R108.bed
- genome_A17.genome_R108.anchors.simple

Command:
python -m jcvi.graphics.karyotype seqids layout --basepair --shadestyle=line --figsize=12x5 --dpi=450 --format pdf/png/svg --notex
EOF

find "${OUT}" -maxdepth 1 -type f | sort
