#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/10.centromere_analysis
SRC=${ROOT}/output_synthen/01_JCVI
OUT=${ROOT}/output_key/10_jcvi_mpo_synteny

mkdir -p "${OUT}"/global_R108_Mpo_Msa "${OUT}"/focus_CEN4_CEN5_MpoChr4 "${OUT}"/logs

prepare_global() {
    local d="${OUT}/global_R108_Mpo_Msa"
    cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_R108.bed" "${d}/"
    cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_Mpo.bed" "${d}/"
    cp "${SRC}/genome_Msa_genome_Mpo_synteny_lz10_run/genome_Msa.bed" "${d}/"
    cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_R108.genome_Mpo.anchors.simple" "${d}/"
    cp "${SRC}/genome_Msa_genome_Mpo_synteny_lz10_run/genome_Msa.genome_Mpo.anchors.simple" "${d}/"

    cat > "${d}/seqids" <<'EOF'
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
EOF

    cat > "${d}/layout" <<'EOF'
# y, xstart, xend, rotation, color, label, va, bed
.74, .08, .92, 0, #4C78A8, R108, top, genome_R108.bed
.52, .08, .92, 0, #F58518, Mpo, top, genome_Mpo.bed
.30, .08, .92, 0, #54A24B, Msa, top, genome_Msa.bed
# edges
e, 0, 1, genome_R108.genome_Mpo.anchors.simple
e, 2, 1, genome_Msa.genome_Mpo.anchors.simple
EOF
}

prepare_focus() {
    local d="${OUT}/focus_CEN4_CEN5_MpoChr4"
    cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_R108.bed" "${d}/"
    cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_Mpo.bed" "${d}/"
    cp "${SRC}/genome_Msa_genome_Mpo_synteny_lz10_run/genome_Msa.bed" "${d}/"
    cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_R108.genome_Mpo.anchors.simple" "${d}/"
    cp "${SRC}/genome_Msa_genome_Mpo_synteny_lz10_run/genome_Msa.genome_Mpo.anchors.simple" "${d}/"

    cat > "${d}/seqids" <<'EOF'
Chr4,Chr5
Chr4
Chr4,Chr5
EOF

    cat > "${d}/layout" <<'EOF'
# y, xstart, xend, rotation, color, label, va, bed
.74, .10, .90, 0, #4C78A8, R108 Chr4/Chr5, top, genome_R108.bed
.52, .10, .90, 0, #F58518, Mpo Chr4, top, genome_Mpo.bed
.30, .10, .90, 0, #54A24B, Msa Chr4/Chr5, top, genome_Msa.bed
# edges
e, 0, 1, genome_R108.genome_Mpo.anchors.simple
e, 2, 1, genome_Msa.genome_Mpo.anchors.simple
EOF
}

run_jcvi_plot() {
    local d=$1
    local prefix=$2
    (
        cd "${d}"
        for fmt in pdf png svg; do
            python -m jcvi.graphics.karyotype seqids layout \
                --basepair \
                --shadestyle=line \
                --figsize=12x7 \
                --dpi=450 \
                --format="${fmt}" \
                --notex
            mv "karyotype.${fmt}" "${prefix}.${fmt}"
        done
    )
}

prepare_global
prepare_focus
run_jcvi_plot "${OUT}/global_R108_Mpo_Msa" "Mpo_R108_Msa_global_JCVI_karyotype"
run_jcvi_plot "${OUT}/focus_CEN4_CEN5_MpoChr4" "Mpo_Chr4_R108_Msa_CEN4_CEN5_focus_JCVI_karyotype"

cat > "${OUT}/README.md" <<'EOF'
# Mpo-centric JCVI synteny plots

These plots were drawn with `python -m jcvi.graphics.karyotype` using existing direct JCVI anchors:

- `genome_R108.genome_Mpo.anchors.simple`
- `genome_Msa.genome_Mpo.anchors.simple`

The `global_R108_Mpo_Msa` panel shows whole-genome macrosynteny among R108, Mpo and Msa.
The `focus_CEN4_CEN5_MpoChr4` panel restricts the view to R108 Chr4/Chr5, Mpo Chr4 and Msa Chr4/Chr5, to inspect the syntenic context around the Mpo x=8 to x=7 chromosome-number reduction model.

Direct Mpo-A17 and Mpo-474 anchors were not present in the existing JCVI result folders at generation time; those can be added as separate pairwise runs if needed.
EOF

find "${OUT}" -maxdepth 2 -type f | sort
