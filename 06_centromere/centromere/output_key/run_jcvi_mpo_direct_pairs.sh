#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/10.centromere_analysis
SRC=${ROOT}/output_synthen/01_JCVI
OUT=${ROOT}/output_key/10_jcvi_mpo_synteny
THREADS=${THREADS:-20}

mkdir -p "${OUT}/pairwise_A17_Mpo" "${OUT}/pairwise_474_Mpo" "${OUT}/multi_species_Mpo_center"

copy_pair_inputs() {
    local d=$1
    local q=$2
    local s=$3
    local qbed=$4
    local qcds=$5
    local sbed=$6
    local scds=$7
    cp "${qbed}" "${d}/${q}.bed"
    cp "${qcds}" "${d}/${q}.cds"
    cp "${sbed}" "${d}/${s}.bed"
    cp "${scds}" "${d}/${s}.cds"
}

run_pair() {
    local d=$1
    local q=$2
    local s=$3
    (
        cd "${d}"
        if [[ ! -s "${q}.${s}.anchors.simple" ]]; then
            rm -f "${q}.${s}.last" "${q}.${s}.last.filtered" "${q}.${s}.anchors" "${q}.${s}.anchors.new" "${q}.${s}.anchors.simple"
            makeblastdb -in "${s}.cds" -dbtype nucl -out "${s}.cds"
            blastn \
                -query "${q}.cds" \
                -db "${s}.cds" \
                -out "${q}.${s}.blast" \
                -evalue 1e-5 \
                -outfmt 6 \
                -max_target_seqs 20 \
                -num_threads "${THREADS}"
            python -m jcvi.compara.blastfilter "${q}.${s}.blast" \
                --qbed="${q}.bed" \
                --sbed="${s}.bed" \
                --no_strip_names \
                --cscore=.7
            python -m jcvi.compara.synteny scan "${q}.${s}.blast.filtered" "${q}.${s}.anchors" \
                --qbed="${q}.bed" \
                --sbed="${s}.bed" \
                --no_strip_names \
                --dist=20
            python -m jcvi.compara.synteny screen --minspan=30 --simple "${q}.${s}.anchors" "${q}.${s}.anchors.new"
        fi
        cat > seqids <<'EOF'
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7
EOF
        cat > layout <<EOF
# y, xstart, xend, rotation, color, label, va, bed
.64, .12, .92, 0, #4C78A8, ${q}, top, ${q}.bed
.40, .12, .92, 0, #F58518, Mpo, top, ${s}.bed
# edges
e, 0, 1, ${q}.${s}.anchors.simple
EOF
        for fmt in pdf png svg; do
            python -m jcvi.graphics.karyotype seqids layout \
                --basepair \
                --shadestyle=line \
                --figsize=12x5 \
                --dpi=450 \
                --format="${fmt}" \
                --notex
            mv "karyotype.${fmt}" "${q}_vs_Mpo_JCVI_karyotype.${fmt}"
        done
    )
}

copy_pair_inputs \
    "${OUT}/pairwise_A17_Mpo" genome_A17 genome_Mpo \
    "${SRC}/genome_A17_genome_R108/genome_A17.bed" \
    "${SRC}/genome_A17_genome_R108/genome_A17.cds" \
    "${SRC}/R108_Mpo_synteny_lz10_run/genome_Mpo.bed" \
    "${SRC}/R108_Mpo_synteny_lz10_run/genome_Mpo.cds"

copy_pair_inputs \
    "${OUT}/pairwise_474_Mpo" genome_474 genome_Mpo \
    "${SRC}/genome_474.bed" \
    "${SRC}/genome_474.cds" \
    "${SRC}/R108_Mpo_synteny_lz10_run/genome_Mpo.bed" \
    "${SRC}/R108_Mpo_synteny_lz10_run/genome_Mpo.cds"

run_pair "${OUT}/pairwise_A17_Mpo" genome_A17 genome_Mpo
run_pair "${OUT}/pairwise_474_Mpo" genome_474 genome_Mpo

# Multi-species Mpo-centered summary using all direct Mpo anchors.
MS="${OUT}/multi_species_Mpo_center"
cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_R108.bed" "${MS}/"
cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_R108.genome_Mpo.anchors.simple" "${MS}/"
cp "${SRC}/genome_Msa_genome_Mpo_synteny_lz10_run/genome_Msa.bed" "${MS}/"
cp "${SRC}/genome_Msa_genome_Mpo_synteny_lz10_run/genome_Msa.genome_Mpo.anchors.simple" "${MS}/"
cp "${OUT}/pairwise_A17_Mpo/genome_A17.bed" "${MS}/"
cp "${OUT}/pairwise_A17_Mpo/genome_A17.genome_Mpo.anchors.simple" "${MS}/"
cp "${OUT}/pairwise_474_Mpo/genome_474.bed" "${MS}/"
cp "${OUT}/pairwise_474_Mpo/genome_474.genome_Mpo.anchors.simple" "${MS}/"
cp "${SRC}/R108_Mpo_synteny_lz10_run/genome_Mpo.bed" "${MS}/"

cat > "${MS}/seqids" <<'EOF'
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8
EOF

cat > "${MS}/layout" <<'EOF'
# y, xstart, xend, rotation, color, label, va, bed
.84, .08, .92, 0, #4C78A8, R108, top, genome_R108.bed
.68, .08, .92, 0, #72B7B2, A17, top, genome_A17.bed
.52, .08, .92, 0, #F58518, Mpo, top, genome_Mpo.bed
.36, .08, .92, 0, #54A24B, Msa, top, genome_Msa.bed
.20, .08, .92, 0, #B279A2, 474, top, genome_474.bed
# edges
e, 0, 2, genome_R108.genome_Mpo.anchors.simple
e, 1, 2, genome_A17.genome_Mpo.anchors.simple
e, 3, 2, genome_Msa.genome_Mpo.anchors.simple
e, 4, 2, genome_474.genome_Mpo.anchors.simple
EOF

(
    cd "${MS}"
    for fmt in pdf png svg; do
        python -m jcvi.graphics.karyotype seqids layout \
            --basepair \
            --shadestyle=line \
            --figsize=13x8 \
            --dpi=450 \
            --format="${fmt}" \
            --notex
        mv "karyotype.${fmt}" "Mpo_centered_R108_A17_Msa_474_JCVI_karyotype.${fmt}"
    done
)

find "${OUT}" -maxdepth 2 -type f | sort
