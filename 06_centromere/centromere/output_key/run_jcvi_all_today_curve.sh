#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/10.centromere_analysis/output_key

run_curve() {
    local rel_dir=$1
    local basename=$2
    local figsize=$3
    local d="${ROOT}/${rel_dir}"

    if [[ ! -f "${d}/seqids" || ! -f "${d}/layout" ]]; then
        echo "[SKIP] ${rel_dir}: missing seqids or layout" >&2
        return 0
    fi

    (
        cd "${d}"
        for fmt in pdf png svg; do
            python -m jcvi.graphics.karyotype seqids layout \
                --basepair \
                --shadestyle=curve \
                --figsize="${figsize}" \
                --dpi=450 \
                --format="${fmt}" \
                --notex
            mv "karyotype.${fmt}" "${basename}_curve.${fmt}"
        done
    )

    echo "[OK] ${rel_dir}/${basename}_curve.{pdf,png,svg}"
}

run_curve "10_jcvi_mpo_synteny/global_R108_Mpo_Msa" \
    "Mpo_R108_Msa_global_JCVI_karyotype" "12x7"

run_curve "10_jcvi_mpo_synteny/focus_CEN4_CEN5_MpoChr4" \
    "Mpo_Chr4_R108_Msa_CEN4_CEN5_focus_JCVI_karyotype" "12x7"

run_curve "10_jcvi_mpo_synteny/pairwise_A17_Mpo" \
    "genome_A17_vs_Mpo_JCVI_karyotype" "12x5"

run_curve "10_jcvi_mpo_synteny/pairwise_474_Mpo" \
    "genome_474_vs_Mpo_JCVI_karyotype" "12x5"

run_curve "10_jcvi_mpo_synteny/multi_species_Mpo_center" \
    "Mpo_centered_R108_A17_Msa_474_JCVI_karyotype" "13x8"

run_curve "11_jcvi_A17_R108_only" \
    "A17_R108_only_JCVI_karyotype" "12x5"

find "${ROOT}/10_jcvi_mpo_synteny" "${ROOT}/11_jcvi_A17_R108_only" \
    -maxdepth 2 -type f -name "*_curve.*" | sort
