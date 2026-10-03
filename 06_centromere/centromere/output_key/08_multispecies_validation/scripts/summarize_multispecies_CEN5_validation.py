#!/usr/bin/env python3
"""Summarize multispecies evidence for the Mpo ancestral CEN5 fate story."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cen_id_from_label(label: str) -> str:
    if not label:
        return ""
    return label.split("|", 1)[0]


def supported_set_from_mpo_recip(rows: list[dict[str, str]], col: str) -> set[str]:
    out = set()
    for row in rows:
        cid = row.get(col, "")
        if not cid:
            if "nearest_R108_CEN" in row:
                cid = cen_id_from_label(row.get("nearest_R108_CEN", ""))
            elif "nearest_Msa_CEN" in row:
                cid = cen_id_from_label(row.get("nearest_Msa_CEN", ""))
        if cid and row.get("classification", "").startswith("reciprocal_support_near"):
            out.add(cid)
    return out


def get_summary_row(rows: list[dict[str, str]], species: str) -> dict[str, str]:
    for row in rows:
        if row.get("species") == species:
            return row
    return {}


def find_row(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    for row in rows:
        if all(row.get(k) == v for k, v in criteria.items()):
            return row
    return {}


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default="path/to/project/10.centromere_analysis")
    args = parser.parse_args()

    root = Path(args.project_root)
    run = root / "output_key"
    out_dir = run / "08_multispecies_validation"
    script_dir = out_dir / "scripts"
    script_dir.mkdir(parents=True, exist_ok=True)

    r108_recip = read_tsv(run / "04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_R108.tsv")
    msa_recip = read_tsv(run / "04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_Msa_curated.pm3Mb.tsv")
    msa_to_mpo = read_tsv(run / "03_project_other_x8_CENs_to_Mpo/Msa_curated_CEN1_CEN8_sequence_projection_to_Mpo.pm3Mb.tsv")
    x8_core = read_tsv(run / "03_project_other_x8_CENs_to_Mpo/x8_active_CENH3_cores_to_R108_candidate_groups.core_only.summary.tsv")
    x8_pm3 = read_tsv(run / "03_project_other_x8_CENs_to_Mpo/x8_active_CENH3_cores_to_R108_candidate_groups.pm3Mb.summary.tsv")
    fusion = read_tsv(run / "05_fusion_chr_CEN_fate/Mpo_Chr4_CEN4_CEN5_fusion_context.tsv")
    coverage = read_tsv(run / "06_repeat_gene_mappability_controls/Mpo_CEN5_candidate_CENH3_Input_coverage_control.tsv")
    chromatin = read_tsv(run / "06_repeat_gene_mappability_controls/Mpo_CEN5_candidate_orthogonal_chromatin_features.tsv")

    all_cens = [f"CEN{i}" for i in range(1, 9)]
    r108_supported = supported_set_from_mpo_recip(r108_recip, "nearest_R108_CEN_id")
    msa_supported = supported_set_from_mpo_recip(msa_recip, "nearest_Msa_CEN_id")

    matrix_rows = [
        {
            "evidence_module": "Mpo active cores -> R108",
            "comparison": "Mpo x7 active CENH3 cores projected to R108 x8 ancestral CENs",
            "key_result": "supported=" + ",".join(c for c in all_cens if c in r108_supported) + "; missing=" + ",".join(c for c in all_cens if c not in r108_supported),
            "CEN5_status": "missing_from_primary_active_Mpo_CEN_set" if "CEN5" not in r108_supported else "supported",
            "support_level": "strong",
            "key_file": "04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_R108.tsv",
        },
        {
            "evidence_module": "Mpo active cores -> Msa curated",
            "comparison": "Mpo x7 active CENH3 cores projected to Msa x8 curated CEN domains",
            "key_result": "supported=" + ",".join(c for c in all_cens if c in msa_supported) + "; missing=" + ",".join(c for c in all_cens if c not in msa_supported),
            "CEN5_status": "missing_from_primary_active_Mpo_CEN_set" if "CEN5" not in msa_supported else "supported",
            "support_level": "strong",
            "key_file": "04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_Msa_curated.pm3Mb.tsv",
        },
    ]

    msa_cen5 = find_row(msa_to_mpo, centromere_id="CEN5")
    matrix_rows.append(
        {
            "evidence_module": "Msa curated CEN5 -> Mpo",
            "comparison": "Msa x8 curated CEN5 projected to Mpo x7 active CENH3 core set",
            "key_result": f"Msa CEN5 best target {msa_cen5.get('Mpo_chr','')}:{msa_cen5.get('Mpo_projected_start','')}-{msa_cen5.get('Mpo_projected_end','')}; {msa_cen5.get('classification','')}",
            "CEN5_status": "CEN5_sequence_projects_without_active_Mpo_core" if msa_cen5.get("classification") == "projects_without_active_Mpo_CEN" else msa_cen5.get("classification", ""),
            "support_level": "strong",
            "key_file": "03_project_other_x8_CENs_to_Mpo/Msa_curated_CEN1_CEN8_sequence_projection_to_Mpo.pm3Mb.tsv",
        }
    )

    for species, level, source in [
        ("R108", "strong", "self CENH3 core"),
        ("A17", "moderate", "core-only sequence projection to R108 CEN5"),
        ("474", "strong", "core-only sequence projection to R108 CEN5"),
        ("Msa", "moderate", "curated CEN5 plus Msa-Mpo reciprocal absence of Mpo active counterpart"),
    ]:
        core_row = get_summary_row(x8_core, species)
        pm3_row = get_summary_row(x8_pm3, species)
        matrix_rows.append(
            {
                "evidence_module": f"{species} x8 CEN5 retention",
                "comparison": f"{species} active CEN evidence relative to R108/Mpo",
                "key_result": f"core_CEN5_supported={core_row.get('CEN5_supported','NA')}; core_best={core_row.get('best_CEN5_support','')}; pm3_CEN5_supported={pm3_row.get('CEN5_supported','NA')}; pm3_best={pm3_row.get('best_CEN5_support','')}; source={source}",
                "CEN5_status": "active_CEN5_retained_in_x8_reference" if species != "Mpo" else "",
                "support_level": level,
                "key_file": "03_project_other_x8_CENs_to_Mpo/x8_active_CENH3_cores_to_R108_candidate_groups.*.summary.tsv",
            }
        )

    cen5_feature = find_row(fusion, feature="R108_CEN5_top_projection_cluster_on_Mpo_Chr4")
    active_feature = find_row(fusion, feature="active_Mpo_Chr4_CENH3_core")
    distance_feature = find_row(fusion, feature="distance_CEN5_top_cluster_to_active_Chr4_CENH3_core")
    matrix_rows.append(
        {
            "evidence_module": "Mpo CEN5 remnant",
            "comparison": "R108 CEN5 sequence-block cluster and active Mpo Chr4 CEN4-like core",
            "key_result": f"CEN5 cluster {cen5_feature.get('Mpo_chr','')}:{cen5_feature.get('start','')}-{cen5_feature.get('end','')}; active Chr4 core {active_feature.get('start','')}-{active_feature.get('end','')}; distance={distance_feature.get('note','')} bp",
            "CEN5_status": "local_CEN5_associated_sequence_remnant",
            "support_level": "strong_for_sequence_remnant_not_functional_loss",
            "key_file": "05_fusion_chr_CEN_fate/Mpo_Chr4_CEN4_CEN5_fusion_context.tsv",
        }
    )

    cov_cen5 = find_row(coverage, label="R108_CEN5_projection_cluster_top")
    cov_active = find_row(coverage, label="Mpo_active_CENH3_core_Chr4_primary_compare")
    chrom_cen5 = find_row(chromatin, label="R108_CEN5_projection_cluster_top")
    chrom_active = find_row(chromatin, label="Mpo_active_CENH3_core_Chr4_primary_compare")
    matrix_rows.append(
        {
            "evidence_module": "coverage and chromatin controls",
            "comparison": "CEN5 candidate interval versus active Mpo Chr4 CENH3 core",
            "key_result": f"q20_log2_CENH3/Input CEN5={cov_cen5.get('log2_CENH3_over_Input_q20_rpkm','')}, activeChr4={cov_active.get('log2_CENH3_over_Input_q20_rpkm','')}; q0 CEN5={cov_cen5.get('log2_CENH3_over_Input_q0_rpkm','')}, activeChr4={cov_active.get('log2_CENH3_over_Input_q0_rpkm','')}; CpG CEN5={chrom_cen5.get('CpG_methylation_50kb_bw_mean','')}, activeChr4={chrom_active.get('CpG_methylation_50kb_bw_mean','')}; TAD40 CEN5={chrom_cen5.get('TAD_score_40kb_mean','')}, activeChr4={chrom_active.get('TAD_score_40kb_mean','')}",
            "CEN5_status": "residual_CENH3_signal_with_mapping_stringency_caveat",
            "support_level": "control_supports_cautious_interpretation",
            "key_file": "06_repeat_gene_mappability_controls/Mpo_CEN5_candidate_*",
        }
    )

    species_rows = [
        {
            "species": "R108",
            "chromosome_number": 8,
            "active_CEN_count_used": 8,
            "CEN5_active_or_equivalent": "yes",
            "CEN5_evidence": get_summary_row(x8_core, "R108").get("best_CEN5_support", "self"),
            "Mpo_primary_active_counterpart": "no; Mpo active cores back to R108 miss CEN5",
            "interpretation": "x8 reference retains active CEN5; Mpo lacks a primary active CEN5 counterpart.",
        },
        {
            "species": "A17",
            "chromosome_number": 8,
            "active_CEN_count_used": 8,
            "CEN5_active_or_equivalent": "yes",
            "CEN5_evidence": get_summary_row(x8_core, "A17").get("best_CEN5_support", ""),
            "Mpo_primary_active_counterpart": "not direct; validates R108 CEN5 identity",
            "interpretation": "A17 active CEN5 core projects to R108 CEN5 at core level.",
        },
        {
            "species": "474",
            "chromosome_number": 8,
            "active_CEN_count_used": 8,
            "CEN5_active_or_equivalent": "yes",
            "CEN5_evidence": get_summary_row(x8_core, "474").get("best_CEN5_support", ""),
            "Mpo_primary_active_counterpart": "not direct; validates x8 CEN5 retention",
            "interpretation": "474 active CEN5 core strongly projects to R108 CEN5.",
        },
        {
            "species": "Msa",
            "chromosome_number": 8,
            "active_CEN_count_used": 8,
            "CEN5_active_or_equivalent": "yes",
            "CEN5_evidence": "curated CEN5 domain; pm3 R108 CEN5 candidate support; direct Msa CEN5-to-Mpo lacks active target",
            "Mpo_primary_active_counterpart": "no; Mpo active cores back to Msa miss CEN5",
            "interpretation": "Msa provides independent x8 comparison: CEN5 domain exists, but Mpo has no primary active counterpart.",
        },
        {
            "species": "Mpo",
            "chromosome_number": 7,
            "active_CEN_count_used": 7,
            "CEN5_active_or_equivalent": "no primary active CEN5",
            "CEN5_evidence": "local R108 CEN5-associated sequence cluster on Mpo Chr4 with residual/minor CENH3-associated signal",
            "Mpo_primary_active_counterpart": "not applicable",
            "interpretation": "Mpo likely retained CEN5-associated sequence locally but did not retain CEN5 as one of seven primary active CENH3 cores.",
        },
    ]

    write_tsv(out_dir / "multispecies_CEN5_validation_matrix.tsv", matrix_rows)
    write_tsv(out_dir / "multispecies_CEN5_species_summary.tsv", species_rows)

    summary_md = out_dir / "multispecies_CEN5_validation_summary.md"
    with summary_md.open("w") as handle:
        handle.write("# Multispecies CEN5 Validation Summary\n\n")
        handle.write("## Main Result\n\n")
        handle.write(
            "Across the available x=8 references and the x=7 Mpo genome, the data support a model in which ancestral CEN5 was not retained as a primary active CENH3 centromere in Mpo. This does not mean that the entire ancestral chromosome 5 was lost; rather, CEN5-associated sequence is locally retained near Mpo Chr4 while the primary active centromere set contains seven cores corresponding to CEN1, CEN2, CEN3, CEN4, CEN6, CEN7 and CEN8.\n\n"
        )
        handle.write("## Evidence\n\n")
        handle.write("- Mpo active cores projected to R108 support CEN1, CEN2, CEN3, CEN4, CEN6, CEN7 and CEN8, but not CEN5.\n")
        handle.write("- Mpo active cores projected to curated Msa CEN domains show the same pattern: CEN1, CEN2, CEN3, CEN4, CEN6, CEN7 and CEN8 are supported, but CEN5 is missing.\n")
        handle.write("- Curated Msa CEN5 projects to Mpo sequence without overlapping a primary active Mpo CENH3 core.\n")
        handle.write("- A17 and 474 independently retain active CEN5-like cores that project to R108 CEN5; 474 support is strong at core level, and A17 support is modest but direct.\n")
        handle.write("- The strongest R108 CEN5-associated Mpo sequence cluster is on Mpo Chr4 upstream of the active CEN4-like core, consistent with a local CEN5-associated remnant.\n\n")
        handle.write("## Writing Boundary\n\n")
        handle.write("Use cautious wording: ancestral CEN5 was not retained as a primary active centromere in Mpo; CEN5-associated sequence persists as a local remnant. Avoid saying the entire ancestral chromosome 5 was lost or that centromere inactivation is proven as the causal mechanism unless breakpoint-level and chromosome-painting evidence are added.\n\n")
        handle.write("## Key Outputs\n\n")
        handle.write("- `08_multispecies_validation/multispecies_CEN5_validation_matrix.tsv`\n")
        handle.write("- `08_multispecies_validation/multispecies_CEN5_species_summary.tsv`\n")
        handle.write("- `04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_Msa_curated.pm3Mb.tsv`\n")
        handle.write("- `03_project_other_x8_CENs_to_Mpo/Msa_curated_CEN1_CEN8_sequence_projection_to_Mpo.pm3Mb.tsv`\n")

    print(out_dir / "multispecies_CEN5_validation_matrix.tsv")
    print(out_dir / "multispecies_CEN5_species_summary.tsv")
    print(summary_md)


if __name__ == "__main__":
    main()
