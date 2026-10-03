#!/usr/bin/env python3
import csv
import os


RUN_ROOT = "path/to/project/10.centromere_analysis/output_key"


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    ref = read_tsv(os.path.join(RUN_ROOT, "01_define_reference_CEN1_CEN8/results/R108_reference_CEN1_CEN8.tsv"))
    gene_proj = read_tsv(os.path.join(RUN_ROOT, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN1_CEN8_projected_to_Mpo.tsv"))
    peri_proj = read_tsv(os.path.join(RUN_ROOT, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN1_CEN8_sequence_block_projection_to_Mpo.pm3Mb.tsv"))
    core_proj = read_tsv(os.path.join(RUN_ROOT, "02_project_A17_R108_CENs_to_Mpo/results/R108_CEN1_CEN8_sequence_block_projection_to_Mpo.core_only.tsv"))
    reciprocal = read_tsv(os.path.join(RUN_ROOT, "04_reciprocal_project_Mpo_CENs_to_x8/Mpo_active_CENs_back_projected_to_R108.tsv"))

    ref_by_id = {r["ancestral_centromere_id"]: r for r in ref}
    gene_by_id = {r["ancestral_centromere_id"]: r for r in gene_proj}
    peri_by_id = {r["ancestral_centromere_id"]: r for r in peri_proj}
    core_by_id = {r["ancestral_centromere_id"]: r for r in core_proj}
    recip_by_cen = {}
    for r in reciprocal:
        cen = r["nearest_R108_CEN"].split("|")[0] if r.get("nearest_R108_CEN") else ""
        if cen:
            recip_by_cen.setdefault(cen, []).append(r)

    rows = []
    for cen_id in [f"CEN{i}" for i in range(1, 9)]:
        r = ref_by_id.get(cen_id, {})
        g = gene_by_id.get(cen_id, {})
        p = peri_by_id.get(cen_id, {})
        c = core_by_id.get(cen_id, {})
        recip_hits = recip_by_cen.get(cen_id, [])
        reciprocal_support = "yes" if recip_hits else "no"
        active_mpo = ";".join(x["Mpo_CENH3_core"] for x in recip_hits)
        if reciprocal_support == "yes":
            classification = "retained_active_centromere"
            interpretation = "At least one active Mpo CENH3 core reciprocally projects near this R108 ancestral centromere."
        else:
            classification = "missing_from_active_Mpo_CENH3_set"
            interpretation = "No active Mpo CENH3 core currently reciprocally projects near this R108 ancestral centromere; this is the primary candidate for centromere loss, inactivation or unresolved repositioning during the x=8 to x=7 transition."
        rows.append({
            "ancestral_centromere_id": cen_id,
            "R108_chr": r.get("reference_chr", ""),
            "R108_start": r.get("reference_start", ""),
            "R108_end": r.get("reference_end", ""),
            "A17_chr": r.get("A17_chr", ""),
            "A17_start": r.get("A17_start", ""),
            "A17_end": r.get("A17_end", ""),
            "reciprocal_active_Mpo_support": reciprocal_support,
            "reciprocal_active_Mpo_CENH3_core": active_mpo,
            "gene_anchor_projection": f"{g.get('Mpo_chr','')}:{g.get('Mpo_projected_start','')}-{g.get('Mpo_projected_end','')}",
            "gene_anchor_candidate_chr_summary": g.get("candidate_chr_summary", ""),
            "pericentromere_sequence_projection": f"{p.get('Mpo_chr','')}:{p.get('Mpo_projected_start','')}-{p.get('Mpo_projected_end','')}",
            "pericentromere_candidate_target_summary": p.get("candidate_target_summary", ""),
            "core_sequence_projection_class": c.get("classification", ""),
            "core_sequence_candidate_target_summary": c.get("candidate_target_summary", ""),
            "classification": classification,
            "interpretation": interpretation,
        })

    fields = [
        "ancestral_centromere_id", "R108_chr", "R108_start", "R108_end", "A17_chr", "A17_start", "A17_end",
        "reciprocal_active_Mpo_support", "reciprocal_active_Mpo_CENH3_core",
        "gene_anchor_projection", "gene_anchor_candidate_chr_summary",
        "pericentromere_sequence_projection", "pericentromere_candidate_target_summary",
        "core_sequence_projection_class", "core_sequence_candidate_target_summary",
        "classification", "interpretation",
    ]
    out = os.path.join(RUN_ROOT, "05_fusion_chr_CEN_fate/Mpo_ancestral_centromere_fate_integrated.tsv")
    write_tsv(out, rows, fields)
    print(out)


if __name__ == "__main__":
    main()
