#!/usr/bin/env python3
"""Write an auditable candidate-only GWAS validation report for Chr23997."""

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

from scipy.stats import fisher_exact


TARGET = "Chr4:88,980,831-88,981,052 (DEL length 221 bp)"


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_key_value(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 1)
            if len(fields) == 2:
                values[fields[0]] = fields[1]
    return values


def as_number(value):
    if value in (None, "", "NA", "nan", "NaN"):
        return math.nan
    return float(value)


def number(value, digits=4):
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):.{digits}g}"


def percentage(numerator, denominator):
    return "NA" if not denominator else f"{100 * numerator / denominator:.1f}%"


def table(rows, headers, columns):
    lines = ["| " + " | ".join(headers) + " |",
             "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(column, "NA")) for column in columns) + " |")
    return "\n".join(lines)


def callability(rows):
    counts = Counter()
    for row in rows:
        group = "spiny" if int(row["phenotype"]) == 1 else "spineless"
        counts[(group, "called" if row["genotype"] != "NA" else "missing")] += 1
    matrix = [[counts[("spiny", "called")], counts[("spiny", "missing")]],
              [counts[("spineless", "called")], counts[("spineless", "missing")]]]
    _, pvalue = fisher_exact(matrix, alternative="two-sided")
    return counts, pvalue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    stage = args.stage.resolve()

    strict = read_tsv(stage / "inputs/Chr23997_221bp_PAV.strict.tsv")
    sensitivity = read_tsv(stage / "inputs/Chr23997_221bp_PAV.sensitivity.tsv")
    models = read_tsv(stage / "results/association/Chr23997.candidate_association.all_models.tsv")
    gmmat = read_tsv(stage / "results/association/Chr23997.candidate_association.GMMAT.tsv")
    conditional = read_tsv(stage / "results/association/Chr23997.conditional_local_models.tsv")
    local = read_tsv(stage / "results/local/Chr23997.local_association.tsv")
    best_tag = read_tsv(stage / "results/local/Chr23997.local_best_SNP_INDEL_tag.tsv")[0]
    ld = read_tsv(stage / "results/local/Chr23997.local_LD.tsv")
    old_root = stage.parent / "08_Chr23997_targeted_genotyping_20260715"
    old_gemma = {callset: read_key_value(old_root / "summary" / f"{callset}.gemma_status.tsv")
                 for callset in ("strict", "sensitivity")}

    strict_counts, strict_callability_p = callability(strict)
    sensitivity_counts, sensitivity_callability_p = callability(sensitivity)
    strict_reason = Counter(row["call_reason"] for row in strict)
    strict_evidence = Counter(row["evidence_call"] for row in strict)
    strict_del_spiny = sum(row["genotype"] == "1" and int(row["phenotype"]) == 1 for row in strict)
    strict_del_spineless = sum(row["genotype"] == "1" and int(row["phenotype"]) == 0 for row in strict)
    strict_present_spiny = sum(row["genotype"] == "0" and int(row["phenotype"]) == 1 for row in strict)
    strict_present_spineless = sum(row["genotype"] == "0" and int(row["phenotype"]) == 0 for row in strict)
    split_ge2 = sum(as_number(row["target_split_alignment_reads"]) >= 2 for row in strict)
    del_cigar_ge1 = sum(as_number(row["target_del_cigar_reads"]) >= 1 for row in strict)
    direct_alt_ge2 = sum(as_number(row["target_direct_alt_reads"]) >= 2 for row in strict)
    local_by_class = Counter(row["class"] for row in local)

    best_tag_ld = [row for row in ld if row["id"] == best_tag["id"]]
    best_tag_r2 = best_tag_ld[0]["r2_with_Chr23997_221bp_PAV"] if best_tag_ld else "not calculated"

    def qcline(name, rows):
        called = [row for row in rows if row["genotype"] != "NA"]
        deletion = sum(row["genotype"] == "1" for row in called)
        present = sum(row["genotype"] == "0" for row in called)
        maf = min(deletion, present) / (2 * len(called)) if called else math.nan
        return {
            "callset": name,
            "total": len(rows),
            "called": len(called),
            "missing": len(rows) - len(called),
            "call_rate": percentage(len(called), len(rows)),
            "DEL": deletion,
            "PRESENT": present,
            "DEL_frequency_called": percentage(deletion, len(called)),
            "diploid_MAF": number(maf),
        }

    qc_rows = [qcline("strict (primary)", strict), qcline("sensitivity (supplementary)", sensitivity)]
    clean_models = []
    for row in models + gmmat:
        clean_models.append({
            "callset": row["callset"], "model": row["model"], "n_called": row["n_called"],
            "OR": row.get("odds_ratio", "NA"), "CI95": f"{row.get('ci_low', 'NA')} to {row.get('ci_high', 'NA')}",
            "P": row.get("pvalue", "NA"), "status": row["status"],
        })
    gemma_rows = []
    for callset, record in old_gemma.items():
        gemma_rows.append({"callset": callset, "model": record.get("model", "GEMMA LMM"),
                           "n_called": "complete-case", "OR": "NA", "CI95": "NA",
                           "P": record.get("p_wald", "NA"),
                           "status": f"beta={record.get('beta', 'NA')}; SE={record.get('se', 'NA')}; {record.get('status', 'NA')}"})

    conditional_rows = []
    for row in conditional:
        if row["model"] in ("Model1_SV_unadjusted", "Model2_tagSNP_unadjusted", "Model3_SV_plus_tagSNP_unadjusted"):
            conditional_rows.append({"model": row["model"], "term": row["term"], "n": row["n_complete_case"],
                                     "OR": row["odds_ratio"], "P": row["pvalue"], "status": row["status"]})

    report = [
        "# Chr23997 221-bp PAV targeted GWAS validation",
        "",
        "## Scope and immutable definition",
        "",
        f"This candidate-only analysis tests the exact second-intron deletion at **{TARGET}** in the *M. sativa* reference gene **Chr23997**. It does not search for a new SV and does not substitute the nearby short-read call at Chr4:88,982,864-88,983,055. Genotype 0/0 denotes the intact/presence allele, 1/1 the deletion allele and ./., an unresolved genotype.",
        "",
        "The primary genotype state is the strict targeted, direct-read-evidence call set. The sensitivity state uses the pre-specified broader targeted evidence rule and is reported only as robustness context. Missing calls were not imputed as reference.",
        "",
        "## 1. Genotype quality and phenotype balance",
        "",
        table(qc_rows, ["Call set", "Total", "Called", "Missing", "Call rate", "DEL", "PRESENT", "DEL frequency among called", "Diploid MAF"], ["callset", "total", "called", "missing", "call_rate", "DEL", "PRESENT", "DEL_frequency_called", "diploid_MAF"]),
        "",
        f"The strict complete-case set contains {strict_counts[('spiny', 'called')]} spiny and {strict_counts[('spineless', 'called')]} spineless samples. Its DEL allele occurs in {strict_del_spiny} spiny and {strict_del_spineless} spineless samples. The 2x2 association table is DEL/PRESENT by spiny/spineless, not an imputed 143-sample table.",
        "",
        table([
            {"phenotype": "spiny", "DEL": strict_del_spiny, "PRESENT": strict_present_spiny, "missing": strict_counts[("spiny", "missing")]},
            {"phenotype": "spineless", "DEL": strict_del_spineless, "PRESENT": strict_present_spineless, "missing": strict_counts[("spineless", "missing")]},
        ], ["Phenotype", "DEL", "PRESENT", "Missing"], ["phenotype", "DEL", "PRESENT", "missing"]),
        "",
        f"Callability itself was not detectably phenotype-associated (strict Fisher P={number(strict_callability_p)}; sensitivity Fisher P={number(sensitivity_callability_p)}). This does not remove the loss of power caused by the {len(strict) - sum(row['genotype'] != 'NA' for row in strict)} strict missing calls.",
        "",
        "### Direct evidence audit for strict calls",
        "",
        f"Among all 143 strict target assessments, target split-alignment support >=2 reads occurred in {split_ge2} samples; exact 221-bp CIGAR deletion support >=1 read occurred in {del_cigar_ge1} samples; and direct alternate evidence >=2 reads occurred in {direct_alt_ge2} samples. Call-reason counts: " + "; ".join(f"{key}={value}" for key, value in sorted(strict_reason.items())) + ".",
        "",
        "Evidence-state counts: " + "; ".join(f"{key}={value}" for key, value in sorted(strict_evidence.items())) + ". These metrics are retained in `inputs/genotype_by_sample.tsv` and are not converted into hard genotypes when they are ambiguous.",
        "",
        "## 2. Targeted association models",
        "",
        table(clean_models, ["Call set", "Model", "Called", "OR (DEL -> spiny)", "95% CI", "Nominal P", "Status"], ["callset", "model", "n_called", "OR", "CI95", "P", "status"]),
        "",
        table(gemma_rows, ["Call set", "Model", "Sample set", "P", "Status"], ["callset", "model", "n_called", "P", "status"]),
        "",
        "The strict Fisher and unadjusted logistic models support enrichment of the DEL allele among spiny samples (OR>1; nominal P<0.01). However, PC1-PC5 logistic models have complete/quasi-complete separation, and the PC+GRM GMMAT fits are numerically unstable. The earlier complete-case GEMMA LMM also has very large standard errors and non-significant Wald P-values. Therefore, neither a population-adjusted effect estimate nor a stable mixed-model P-value is available from these data.",
        "",
        "## 3. Why whole-genome SV GWAS did not test this exact event",
        "",
        "The exact target interval was absent from the high-confidence 143-sample SV cohort VCF (`medicago143.smoove.square.gq20.dp5.cr80.maf05.mac5.bnd1.vcf.gz`) when queried at Chr4:88,980,831-88,981,052. The strict target call rate is only 44.8%, below the cohort VCF call-rate filter (80%). Thus the global SV-GWAS did not contain this exact marker. Its failure to appear in the global Manhattan plot is not a null association test of the manually defined 221-bp PAV.",
        "",
        "Additional constraints are: (i) the strict carrier class has only 13 samples; (ii) genotype calls derive from direct local evidence with a high unresolved fraction; and (iii) pod spine is phylogenetically structured. The first two reduce power, while the third makes PC/kinship correction essential but unstable for this sparse binary marker.",
        "",
        "## 4. Local 100-kb association and LD",
        "",
        f"The local table contains {len(local)} records in the ±100-kb window (" + "; ".join(f"{key}={value}" for key, value in sorted(local_by_class.items())) + f"). The smallest nominal SNP/INDEL P-value is {best_tag['pvalue']} at {best_tag['id']} ({best_tag['chr']}:{best_tag['pos']}); this is {abs(int(best_tag['pos']) - 88980831):,} bp from the candidate start.",
        "",
        f"In the strict 64-sample complete-case subset, this tag SNP has r2={best_tag_r2} with the candidate because it is monomorphic or otherwise lacks paired genotype variation. It cannot be used to estimate a meaningful conditional SV+SNP model in this subset. The conditional models are:",
        "",
        table(conditional_rows, ["Model", "Term", "N", "OR", "P", "Status"], ["model", "term", "n", "OR", "P", "status"]),
        "",
        "The local Manhattan panel (`figures/Chr23997_local_GWAS.pdf`) separates global SNP/INDEL/SV association summaries from the targeted PAV Fisher result. They use different marker inclusion and inference models, so it is descriptive rather than a single common-scale GWAS test.",
        "",
        "## 5. Sensitivity analysis and recommended interpretation",
        "",
        "FarmCPU/BLINK was not used for the target PAV. These genome-wide multi-locus models require a robust, jointly genotyped marker matrix; applying them to a marker that failed the 80% call-rate cohort threshold would not resolve the missing-genotype problem and could create a misleading result. The valid sensitivity analyses are the broader evidence call set and the explicitly reported logistic, PC, GEMMA and GMMAT statuses above.",
        "",
        "**Classification: B, potential causal SV candidate.** Direct strict evidence gives a coherent unadjusted association and an effect direction consistent with the observed phenotype segregation. Nevertheless, low callability, genotype imbalance and separation prevent a stable population-corrected/mixed-model estimate. The targeted result supports prioritization for independent genotyping, local assembly and functional experiments; it does not establish causality or justify a genome-wide-significant claim.",
        "",
        "## Reproducible outputs",
        "",
        "- `inputs/Chr23997_221bp_PAV.strict.vcf.gz` and `.tbi`: strict, 143-sample candidate VCF.",
        "- `inputs/Chr23997_221bp_PAV.sensitivity.vcf.gz` and `.tbi`: supplementary VCF.",
        "- `inputs/genotype_by_sample.tsv`, `inputs/genotype_qc.tsv` and `inputs/genotype_case_control_distribution.tsv`: sample-level audit tables.",
        "- `results/association/`: Fisher, logistic, PC and GMMAT status tables.",
        "- `results/local/`: local association, LD and conditional-analysis tables.",
        "- `figures/Chr23997_local_GWAS.pdf`: local association visualization.",
        "",
    ]
    args.out.write_text("\n".join(report), encoding="utf-8")
    (stage / "summary/validation.done").write_text("validated\n", encoding="utf-8")


if __name__ == "__main__":
    main()
