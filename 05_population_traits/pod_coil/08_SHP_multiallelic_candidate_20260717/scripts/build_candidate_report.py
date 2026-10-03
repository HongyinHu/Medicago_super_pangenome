#!/usr/bin/env python3
"""Join audited functional alleles to low-frequency GMMAT score outputs."""

import argparse
import csv
from pathlib import Path


REPORT_FIELDS = [
    "site", "gene", "gene_class", "aa_change", "chrom", "pos", "ref", "alt",
    "analysis_called_samples", "alt_allele_count", "alt_frequency", "maf",
    "alt_case_carriers", "alt_control_carriers", "ref_case_samples", "ref_control_samples",
    "raw_fisher_p", "plink_snp", "pc1_5_n", "pc1_5_effect_allele", "pc1_5_af",
    "pc1_5_p", "pc1_2_n", "pc1_2_effect_allele", "pc1_2_af", "pc1_2_p",
    "interpretation",
]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_bim(path):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split()
            if len(fields) != 6:
                raise ValueError("Malformed BIM row: {0}".format(line.rstrip()))
            chrom, snp, _, pos, a1, a2 = fields
            rows.append({"chrom": chrom, "snp": snp, "pos": pos, "a1": a1, "a2": a2})
    return rows


def canonical_chrom(chrom):
    chrom = str(chrom)
    return chrom if chrom.startswith("Chr") else "Chr" + chrom


def find_plink_variant(functional_row, bim_rows):
    target_alleles = {functional_row["ref"].upper(), functional_row["alt"].upper()}
    matches = [
        row for row in bim_rows
        if canonical_chrom(row["chrom"]) == functional_row["chrom"]
        and row["pos"] == str(functional_row["pos"])
        and {row["a1"].upper(), row["a2"].upper()} == target_alleles
    ]
    if len(matches) != 1:
        raise ValueError(
            "Expected one PLINK record for {0} {1}:{2} {3}>{4}; found {5}".format(
                functional_row["site"], functional_row["chrom"], functional_row["pos"],
                functional_row["ref"], functional_row["alt"], len(matches)
            )
        )
    return matches[0]["snp"]


def score_by_snp(path):
    rows = read_tsv(path)
    scores = {row["SNP"]: row for row in rows}
    if len(scores) != len(rows):
        raise ValueError("Duplicate SNP IDs in score table: {0}".format(path))
    return scores


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def interpretation(pc5, pc2):
    if as_float(pc5["PVAL"]) < 0.05 and as_float(pc2["PVAL"]) < 0.05:
        return "supported in both structure-adjusted candidate models"
    return "not supported after structure correction"


def build_rows(functional_summary, bim, pc1_5, pc1_2):
    bim_rows = read_bim(bim)
    scores_pc5 = score_by_snp(pc1_5)
    scores_pc2 = score_by_snp(pc1_2)
    report_rows = []
    for row in read_tsv(functional_summary):
        snp = find_plink_variant(row, bim_rows)
        if snp not in scores_pc5 or snp not in scores_pc2:
            raise ValueError("Functional SNP missing from a GMMAT result: {0}".format(snp))
        pc5 = scores_pc5[snp]
        pc2 = scores_pc2[snp]
        result = {field: row.get(field, "") for field in REPORT_FIELDS}
        result.update({
            "plink_snp": snp,
            "pc1_5_n": pc5["N"], "pc1_5_effect_allele": pc5["A1"],
            "pc1_5_af": pc5["AF"], "pc1_5_p": pc5["PVAL"],
            "pc1_2_n": pc2["N"], "pc1_2_effect_allele": pc2["A1"],
            "pc1_2_af": pc2["AF"], "pc1_2_p": pc2["PVAL"],
            "interpretation": interpretation(pc5, pc2),
        })
        report_rows.append(result)
    return report_rows


def write_tsv(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=REPORT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path, rows):
    lines = [
        "# Candidate multiallelic association report",
        "",
        "Input: hard-filtered, non-biallelic SNP VCF; alleles were split before PLINK import.",
        "Association: 143 samples (103 spiral, 40 nonspiral), existing SNP GRM, MAF >= 0.005, two GMMAT models (PC1-5 and PC1-2).",
        "Raw Fisher P values are descriptive only; inference should use the structure-adjusted GMMAT P values.",
        "",
        "|Site|Amino-acid change|ALT MAF|Raw Fisher P|PC1-5 P|PC1-2 P|Conclusion|",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            "|{0}|{1}|{2}|{3}|{4}|{5}|{6}|".format(
                row["site"], row["aa_change"], row["maf"], row["raw_fisher_p"],
                row["pc1_5_p"], row["pc1_2_p"], row["interpretation"]
            )
        )
    lines.extend([
        "",
        "## Interpretation",
        "",
        "A functional allele can be retained in this candidate scan even when it was excluded from the standard MAF>=0.05 GWAS. A non-significant GMMAT score must not be presented as a GWAS peak; it instead motivates independent validation, within-lineage analysis, or a pre-specified gene/haplotype burden test.",
    ])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--functional-summary", required=True)
    parser.add_argument("--bim", required=True)
    parser.add_argument("--pc1-5", required=True)
    parser.add_argument("--pc1-2", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()
    rows = build_rows(args.functional_summary, args.bim, args.pc1_5, args.pc1_2)
    if not rows:
        raise ValueError("No functional alleles were available for reporting")
    write_tsv(args.out_tsv, rows)
    write_markdown(args.out_md, rows)
    print("Wrote {0} functional association rows".format(len(rows)))


if __name__ == "__main__":
    main()
