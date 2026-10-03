#!/usr/bin/env python3
"""Run local tag-SNP conditional logistic models on the common complete-case set."""

import argparse
import csv
import gzip
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_candidate_association import logistic_irls, sample_order, pc_map


def dosage(gt):
    allele = gt.split(":", 1)[0]
    if allele in {".", "./.", ".|."}:
        return None
    tokens = allele.replace("|", "/").split("/")
    if len(tokens) != 2 or any(token not in {"0", "1"} for token in tokens):
        return None
    return int(tokens[0]) + int(tokens[1])


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def vcf_samples(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError("missing VCF sample header")


def query_tag_genotypes(tabix, vcf, chrom, pos, marker_id):
    samples = vcf_samples(vcf)
    completed = subprocess.run([tabix, vcf, "{}:{}-{}".format(chrom, pos, pos)], check=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 10 and int(fields[1]) == pos and fields[2] == marker_id:
            return {sample: dosage(gt) for sample, gt in zip(samples, fields[9:])}
    raise ValueError("tag marker is absent from input VCF")


def pvalue_from_beta(beta, se):
    if se <= 0.0 or not np.isfinite(se):
        return math.nan
    return math.erfc(abs(beta / se) / math.sqrt(2.0))


def model_result(name, outcome, predictors, predictor_names):
    x = np.column_stack([np.ones(len(outcome)), predictors])
    try:
        beta, covariance = logistic_irls(x, outcome)
        rows = []
        for index, predictor in enumerate(predictor_names, start=1):
            estimate = float(beta[index])
            se = float(math.sqrt(covariance[index, index]))
            rows.append({
                "model": name, "term": predictor, "n_complete_case": len(outcome),
                "beta": "{:.10g}".format(estimate), "se": "{:.10g}".format(se),
                "odds_ratio": "{:.10g}".format(math.exp(estimate)),
                "ci_low": "{:.10g}".format(math.exp(estimate - 1.96 * se)),
                "ci_high": "{:.10g}".format(math.exp(estimate + 1.96 * se)),
                "pvalue": "{:.10g}".format(pvalue_from_beta(estimate, se)), "status": "ok",
            })
        return rows
    except (ValueError, np.linalg.LinAlgError, OverflowError) as exc:
        return [{"model": name, "term": term, "n_complete_case": len(outcome), "beta": "NA", "se": "NA",
                 "odds_ratio": "NA", "ci_low": "NA", "ci_high": "NA", "pvalue": "NA",
                 "status": "not_estimable: {}".format(exc)} for term in predictor_names]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict-target", required=True)
    parser.add_argument("--sample-order", required=True)
    parser.add_argument("--pc-file", required=True)
    parser.add_argument("--snp-vcf", required=True)
    parser.add_argument("--tabix", required=True)
    parser.add_argument("--tag-chrom", default="4")
    parser.add_argument("--tag-pos", required=True, type=int)
    parser.add_argument("--tag-id", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    order = sample_order(args.sample_order)
    pcs = pc_map(args.pc_file, order)
    with open(args.strict_target, "r", newline="") as handle:
        target_rows = list(csv.DictReader(handle, delimiter="\t"))
    target = {row["sample"]: row for row in target_rows}
    if list(target) != order:
        raise ValueError("target sample order differs from analysis sample order")
    tag = query_tag_genotypes(args.tabix, args.snp_vcf, args.tag_chrom, args.tag_pos, args.tag_id)

    common = []
    for sample in order:
        row = target[sample]
        sv = None if row["genotype"] == "NA" else float(row["genotype"])
        if sv is None or tag.get(sample) is None:
            continue
        common.append((sample, float(row["phenotype"]), sv, float(tag[sample])))
    outcome = np.asarray([row[1] for row in common], dtype=float)
    sv = np.asarray([row[2] for row in common], dtype=float)
    tag_dosage = np.asarray([row[3] for row in common], dtype=float)
    pc_values = np.asarray([pcs[row[0]] for row in common], dtype=float)
    results = []
    results.extend(model_result("Model1_SV_unadjusted", outcome, sv[:, None], ["Chr23997_221bp_DEL"]))
    results.extend(model_result("Model2_tagSNP_unadjusted", outcome, tag_dosage[:, None], ["tag_SNP"]))
    results.extend(model_result("Model3_SV_plus_tagSNP_unadjusted", outcome, np.column_stack([sv, tag_dosage]),
                                ["Chr23997_221bp_DEL", "tag_SNP"]))
    results.extend(model_result("Model1_SV_plus_PC", outcome, np.column_stack([sv, pc_values]),
                                ["Chr23997_221bp_DEL", "PC1", "PC2", "PC3", "PC4", "PC5"]))
    results.extend(model_result("Model2_tagSNP_plus_PC", outcome, np.column_stack([tag_dosage, pc_values]),
                                ["tag_SNP", "PC1", "PC2", "PC3", "PC4", "PC5"]))
    results.extend(model_result("Model3_SV_plus_tagSNP_plus_PC", outcome,
                                np.column_stack([sv, tag_dosage, pc_values]),
                                ["Chr23997_221bp_DEL", "tag_SNP", "PC1", "PC2", "PC3", "PC4", "PC5"]))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out, ["model", "term", "n_complete_case", "beta", "se", "odds_ratio", "ci_low", "ci_high", "pvalue", "status"], results)


if __name__ == "__main__":
    main()
