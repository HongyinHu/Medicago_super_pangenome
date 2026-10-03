#!/usr/bin/env python3
"""ALT-carrier Fisher screen and plotting for the supplementary SV analysis."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact


def parse_gt(value):
    """Return ALT dosage for a biallelic GT, or None for a missing/invalid call."""
    gt = value.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return None
    return sum(int(allele) for allele in alleles)


def lambda_gc(pvalues):
    values = np.asarray([value for value in pvalues if 0 < value <= 1], dtype=float)
    if len(values) == 0:
        return float("nan")
    return float(np.median(chi2.isf(values, 1)) / 0.4549364)


def markdown_table(frame):
    table = frame.fillna("NA").astype(str)
    columns = list(table.columns)
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    lines.extend("| " + " | ".join(row) + " |" for row in table.itertuples(index=False, name=None))
    return "\n".join(lines)


def read_phenotype(path):
    frame = pd.read_csv(path, sep="\t")
    if not {"id", "y"}.issubset(frame.columns):
        raise ValueError(f"Expected id/y columns in {path}; found {list(frame.columns)}")
    values = dict(zip(frame["id"].astype(str), frame["y"].astype(int)))
    if set(values.values()) - {0, 1}:
        raise ValueError("Phenotype must be encoded 0/1")
    return values


def read_lines(path):
    with open(path, encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def run_fisher(
    query_path,
    all_sample_path,
    phenotype_path,
    output_path,
    min_maf,
    max_initial_missing,
    max_final_missing,
):
    samples = read_lines(all_sample_path)
    phenotype = read_phenotype(phenotype_path)
    final_indices = [index for index, sample in enumerate(samples) if sample in phenotype]
    if len(final_indices) != len(phenotype):
        observed = {samples[index] for index in final_indices}
        absent = sorted(set(phenotype) - observed)
        raise ValueError(f"Phenotype samples absent from VCF: {absent[:10]}")
    labels = np.array([phenotype[samples[index]] for index in final_indices], dtype=int)
    records = []
    with open(query_path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5 + len(samples):
                raise ValueError(f"Malformed query row {line_number}: expected {5 + len(samples)} fields")
            dosage_all = np.array([parse_gt(value) for value in fields[5:]], dtype=object)
            called_all = np.array([value is not None for value in dosage_all], dtype=bool)
            n_initial_called = int(np.sum(called_all))
            n_initial_missing = int(np.sum(~called_all))
            if n_initial_called == 0 or n_initial_missing / len(samples) > max_initial_missing:
                continue
            dosage = dosage_all[final_indices]
            called = np.array([value is not None for value in dosage], dtype=bool)
            n_called = int(np.sum(called))
            n_missing = int(np.sum(~called))
            if n_called == 0 or n_missing / len(final_indices) > max_final_missing:
                continue
            alt_count = int(sum(value for value in dosage if value is not None))
            allele_number = 2 * n_called
            maf = min(alt_count / allele_number, 1 - alt_count / allele_number)
            if maf < min_maf:
                continue
            case = labels == 1
            control = labels == 0
            alt_carrier = np.array([False if value is None else value > 0 for value in dosage], dtype=bool)
            case_alt = int(np.sum(case & called & alt_carrier))
            case_ref = int(np.sum(case & called & ~alt_carrier))
            control_alt = int(np.sum(control & called & alt_carrier))
            control_ref = int(np.sum(control & called & ~alt_carrier))
            odds_ratio, pvalue = fisher_exact(
                [[case_alt, case_ref], [control_alt, control_ref]], alternative="two-sided"
            )
            records.append(
                {
                    "CHROM": fields[0],
                    "POS": int(fields[1]),
                    "SV_ID": fields[2],
                    "REF": fields[3],
                    "ALT": fields[4],
                    "N_nonmissing": n_called,
                    "N_missing": n_missing,
                    "N_initial_nonmissing": n_initial_called,
                    "N_initial_missing": n_initial_missing,
                    "call_rate": n_called / len(samples),
                    "ALT_allele_count": alt_count,
                    "MAF": maf,
                    "case_ALT_carrier": case_alt,
                    "case_REF_only": case_ref,
                    "control_ALT_carrier": control_alt,
                    "control_REF_only": control_ref,
                    "ALT_carrier_OR": odds_ratio,
                    "Fisher_P": pvalue,
                    "minus_log10_Fisher_P": -math.log10(max(pvalue, np.finfo(float).tiny)),
                }
            )
    return pd.DataFrame(records)


def chromosome_order(values):
    def key(value):
        text = str(value)
        digits = "".join(char for char in text if char.isdigit())
        return (int(digits) if digits else 10**9, text)

    return sorted(pd.unique(values), key=key)


def manhattan(frame, p_column, output_prefix, title):
    data = frame.dropna(subset=[p_column, "POS"]).copy()
    data = data[(data[p_column] > 0) & (data[p_column] <= 1)]
    data["POS"] = pd.to_numeric(data["POS"])
    chromosomes = chromosome_order(data["CHROM"])
    offset = 0
    ticks, labels, chunks = [], [], []
    for index, chromosome in enumerate(chromosomes):
        chunk = data[data["CHROM"].astype(str) == str(chromosome)].sort_values("POS").copy()
        if chunk.empty:
            continue
        chunk["x"] = chunk["POS"] + offset
        offset = int(chunk["x"].max()) + 1_000_000
        ticks.append(float(chunk["x"].median()))
        labels.append(str(chromosome).replace("Chr", ""))
        chunk["colour"] = "#78bce6" if index % 2 == 0 else "#e84a9b"
        chunks.append(chunk)
    plot = pd.concat(chunks, ignore_index=True)
    fig, axis = plt.subplots(figsize=(11, 4.6))
    axis.scatter(plot["x"], -np.log10(plot[p_column]), c=plot["colour"], s=2.4, linewidths=0, rasterized=True)
    axis.axhline(-math.log10(0.05), linestyle="--", color="#777777", linewidth=0.8, label="P = 0.05")
    axis.axhline(-math.log10(0.01), linestyle="--", color="#333333", linewidth=0.8, label="P = 0.01")
    axis.set_xticks(ticks)
    axis.set_xticklabels(labels)
    axis.set_xlabel("Chromosome")
    axis.set_ylabel(r"$-\log_{10}(P)$")
    axis.set_title(title)
    axis.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    for suffix in (".pdf", ".png"):
        fig.savefig(str(output_prefix) + suffix, dpi=360, bbox_inches="tight")
    plt.close(fig)


def qq(frame, p_column, output_prefix, title):
    values = np.asarray(frame[p_column].dropna(), dtype=float)
    values = np.sort(values[(values > 0) & (values <= 1)])
    expected = -np.log10((np.arange(1, len(values) + 1) - 0.5) / len(values))
    observed = -np.log10(values)
    maximum = max(float(expected.max(initial=0)), float(observed.max(initial=0)), 1.0)
    fig, axis = plt.subplots(figsize=(4.8, 4.6))
    axis.scatter(expected, observed, s=5, color="#2f80c1", alpha=0.72, linewidths=0, rasterized=True)
    axis.plot([0, maximum], [0, maximum], linestyle="--", color="#777777", linewidth=0.9)
    axis.set_xlim(0, maximum)
    axis.set_ylim(0, maximum)
    axis.set_xlabel(r"Expected $-\log_{10}(P)$")
    axis.set_ylabel(r"Observed $-\log_{10}(P)$")
    axis.set_title(f"{title}, $\\lambda_{{GC}}$={lambda_gc(values):.3f}")
    fig.tight_layout()
    for suffix in (".pdf", ".png"):
        fig.savefig(str(output_prefix) + suffix, dpi=360, bbox_inches="tight")
    plt.close(fig)


def load_logistic(path):
    frame = pd.read_csv(path, sep=r"\s+")
    frame = frame[frame["TEST"].astype(str).eq("ADD")].copy()
    frame["P"] = pd.to_numeric(frame["P"], errors="coerce")
    frame["BP"] = pd.to_numeric(frame["BP"], errors="coerce")
    return frame


def write_report(outdir, fisher, logistic, sample_audit):
    summary = pd.DataFrame(
        [
            {
                "method": "Fisher_exact_ALT_carrier",
                "tested_SVs": len(fisher),
                "lambdaGC_descriptive": lambda_gc(fisher["Fisher_P"]),
                "minimum_P": fisher["Fisher_P"].min(),
                "sample_definition": sample_audit,
            },
            {
                "method": "PLINK_additive_logistic_PC1_PC5",
                "tested_SVs": len(logistic),
                "lambdaGC": lambda_gc(logistic["P"]),
                "minimum_P": logistic["P"].min(),
                "sample_definition": sample_audit,
            },
        ]
    )
    summary.to_csv(outdir / "results" / "SV_fisher_logistic_summary.tsv", sep="\t", index=False)
    fisher.nsmallest(20, "Fisher_P").to_csv(outdir / "results" / "SV_fisher_top20.tsv", sep="\t", index=False)
    logistic.nsmallest(20, "P").to_csv(outdir / "results" / "SV_logistic_top20.tsv", sep="\t", index=False)
    with open(outdir / "results" / "SV_fisher_logistic_report.md", "w", encoding="utf-8") as handle:
        handle.write("# Supplementary SV association screen\n\n")
        handle.write("The analysis reproduced the SV GWAS QC: GQ >=20 genotype filtering, marker missingness <=20% in all 144 samples, then marker missingness <=20% and MAF >=0.05 after retaining the common 137-sample GWAS cohort.\n\n")
        handle.write("Fisher's exact test was two-sided and defined event presence directly from the raw VCF: any non-reference ALT genotype (0/1, 1/0, or 1/1) was an ALT carrier. This avoids interpreting PLINK's potentially reordered A1 allele as structural-event presence.\n\n")
        handle.write("PLINK additive logistic regression used PC1-PC5 as fixed covariates and no GRM, so it is a sensitivity analysis rather than the primary mixed-model GWAS.\n\n")
        handle.write(markdown_table(summary))
        handle.write("\n\nBoth tests are exploratory under the strong lineage-phenotype structure; unadjusted Fisher P values must not be interpreted as causal evidence.\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True)
    parser.add_argument("--all-samples", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--fisher-output", required=True)
    parser.add_argument("--logistic", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-audit", required=True)
    parser.add_argument("--maf", type=float, required=True)
    parser.add_argument("--max-initial-missing", type=float, required=True)
    parser.add_argument("--max-final-missing", type=float, required=True)
    parser.add_argument("--expected-records", type=int, required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    (outdir / "results").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    fisher = run_fisher(
        args.query,
        args.all_samples,
        args.phenotype,
        args.fisher_output,
        args.maf,
        args.max_initial_missing,
        args.max_final_missing,
    )
    if len(fisher) != args.expected_records:
        raise ValueError(
            f"Two-stage VCF QC retained {len(fisher)} records; expected {args.expected_records} PLINK markers"
        )
    fisher.to_csv(args.fisher_output, sep="\t", index=False)
    logistic = load_logistic(args.logistic)
    manhattan(fisher, "Fisher_P", outdir / "figures" / "SV_Fisher_manhattan", "SV Fisher exact test")
    qq(fisher, "Fisher_P", outdir / "figures" / "SV_Fisher_QQ", "SV Fisher exact test")
    manhattan(logistic.rename(columns={"CHR": "CHROM", "BP": "POS"}), "P", outdir / "figures" / "SV_logistic_manhattan", "SV logistic regression (PC1-PC5)")
    qq(logistic, "P", outdir / "figures" / "SV_logistic_QQ", "SV logistic regression (PC1-PC5)")
    write_report(outdir, fisher, logistic, args.sample_audit)
    print(f"Fisher tested {len(fisher)} ALT-carrier SVs; logistic tested {len(logistic)} SVs")


if __name__ == "__main__":
    main()
