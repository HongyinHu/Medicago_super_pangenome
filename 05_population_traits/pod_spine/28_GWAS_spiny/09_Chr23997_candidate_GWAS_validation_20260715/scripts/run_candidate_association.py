#!/usr/bin/env python3
"""Run direct candidate-SV association analyses without imputing missing genotypes."""

import argparse
import csv
import math
import warnings
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def contingency_table(rows):
    """Return [[DEL spiny, PRESENT spiny], [DEL spineless, PRESENT spineless]]."""
    table = [[0, 0], [0, 0]]
    for row in rows:
        if row["genotype"] == "NA":
            continue
        phenotype = int(row["phenotype"])
        if phenotype not in (0, 1):
            raise ValueError("phenotype must be 0 or 1")
        genotype_index = 0 if row["genotype"] == "1" else 1
        table[0 if phenotype == 1 else 1][genotype_index] += 1
    return table


def odds_ratio(table):
    """DEL odds of spiny phenotype relative to PRESENT, with no continuity correction."""
    numerator = table[0][0] * table[1][1]
    denominator = table[0][1] * table[1][0]
    if denominator == 0:
        return math.inf if numerator else math.nan
    return numerator / float(denominator)


def sample_order(path):
    values = []
    with open(path, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("invalid keep row")
            values.append(fields[1])
    return values


def pc_map(path, order):
    values = []
    with open(path, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split()
            if not fields:
                continue
            if len(fields) < 6:
                raise ValueError("PC file needs intercept plus PC1-PC5")
            values.append([float(value) for value in fields[1:6]])
    if len(values) != len(order):
        raise ValueError("PC row count does not match analysis sample order")
    return {sample: vector for sample, vector in zip(order, values)}


def exp_or_na(value):
    if not np.isfinite(value):
        return "NA"
    try:
        result = math.exp(value)
    except OverflowError:
        return "Inf"
    return "{:.10g}".format(result) if np.isfinite(result) else "Inf"


def logistic_irls(x, y, max_iter=200, tolerance=1e-9):
    """Fit a logistic model by Newton-Raphson without an external statistics package."""
    beta = np.zeros(x.shape[1], dtype=float)
    for _ in range(max_iter):
        linear = np.clip(x @ beta, -35.0, 35.0)
        probability = 1.0 / (1.0 + np.exp(-linear))
        weight = np.maximum(probability * (1.0 - probability), 1e-12)
        information = x.T @ (x * weight[:, None])
        score = x.T @ (y - probability)
        try:
            step = np.linalg.solve(information, score)
        except np.linalg.LinAlgError:
            raise ValueError("singular logistic information matrix")
        beta_next = beta + step
        if np.max(np.abs(beta_next - beta)) < tolerance:
            covariance = np.linalg.inv(information)
            return beta_next, covariance
        beta = beta_next
        if np.max(np.abs(beta)) > 30.0:
            raise ValueError("complete or quasi-complete separation")
    raise ValueError("logistic model did not converge")


def logistic_result(rows, pcs=None):
    called = [row for row in rows if row["genotype"] != "NA"]
    y = np.asarray([int(row["phenotype"]) for row in called], dtype=float)
    geno = np.asarray([float(row["genotype"]) for row in called], dtype=float)
    columns = [np.ones(len(called)), geno]
    if pcs is not None:
        columns.extend(np.asarray([pcs[row["sample"]] for row in called], dtype=float).T)
    x = np.column_stack(columns)
    if len(np.unique(y)) < 2 or len(np.unique(geno)) < 2:
        return {"status": "not_estimable: no phenotype or genotype variation"}
    try:
        fit_beta, covariance = logistic_irls(x, y)
        beta = float(fit_beta[1])
        se = float(math.sqrt(covariance[1, 1]))
        z_score = beta / se
        pvalue = float(math.erfc(abs(z_score) / math.sqrt(2.0)))
        low = beta - 1.96 * se
        high = beta + 1.96 * se
        return {
            "status": "ok",
            "beta": beta,
            "se": se,
            "odds_ratio": math.exp(beta),
            "ci_low": math.exp(low),
            "ci_high": math.exp(high),
            "pvalue": pvalue,
        }
    except (np.linalg.LinAlgError, ValueError) as exc:
        return {"status": "not_estimable: {}".format(exc)}


def format_number(value):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "NA"
    return "{:.10g}".format(value)


def summarize_callset(callset, rows, pcs):
    table = contingency_table(rows)
    called = [row for row in rows if row["genotype"] != "NA"]
    common = {
        "callset": callset,
        "n_total": len(rows),
        "n_called": len(called),
        "n_missing": len(rows) - len(called),
        "del_spiny": table[0][0],
        "present_spiny": table[0][1],
        "del_spineless": table[1][0],
        "present_spineless": table[1][1],
        "n_del": table[0][0] + table[1][0],
        "n_present": table[0][1] + table[1][1],
    }
    fisher_or, fisher_p = fisher_exact(table, alternative="two-sided")
    fisher_row = dict(common)
    fisher_row.update({
        "model": "Fisher_exact_two_sided",
        "beta": format_number(math.log(fisher_or)) if fisher_or > 0 and np.isfinite(fisher_or) else "NA",
        "se": "NA",
        "odds_ratio": format_number(fisher_or),
        "ci_low": "NA",
        "ci_high": "NA",
        "pvalue": format_number(fisher_p),
        "status": "ok",
    })
    output = [fisher_row]
    for model, pc_values in (("Logistic_unadjusted", None), ("Logistic_PC1_PC5", pcs)):
        result = logistic_result(rows, pc_values)
        row = dict(common)
        row.update({
            "model": model,
            "beta": format_number(result.get("beta")),
            "se": format_number(result.get("se")),
            "odds_ratio": format_number(result.get("odds_ratio")),
            "ci_low": format_number(result.get("ci_low")),
            "ci_high": format_number(result.get("ci_high")),
            "pvalue": format_number(result.get("pvalue")),
            "status": result["status"],
        })
        output.append(row)
    return output, table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", required=True)
    parser.add_argument("--sensitivity", required=True)
    parser.add_argument("--sample-order", required=True)
    parser.add_argument("--pc-file", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    order = sample_order(args.sample_order)
    pcs = pc_map(args.pc_file, order)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_results = []
    contingency_rows = []
    for callset, path in (("strict", args.strict), ("sensitivity", args.sensitivity)):
        rows = read_tsv(path)
        actual_order = [row["sample"] for row in rows]
        if actual_order != order:
            raise ValueError("candidate input sample order differs from analysis sample order")
        model_rows, table = summarize_callset(callset, rows, pcs)
        all_results.extend(model_rows)
        labels = (("spiny", "DEL", table[0][0]), ("spiny", "PRESENT", table[0][1]),
                  ("spineless", "DEL", table[1][0]), ("spineless", "PRESENT", table[1][1]))
        contingency_rows.extend({"callset": callset, "phenotype_group": group, "allele_state": allele, "n": count}
                                for group, allele, count in labels)

    fields = ["callset", "model", "n_total", "n_called", "n_missing", "n_del", "n_present",
              "del_spiny", "present_spiny", "del_spineless", "present_spineless", "beta", "se",
              "odds_ratio", "ci_low", "ci_high", "pvalue", "status"]
    write_tsv(outdir / "Chr23997.candidate_association.all_models.tsv", fields, all_results)
    write_tsv(outdir / "Chr23997.candidate_association.fisher_contingency.tsv",
              ["callset", "phenotype_group", "allele_state", "n"], contingency_rows)


if __name__ == "__main__":
    main()
