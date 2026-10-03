#!/usr/bin/env python3
"""Population-aware diagnostics for the Chr23997 intron-2 core-loss state.

This is a candidate validation analysis, not a genome-wide scan.  It reports
the naive individual-level Fisher test alongside Firth penalized logistic
models with SNP PCs and a section-stratified species-label permutation.  The
latter prevents repeated accessions of one species from being treated as
independent evolutionary observations.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.special import expit
from scipy.stats import fisher_exact, norm


def dosage(gt: str) -> int:
    return {"0/0": 0, "0/1": 1, "1/1": 2}[gt]


def is_reseq_breakpoint_cluster(row, tolerance: int = 60) -> bool:
    """Keep reference calls and the recurrent 265-bp resequencing breakpoint.

    Other large deletions that merely span the intron core are alternate
    structural haplotypes, not observations of this specific allele, and are
    excluded from a biallelic candidate test.
    """
    if row["GT"] == "0/0":
        return True
    match = re.match(r"(\d+)-(\d+)\(\d+\)", row["modal_del_breakpoint_1based"])
    if not match:
        return False
    left, right = map(int, match.groups())
    return abs(left - 88980832) <= tolerance and abs(right - 88981096) <= tolerance


def fisher_from_rows(rows, phenotypes=None):
    """Fisher test for binary DEL presence, optionally with replacement labels."""
    case_del = case_ref = control_del = control_ref = 0
    for row in rows:
        phenotype = row["phenotype"] if phenotypes is None else phenotypes[row["latin_name"]]
        present = dosage(row["GT"]) > 0
        if phenotype == 1 and present:
            case_del += 1
        elif phenotype == 1:
            case_ref += 1
        elif present:
            control_del += 1
        else:
            control_ref += 1
    table = [[case_del, case_ref], [control_del, control_ref]]
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    return table, odds_ratio, p_value


def _penalized_loglik(X, y, beta):
    eta = np.clip(X @ beta, -30, 30)
    prob = expit(eta)
    loglik = float(np.sum(y * eta - np.logaddexp(0, eta)))
    weights = np.clip(prob * (1 - prob), 1e-10, None)
    sign, logdet = np.linalg.slogdet(X.T @ (weights[:, None] * X))
    return loglik + 0.5 * logdet if sign > 0 else -np.inf


def firth_logistic(X, y, max_iter: int = 200, tolerance: float = 1e-8):
    """Bias-reduced Firth logistic fit using the Jeffreys-prior adjusted score."""
    beta = np.zeros(X.shape[1])
    current = _penalized_loglik(X, y, beta)
    for _ in range(max_iter):
        eta = np.clip(X @ beta, -30, 30)
        prob = expit(eta)
        weights = np.clip(prob * (1 - prob), 1e-10, None)
        information = X.T @ (weights[:, None] * X)
        inverse = np.linalg.pinv(information, rcond=1e-10)
        leverage = weights * np.einsum("ij,jk,ik->i", X, inverse, X)
        score = X.T @ (y - prob + leverage * (0.5 - prob))
        step = inverse @ score
        scale = 1.0
        while scale > 1e-6:
            candidate = beta + scale * step
            value = _penalized_loglik(X, y, candidate)
            if value >= current:
                beta = candidate
                current = value
                break
            scale /= 2.0
        if scale <= 1e-6 or np.max(np.abs(scale * step)) < tolerance:
            break
    eta = np.clip(X @ beta, -30, 30)
    prob = expit(eta)
    weights = np.clip(prob * (1 - prob), 1e-10, None)
    covariance = np.linalg.pinv(X.T @ (weights[:, None] * X), rcond=1e-10)
    return beta, covariance


def standardize(matrix):
    values = matrix.astype(float).copy()
    for index in range(values.shape[1]):
        mean = values[:, index].mean()
        std = values[:, index].std()
        values[:, index] = 0 if std == 0 else (values[:, index] - mean) / std
    return values


def load_rows(calls_path: str, metadata_path: str, eigenvec_path: str):
    metadata = {}
    with open(metadata_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metadata[row["sample_id"]] = row
    pcs = {}
    with open(eigenvec_path) as handle:
        for line in handle:
            fields = line.split()
            pcs[fields[0]] = [float(value) for value in fields[2:7]]
    rows = []
    with open(calls_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["GT"] == "./." or row["sample"] not in pcs:
                continue
            meta = metadata.get(row["sample"], {})
            rows.append({
                **row,
                # The 143-sample manifest is the phenotype source used for the
                # current GWAS.  The annotation table is used only for the
                # species/section-aware sensitivity analysis.
                "phenotype": int(float(row["phenotype"])),
                "latin_name": meta.get("latin_name") or f"sample_{row['sample']}",
                "section": meta.get("section") or "unassigned",
                "pcs": pcs[row["sample"]],
            })
    return rows


def section_permutation(rows, n_perm: int, seed: int):
    """Permute species phenotypes only within taxonomic sections."""
    species_labels = {}
    species_sections = {}
    for row in rows:
        species_labels.setdefault(row["latin_name"], row["phenotype"])
        species_sections.setdefault(row["latin_name"], row["section"])
    for species, label in species_labels.items():
        assert all(row["phenotype"] == label for row in rows if row["latin_name"] == species)
    observed_table, observed_or, observed_p = fisher_from_rows(rows)
    observed_stat = abs(math.log(observed_or)) if observed_or > 0 and math.isfinite(observed_or) else math.inf
    groups = defaultdict(list)
    for species, section in species_sections.items():
        groups[section].append(species)
    rng = np.random.default_rng(seed)
    exceed = valid = 0
    for _ in range(n_perm):
        permuted = dict(species_labels)
        for members in groups.values():
            labels = [species_labels[species] for species in members]
            rng.shuffle(labels)
            permuted.update(dict(zip(members, labels)))
        _, value, _ = fisher_from_rows(rows, permuted)
        if value <= 0 or not math.isfinite(value):
            continue
        valid += 1
        if abs(math.log(value)) >= observed_stat - 1e-12:
            exceed += 1
    empirical = (exceed + 1) / (valid + 1) if valid else math.nan
    return observed_table, observed_or, observed_p, empirical, valid, len(species_labels)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--eigenvec", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260716)
    parser.add_argument("--allele-mode", choices=["core_loss", "reseq_breakpoint_cluster"], default="core_loss")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = load_rows(args.calls, args.metadata, args.eigenvec)
    if args.allele_mode == "reseq_breakpoint_cluster":
        rows = [row for row in rows if is_reseq_breakpoint_cluster(row)]
    table, odds_ratio, fisher_p = fisher_from_rows(rows)
    y = np.array([row["phenotype"] for row in rows], dtype=float)
    dosages = np.array([dosage(row["GT"]) for row in rows], dtype=float)
    pc_matrix = np.array([row["pcs"] for row in rows])
    models = [("Firth_PC2", 2), ("Firth_PC5", 5)]
    output = []
    output.append({"model": "Fisher_unadjusted", "n": len(rows), "beta_dosage": ".",
                   "odds_ratio_dosage": odds_ratio, "wald_p": fisher_p,
                   "notes": f"case_DEL={table[0][0]};case_REF={table[0][1]};control_DEL={table[1][0]};control_REF={table[1][1]}"})
    for name, n_pc in models:
        design = np.column_stack([np.ones(len(rows)), standardize(dosages[:, None]), standardize(pc_matrix[:, :n_pc])])
        beta, covariance = firth_logistic(design, y)
        se = float(np.sqrt(max(covariance[1, 1], 0)))
        z = beta[1] / se if se else math.nan
        p_value = 2 * norm.sf(abs(z)) if math.isfinite(z) else math.nan
        output.append({"model": name, "n": len(rows), "beta_dosage": beta[1],
                       "odds_ratio_dosage": math.exp(beta[1]), "wald_p": p_value,
                       "notes": f"Firth_penalized_logistic;PCs=1-{n_pc}"})

    section_rows = [row for row in rows if row["section"] != "unassigned"]
    species_rows = []
    for species in sorted({row["latin_name"] for row in section_rows}):
        subset = [row for row in section_rows if row["latin_name"] == species]
        values = [dosage(row["GT"]) for row in subset]
        modal = Counter(values).most_common(1)[0][0]
        species_rows.append({"latin_name": species, "section": subset[0]["section"], "phenotype": subset[0]["phenotype"],
                             "n_called_samples": len(subset), "mean_DEL_dosage": np.mean(values), "modal_GT": ["0/0", "0/1", "1/1"][modal]})
    prefix = "Chr23997" if args.allele_mode == "core_loss" else "Chr23997_reseq_breakpoint_cluster"
    with (outdir / f"{prefix}_species_collapsed.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(species_rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(species_rows)

    p_table, p_or, p_value, empirical_p, valid, n_species = section_permutation(section_rows, args.permutations, args.seed)
    output.append({"model": "Section_stratified_species_permutation", "n": len(section_rows), "beta_dosage": ".",
                   "odds_ratio_dosage": p_or, "wald_p": empirical_p,
                   "notes": f"species={n_species};permutations={valid};nominal_fisher_p={p_value}"})
    with (outdir / f"{prefix}_candidate_association_models.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(output)
    with (outdir / f"{prefix}_candidate_diagnostic.md").open("w") as handle:
        handle.write("# Chr23997 intron-2 DEL candidate diagnostics\n\n")
        handle.write(f"Allele mode: `{args.allele_mode}`.\n\n")
        handle.write(f"Called samples for individual-level models: {len(rows)}. ")
        handle.write(f"Samples with reliable species/section metadata for the permutation: {len(section_rows)}; species: {n_species}.\n\n")
        handle.write("The Fisher result treats accessions as independent and is descriptive only. ")
        handle.write("Firth-PC models account for SNP PCA; the section-stratified permutation treats species as the exchangeable unit within section.\n")


if __name__ == "__main__":
    main()
