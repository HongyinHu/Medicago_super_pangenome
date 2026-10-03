#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2


P_COLUMNS = ("p_wald", "p_lrt", "p_score")
ALPHAS = (0.5, 0.1, 0.05, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 5e-8)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def lambda_gc(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0) & (values < 1)]
    if values.size == 0:
        return np.nan
    return float(np.median(chi2.isf(values, 1)) / chi2.ppf(0.5, 1))


def association_diagnostics(label, path):
    columns = list(P_COLUMNS) + ["l_remle", "l_mle"]
    frame = pd.read_csv(path, sep="\t", usecols=columns, low_memory=False)
    rows = []
    for test in P_COLUMNS:
        values = pd.to_numeric(frame[test], errors="coerce").to_numpy(dtype=float)
        values = values[np.isfinite(values) & (values > 0) & (values <= 1)]
        row = {
            "dataset": label,
            "test": test,
            "n": values.size,
            "lambda_gc": lambda_gc(values),
            "min_p": float(values.min()),
            "median_p": float(np.median(values)),
            "q01_p": float(np.quantile(values, 0.01)),
            "q05_p": float(np.quantile(values, 0.05)),
            "q10_p": float(np.quantile(values, 0.10)),
        }
        for alpha in ALPHAS:
            key = f"prop_p_le_{alpha:g}".replace("-", "m").replace(".", "p")
            row[key] = float(np.mean(values <= alpha))
        rows.append(row)

    boundary = {
        "dataset": label,
        "n_rows": len(frame),
        "l_remle_min": float(pd.to_numeric(frame["l_remle"], errors="coerce").min()),
        "l_remle_max": float(pd.to_numeric(frame["l_remle"], errors="coerce").max()),
        "l_remle_at_100000": int(np.isclose(frame["l_remle"], 100000).sum()),
        "l_mle_min": float(pd.to_numeric(frame["l_mle"], errors="coerce").min()),
        "l_mle_max": float(pd.to_numeric(frame["l_mle"], errors="coerce").max()),
        "l_mle_at_100000": int(np.isclose(frame["l_mle"], 100000).sum()),
    }
    return rows, boundary


def linear_r2(y, design):
    y = np.asarray(y, dtype=float)
    design = np.asarray(design, dtype=float)
    if design.ndim == 1:
        design = design[:, None]
    x = np.column_stack([np.ones(len(y)), design])
    fitted = x @ np.linalg.lstsq(x, y, rcond=None)[0]
    total = np.sum((y - y.mean()) ** 2)
    residual = np.sum((y - fitted) ** 2)
    return float(1 - residual / total) if total > 0 else np.nan


def kinship_diagnostics(label, kinship_path, phenotype_path, covariate_path):
    kinship = np.loadtxt(kinship_path)
    phenotype = np.loadtxt(phenotype_path)
    covariates = np.loadtxt(covariate_path)
    if covariates.ndim == 1:
        covariates = covariates[:, None]
    if kinship.shape[0] != kinship.shape[1] or kinship.shape[0] != len(phenotype):
        raise RuntimeError(
            f"Dimension mismatch for {label}: K={kinship.shape}, y={phenotype.shape}"
        )

    eigvals, eigvecs = np.linalg.eigh((kinship + kinship.T) / 2)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pcs = covariates[:, 1:] if covariates.shape[1] > 1 else np.empty((len(phenotype), 0))

    kinship_row = {
        "dataset": label,
        "n": len(phenotype),
        "trait_0": int(np.sum(phenotype == 0)),
        "trait_1": int(np.sum(phenotype == 1)),
        "kinship_diag_mean": float(np.mean(np.diag(kinship))),
        "kinship_diag_min": float(np.min(np.diag(kinship))),
        "kinship_diag_max": float(np.max(np.diag(kinship))),
        "kinship_offdiag_mean": float(
            (kinship.sum() - np.trace(kinship)) / (kinship.size - len(kinship))
        ),
        "kinship_asymmetry_max": float(np.max(np.abs(kinship - kinship.T))),
        "eig_min": float(eigvals.min()),
        "eig_max": float(eigvals.max()),
        "eig_negative_count": int(np.sum(eigvals < -1e-8)),
        "trait_r2_pc1_5": linear_r2(phenotype, pcs) if pcs.shape[1] else np.nan,
    }

    for k in (1, 2, 3, 5, 10, 20):
        use = min(k, eigvecs.shape[1])
        kinship_row[f"trait_r2_topK_eigenvectors_{k}"] = linear_r2(
            phenotype, eigvecs[:, :use]
        )

    corr_rows = []
    for index in range(pcs.shape[1]):
        corr = np.corrcoef(phenotype, pcs[:, index])[0, 1]
        corr_rows.append(
            {"dataset": label, "pc": index + 1, "trait_pearson_r": float(corr)}
        )
    return kinship_row, corr_rows


def main():
    args = parse_args()
    root = args.root.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    associations = {
        "SNP_pc5": root / "02_snp_indel_gwas/gemma/output/SNP.coil.pc5.assoc.txt",
        "INDEL_pc5": root / "02_snp_indel_gwas/gemma/output/INDEL.coil.pc5.assoc.txt",
        "SV_pc5": root / "03_sv_gwas/gemma/output/SV.coil.pc5.gq20.twostage.assoc.txt",
        "SNP_kinship_only": root
        / "02_snp_indel_gwas/kinship_only/gemma/output/SNP.coil.kinship_only.assoc.txt",
        "INDEL_kinship_only": root
        / "02_snp_indel_gwas/kinship_only/gemma/output/INDEL.coil.kinship_only.assoc.txt",
        "SV_kinship_only": root
        / "03_sv_gwas/kinship_only/gemma/output/SV.coil.kinship_only.assoc.txt",
    }

    calibration_rows = []
    boundary_rows = []
    for label, path in associations.items():
        if not path.is_file():
            continue
        rows, boundary = association_diagnostics(label, path)
        calibration_rows.extend(rows)
        boundary_rows.append(boundary)

    kinship_inputs = {
        "SNP_INDEL_125": (
            root / "01_structure/gemma/output/coil_SNP_kinship.cXX.txt",
            root / "01_structure/summary/phenotype.0_1.txt",
            root / "01_structure/summary/covariates.pc5.txt",
        ),
        "SV_119": (
            root / "03_sv_gwas/summary/SV_kinship.reordered.cXX.txt",
            root / "03_sv_gwas/summary/SV.phenotype.0_1.txt",
            root / "03_sv_gwas/summary/SV.covariates.pc5.txt",
        ),
    }
    kinship_rows = []
    correlation_rows = []
    for label, paths in kinship_inputs.items():
        row, correlations = kinship_diagnostics(label, *paths)
        kinship_rows.append(row)
        correlation_rows.extend(correlations)

    calibration = pd.DataFrame(calibration_rows)
    boundaries = pd.DataFrame(boundary_rows)
    kinship = pd.DataFrame(kinship_rows)
    correlations = pd.DataFrame(correlation_rows)
    calibration.to_csv(outdir / "pvalue_calibration.tsv", sep="\t", index=False)
    boundaries.to_csv(outdir / "gemma_boundary_parameters.tsv", sep="\t", index=False)
    kinship.to_csv(outdir / "kinship_diagnostics.tsv", sep="\t", index=False)
    correlations.to_csv(outdir / "pc_trait_correlations.tsv", sep="\t", index=False)

    with (outdir / "diagnosis.txt").open("w", encoding="utf-8") as handle:
        handle.write("Observed diagnostic facts; interpretation requires review.\n")
        handle.write("\nGEMMA boundary parameters:\n")
        handle.write(boundaries.to_string(index=False))
        handle.write("\n\nKinship and phenotype diagnostics:\n")
        handle.write(kinship.to_string(index=False))
        handle.write("\n\nP-value calibration:\n")
        handle.write(calibration.to_string(index=False))
        handle.write("\n")

    print(f"Wrote diagnostics to {outdir}")


if __name__ == "__main__":
    main()
