#!/usr/bin/env python3
"""Write an evidence-based recommendation from the completed A/B/C audit."""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2


CHI2_MEDIAN = chi2.ppf(0.5, 1)


def gmmat_paths(root: Path, model: str, label: str) -> list[Path]:
    sens = root / "05_model_sensitivity_20260715"
    if model == "A":
        if label in ("SNP", "INDEL"):
            return [root / "02_snp_indel_gwas/gmmat" / f"{label}.coil_ge1.pc5.gmmat.score.tsv"]
        return [root / "03_sv_gwas/gmmat/SV.coil_ge1.pc5.gmmat.score.tsv"]
    if label in ("SNP", "INDEL"):
        return sorted(Path(p) for p in glob.glob(str(sens / "model_C_gmmat_grm_pc2" / label / "*.gmmat.score.tsv")))
    return [sens / "model_C_gmmat_grm_pc2/SV/SV.grm_pc2.gmmat.score.tsv"]


def common_metrics(paths: list[Path], cutoff: float) -> tuple[int, float, float]:
    pvalues = []
    for path in paths:
        for chunk in pd.read_csv(path, sep=r"\s+", usecols=["AF", "PVAL"], chunksize=200_000):
            af = pd.to_numeric(chunk.AF, errors="coerce").to_numpy()
            p = pd.to_numeric(chunk.PVAL, errors="coerce").to_numpy()
            maf = np.minimum(af, 1.0 - af)
            keep = np.isfinite(p) & (p > 0) & (p <= 1) & (maf >= cutoff)
            pvalues.append(p[keep])
    p = np.concatenate(pvalues)
    lam = float(np.median(chi2.isf(p, 1)) / CHI2_MEDIAN)
    return len(p), lam, float(p.min())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    comparison = pd.read_csv(args.outdir / "ABC_model_comparison.tsv", sep="\t")
    maf_bins = pd.read_csv(args.outdir / "ABC_lambda_by_MAF.tsv", sep="\t")
    diagnostics = pd.read_csv(
        args.root / "05_model_sensitivity_20260715/diagnostics/structure_diagnostics.tsv",
        sep="\t",
    ).set_index("metric").value

    common_rows = []
    for model in ("A", "C"):
        for label in ("SNP", "INDEL", "SV"):
            for cutoff in (0.10, 0.20):
                tested, lam, minimum = common_metrics(gmmat_paths(args.root, model, label), cutoff)
                common_rows.append({
                    "Model": model, "Variant_class": label, "MAF_cutoff": cutoff,
                    "Tested": tested, "lambdaGC": lam, "Best_P": minimum,
                })
    common = pd.DataFrame(common_rows)
    common.to_csv(args.outdir / "AC_common_variant_calibration.tsv", sep="\t", index=False)

    low_maf = maf_bins[(maf_bins.Model == "A") & (maf_bins.MAF_bin == "0.05-0.10")].set_index("Variant_class")
    full_a = comparison[comparison.Model == "A"].set_index("Variant_class")
    percentages = {label: 100.0 * low_maf.loc[label, "Tested"] / full_a.loc[label, "Tested"] for label in ("SNP", "INDEL", "SV")}
    common10 = common[common.MAF_cutoff == 0.10].set_index(["Model", "Variant_class"])

    best_ids = comparison.pivot(index="Variant_class", columns="Model", values="Best_ID")
    stable_ac = [label for label in ("SNP", "INDEL", "SV") if best_ids.loc[label, "A"] == best_ids.loc[label, "C"]]

    with (args.outdir / "FINAL_MODEL_RECOMMENDATION.md").open("w", encoding="utf-8") as handle:
        handle.write("# 144-sample coil GWAS model sensitivity conclusion\n\n")
        handle.write("## Main diagnosis\n\n")
        handle.write(
            "The full-set QQ deflation is not explained primarily by GMMAT or by PC4-PC5. "
            "Variants with MAF 0.05-0.10 account for "
            f"{percentages['SNP']:.1f}% of SNPs, {percentages['INDEL']:.1f}% of INDELs, and "
            f"{percentages['SV']:.1f}% of SVs, and this stratum has the strongest deflation.\n\n"
        )
        handle.write(
            f"PC2 is strongly aligned with the phenotype (point-biserial r={float(diagnostics['PC2_point_biserial_r']):.3f}, "
            f"P={float(diagnostics['PC2_point_biserial_p']):.3g}); therefore a no-GRM model is at substantial risk of "
            "confounding by taxonomy/phylogeny.\n\n"
        )
        handle.write("## Model comparison\n\n")
        handle.write("- Model A (GMMAT, PC1-PC5 + GRM) remains the conservative baseline.\n")
        handle.write("- Model B (PLINK, PC1-PC3, no GRM) is not suitable as the primary model: calibration is inconsistent across marker classes and many tests return no valid P value.\n")
        handle.write("- Model C (GMMAT, PC1-PC2 + GRM) does not improve the unstratified QQ plot, but is best calibrated for common SNPs and INDELs after MAF >= 0.10 filtering.\n\n")
        handle.write("At MAF >= 0.10, Model C lambdaGC values are "
                     f"{common10.loc[('C', 'SNP'), 'lambdaGC']:.3f} (SNP), "
                     f"{common10.loc[('C', 'INDEL'), 'lambdaGC']:.3f} (INDEL), and "
                     f"{common10.loc[('C', 'SV'), 'lambdaGC']:.3f} (SV).\n\n")
        handle.write(f"The genome-wide best marker is identical between Models A and C for {len(stable_ac)}/3 marker classes ({', '.join(stable_ac)}), supporting peak stability to PC count.\n\n")
        handle.write("## Recommended reporting strategy\n\n")
        handle.write("1. Use Model C with MAF >= 0.10 as the primary common-variant SNP/INDEL analysis.\n")
        handle.write("2. Require direction and peak concordance with Model A (PC1-PC5 + GRM) as a sensitivity criterion.\n")
        handle.write("3. Keep Model A as the primary confounding-controlled SV analysis; describe its residual deflation and limited power.\n")
        handle.write("4. Do not promote Model B to primary solely because the SV lambdaGC is close to 1.\n")
        handle.write("5. Treat all associations as exploratory because the panel spans many Medicago taxa and no marker reaches the experiment-wide threshold.\n")


if __name__ == "__main__":
    main()
