#!/usr/bin/env python3
"""Audit GWAS calibration without changing association results."""

from __future__ import annotations

import argparse
import glob
import math
import statistics
from pathlib import Path

from scipy.stats import chi2, mannwhitneyu, pointbiserialr


CHI2_MEDIAN_DF1 = chi2.ppf(0.5, 1)


def valid_p(value: str) -> float | None:
    try:
        p = float(value)
    except (TypeError, ValueError):
        return None
    return p if math.isfinite(p) and 0.0 < p <= 1.0 else None


def metrics(pvalues: list[float]) -> dict[str, float | int]:
    if not pvalues:
        return {"n": 0, "lambda_gc": math.nan, "median_p": math.nan, "min_p": math.nan}
    chisq = [chi2.isf(p, 1) for p in pvalues]
    return {
        "n": len(pvalues),
        "lambda_gc": statistics.median(chisq) / CHI2_MEDIAN_DF1,
        "median_p": statistics.median(pvalues),
        "min_p": min(pvalues),
    }


def read_gmmat(path: Path) -> tuple[list[float], dict[str, list[float]]]:
    all_p: list[float] = []
    bins = {
        "MAF_0.05_0.10": [],
        "MAF_0.10_0.20": [],
        "MAF_0.20_0.30": [],
        "MAF_0.30_0.50": [],
    }
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        p_index = header.index("PVAL")
        af_index = header.index("AF")
        for line in handle:
            fields = line.split()
            p = valid_p(fields[p_index])
            if p is None:
                continue
            try:
                af = float(fields[af_index])
            except ValueError:
                continue
            maf = min(af, 1.0 - af)
            all_p.append(p)
            if 0.05 <= maf < 0.10:
                bins["MAF_0.05_0.10"].append(p)
            elif maf < 0.20:
                bins["MAF_0.10_0.20"].append(p)
            elif maf < 0.30:
                bins["MAF_0.20_0.30"].append(p)
            elif maf <= 0.50:
                bins["MAF_0.30_0.50"].append(p)
    return all_p, bins


def read_plink(paths: list[Path]) -> list[float]:
    pvalues: list[float] = []
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            header = handle.readline().split()
            p_index = header.index("P")
            for line in handle:
                fields = line.split()
                if len(fields) <= p_index:
                    continue
                p = valid_p(fields[p_index])
                if p is not None:
                    pvalues.append(p)
    return pvalues


def write_calibration(root: Path, outdir: Path) -> None:
    rows: list[tuple[str, str, str, dict[str, float | int]]] = []
    for label in ("SNP", "INDEL"):
        gmmat = root / "02_snp_indel_gwas" / "gmmat" / f"{label}.coil_ge1.pc5.gmmat.score.tsv"
        pvalues, maf_bins = read_gmmat(gmmat)
        rows.append((label, "GMMAT_GRM_PC5", "all", metrics(pvalues)))
        rows.extend((label, "GMMAT_GRM_PC5", name, metrics(values)) for name, values in maf_bins.items())
        logistic_paths = [Path(p) for p in sorted(glob.glob(str(root / "02_snp_indel_gwas" / "logistic" / f"{label}.chr*.coil_ge1.pc5.assoc.logistic")))]
        rows.append((label, "PLINK_LOGISTIC_PC5", "all", metrics(read_plink(logistic_paths))))

    gmmat = root / "03_sv_gwas" / "gmmat" / "SV.coil_ge1.pc5.gmmat.score.tsv"
    pvalues, maf_bins = read_gmmat(gmmat)
    rows.append(("SV", "GMMAT_GRM_PC5", "all", metrics(pvalues)))
    rows.extend(("SV", "GMMAT_GRM_PC5", name, metrics(values)) for name, values in maf_bins.items())
    logistic = root / "03_sv_gwas" / "logistic" / "SV.coil_ge1.pc5.assoc.logistic"
    rows.append(("SV", "PLINK_LOGISTIC_PC5", "all", metrics(read_plink([logistic]))))

    out_path = outdir / "model_calibration.tsv"
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write("variant_class\tmodel\tstratum\tn_valid\tlambda_gc\tmedian_p\tmin_p\n")
        for label, model, stratum, result in rows:
            handle.write(
                f"{label}\t{model}\t{stratum}\t{result['n']}\t"
                f"{result['lambda_gc']:.6g}\t{result['median_p']:.6g}\t{result['min_p']:.6g}\n"
            )


def write_pc_audit(root: Path, outdir: Path) -> None:
    model = root / "01_structure" / "model" / "SNP.model.tsv"
    with model.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    y_index = header.index("y")
    y = [int(row[y_index]) for row in rows]

    out_path = outdir / "phenotype_pc_association.tsv"
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write("PC\tcase_mean\tcontrol_mean\tpoint_biserial_r\tpoint_biserial_p\tmann_whitney_p\n")
        for pc in ("PC1", "PC2", "PC3", "PC4", "PC5"):
            index = header.index(pc)
            values = [float(row[index]) for row in rows]
            cases = [value for value, phenotype in zip(values, y) if phenotype == 1]
            controls = [value for value, phenotype in zip(values, y) if phenotype == 0]
            correlation = pointbiserialr(y, values)
            mw = mannwhitneyu(cases, controls, alternative="two-sided")
            handle.write(
                f"{pc}\t{statistics.mean(cases):.8g}\t{statistics.mean(controls):.8g}\t"
                f"{correlation.statistic:.8g}\t{correlation.pvalue:.8g}\t{mw.pvalue:.8g}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_calibration(args.root, args.outdir)
    write_pc_audit(args.root, args.outdir)


if __name__ == "__main__":
    main()
