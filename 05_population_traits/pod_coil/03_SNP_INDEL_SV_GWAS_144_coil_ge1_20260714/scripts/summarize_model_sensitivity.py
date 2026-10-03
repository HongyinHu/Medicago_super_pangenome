#!/usr/bin/env python3
"""Compare Model A/B/C GWAS calibration and peak stability."""

from __future__ import annotations

import argparse
import glob
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2


MODELS = {
    "A": ("GMMAT", "PC1-PC5", "Yes"),
    "B": ("PLINK", "PC1-PC3", "No"),
    "C": ("GMMAT", "PC1-PC2", "Yes"),
}
CLASSES = ("SNP", "INDEL", "SV")
GENES = {
    "Chr14120": (3, 53_196_386, 53_203_316),
    "Chr24229": (4, 91_173_307, 91_178_034),
}
CHI2_MEDIAN = chi2.ppf(0.5, 1)


def choose(columns: list[str], *candidates: str) -> str:
    lookup = {column.upper(): column for column in columns}
    for candidate in candidates:
        if candidate.upper() in lookup:
            return lookup[candidate.upper()]
    raise RuntimeError(f"Missing columns {candidates} from {columns}")


def normalize(frame: pd.DataFrame, model: str, label: str) -> pd.DataFrame:
    chrom = choose(list(frame.columns), "CHR", "CHROM")
    marker = choose(list(frame.columns), "SNP", "ID", "MARKER")
    pos = choose(list(frame.columns), "POS", "BP")
    pvalue = choose(list(frame.columns), "PVAL", "P", "PVALUE")
    keep = [chrom, marker, pos, pvalue]
    if "AF" in frame.columns:
        keep.append("AF")
    frame = frame[keep].rename(columns={chrom: "chrom", marker: "id", pos: "pos", pvalue: "p_value"})
    frame["chrom"] = pd.to_numeric(frame.chrom.astype(str).str.replace("Chr", "", regex=False), errors="coerce")
    frame["pos"] = pd.to_numeric(frame.pos, errors="coerce")
    frame["p_value"] = pd.to_numeric(frame.p_value, errors="coerce")
    frame = frame[
        frame.chrom.between(1, 8) & frame.pos.notna() & frame.p_value.gt(0) & frame.p_value.le(1)
    ].copy()
    frame["chrom"] = frame.chrom.astype(int)
    frame["pos"] = frame.pos.astype(int)
    frame["model"] = model
    frame["class"] = label
    return frame


def read_gmmat(paths: list[Path], model: str, label: str) -> pd.DataFrame:
    return pd.concat(
        [normalize(pd.read_csv(path, sep=r"\s+", low_memory=False), model, label) for path in paths],
        ignore_index=True,
    )


def read_plink(paths: list[Path], model: str, label: str) -> pd.DataFrame:
    frames = []
    for path in paths:
        frame = pd.read_csv(path, sep=r"\s+", low_memory=False, na_values=["NA"])
        if "TEST" in frame.columns:
            frame = frame[frame.TEST == "ADD"]
        frames.append(normalize(frame, model, label))
    return pd.concat(frames, ignore_index=True)


def assoc_paths(root: Path, sens: Path, model: str, label: str) -> tuple[str, list[Path]]:
    if model == "A":
        if label in ("SNP", "INDEL"):
            return "gmmat", [root / "02_snp_indel_gwas/gmmat" / f"{label}.coil_ge1.pc5.gmmat.score.tsv"]
        return "gmmat", [root / "03_sv_gwas/gmmat/SV.coil_ge1.pc5.gmmat.score.tsv"]
    if model == "B":
        if label in ("SNP", "INDEL"):
            paths = sorted(Path(p) for p in glob.glob(str(sens / "model_B_plink_pc3" / label / f"{label}.chr*.pc3.assoc.logistic")))
        else:
            paths = [sens / "model_B_plink_pc3/SV/SV.pc3.assoc.logistic"]
        return "plink", paths
    if label in ("SNP", "INDEL"):
        paths = sorted(Path(p) for p in glob.glob(str(sens / "model_C_gmmat_grm_pc2" / label / f"{label}.chr*.grm_pc2.gmmat.score.tsv")))
    else:
        paths = [sens / "model_C_gmmat_grm_pc2/SV/SV.grm_pc2.gmmat.score.tsv"]
    return "gmmat", paths


def lambda_gc(pvalues: np.ndarray) -> float:
    pvalues = pvalues[np.isfinite(pvalues) & (pvalues > 0) & (pvalues <= 1)]
    return float(np.median(chi2.isf(pvalues, 1)) / CHI2_MEDIAN)


def qq_points(pvalues: np.ndarray, maximum: int = 180_000) -> tuple[np.ndarray, np.ndarray]:
    observed = np.sort(pvalues[np.isfinite(pvalues) & (pvalues > 0) & (pvalues <= 1)])
    expected = (np.arange(1, len(observed) + 1) - 0.5) / len(observed)
    x = -np.log10(expected)
    y = -np.log10(observed)
    if len(x) > maximum:
        tail = min(10_000, len(x))
        bulk = np.linspace(tail, len(x) - 1, maximum - tail, dtype=int)
        index = np.unique(np.concatenate((np.arange(tail), bulk)))
        x, y = x[index], y[index]
    return x, y


def cumulative(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[int, float], dict[int, float]]:
    maxima = frame.groupby("chrom").pos.max().reindex(range(1, 9)).fillna(0)
    offsets: dict[int, float] = {}
    centers: dict[int, float] = {}
    offset = 0.0
    for chrom in range(1, 9):
        offsets[chrom] = offset
        centers[chrom] = offset + maxima.loc[chrom] / 2.0
        offset += maxima.loc[chrom] + 2_000_000
    result = frame.copy()
    result["x"] = result.apply(lambda row: row.pos + offsets[int(row.chrom)], axis=1)
    return result, offsets, centers


def save(fig: plt.Figure, outdir: Path, stem: str) -> None:
    fig.savefig(outdir / f"{stem}.png", dpi=400, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_qq(frames: dict[tuple[str, str], pd.DataFrame], outdir: Path) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(11.2, 10.2))
    for row, label in enumerate(CLASSES):
        row_limit = 0.0
        points = {}
        for column, model in enumerate(MODELS):
            frame = frames[(model, label)]
            x, y = qq_points(frame.p_value.to_numpy())
            points[model] = (x, y)
            row_limit = max(row_limit, float(x.max()), float(y.max()))
        row_limit = math.ceil(row_limit * 10) / 10
        for column, model in enumerate(MODELS):
            ax = axes[row, column]
            x, y = points[model]
            ax.scatter(x, y, s=2.0, color="#4C9ED9", alpha=0.65, rasterized=True)
            ax.plot([0, row_limit], [0, row_limit], linestyle="--", color="#777777", linewidth=0.8)
            ax.set_xlim(0, row_limit)
            ax.set_ylim(0, row_limit)
            ax.set_aspect("equal", adjustable="box")
            ax.set_title(f"{label}, Model {model}, λ={lambda_gc(frames[(model, label)].p_value.to_numpy()):.3f}", fontsize=9)
            if row == 2:
                ax.set_xlabel(r"Expected $-\log_{10}(P)$")
            if column == 0:
                ax.set_ylabel(r"Observed $-\log_{10}(P)$")
            ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("GWAS model calibration: identical phenotype and variant QC", y=1.005)
    fig.tight_layout()
    save(fig, outdir, "ABC_model_comparison_QQ")


def draw_manhattan(frames: dict[tuple[str, str], pd.DataFrame], outdir: Path) -> None:
    palette = ("#75B9E7", "#E83E8C")
    for model in MODELS:
        fig, axes = plt.subplots(3, 1, figsize=(13, 8.5), sharex=False)
        for ax, label in zip(axes, CLASSES):
            frame, _, centers = cumulative(frames[(model, label)])
            for chrom in range(1, 9):
                subset = frame[frame.chrom == chrom]
                ax.scatter(subset.x, -np.log10(subset.p_value), s=1.0, color=palette[(chrom - 1) % 2], alpha=0.72, rasterized=True)
            n = len(frame)
            ax.axhline(-math.log10(1 / n), color="#777777", linestyle="--", linewidth=0.7)
            ax.axhline(-math.log10(0.05 / n), color="#D62728", linestyle=":", linewidth=0.7)
            ax.set_ylabel(f"{label}\n$-\\log_{{10}}(P)$")
            ax.set_xticks([centers[c] for c in range(1, 9)])
            ax.set_xticklabels(range(1, 9))
            ax.spines[["top", "right"]].set_visible(False)
        axes[-1].set_xlabel("Chromosome")
        software, covariates, kinship = MODELS[model]
        fig.suptitle(f"Model {model}: {software}, {covariates}, kinship={kinship}")
        fig.tight_layout()
        save(fig, outdir, f"Model_{model}_SNP_INDEL_SV_Manhattan")


def write_tables(frames: dict[tuple[str, str], pd.DataFrame], outdir: Path, sample_counts: dict[str, int]) -> None:
    rows = []
    strata = []
    candidates = []
    for model in MODELS:
        software, covariates, kinship = MODELS[model]
        for label in CLASSES:
            frame = frames[(model, label)]
            best = frame.loc[frame.p_value.idxmin()]
            rows.append({
                "Model": model, "Software": software, "Covariates": covariates, "Kinship": kinship,
                "Variant_class": label, "N": sample_counts[label], "Tested": len(frame),
                "lambdaGC": lambda_gc(frame.p_value.to_numpy()), "Best_ID": best.id,
                "Best_CHR": int(best.chrom), "Best_POS": int(best.pos), "Best_P": best.p_value,
            })
            if "AF" in frame.columns:
                maf = np.minimum(frame.AF.to_numpy(dtype=float), 1 - frame.AF.to_numpy(dtype=float))
                for name, lower, upper in (
                    ("0.05-0.10", 0.05, 0.10), ("0.10-0.20", 0.10, 0.20),
                    ("0.20-0.30", 0.20, 0.30), ("0.30-0.50", 0.30, 0.5000001),
                ):
                    p = frame.loc[(maf >= lower) & (maf < upper), "p_value"].to_numpy()
                    strata.append({"Model": model, "Variant_class": label, "MAF_bin": name, "Tested": len(p), "lambdaGC": lambda_gc(p) if len(p) else np.nan})
            for gene, (chrom, start, end) in GENES.items():
                for window in (0, 50_000, 500_000):
                    subset = frame[(frame.chrom == chrom) & frame.pos.between(start - window, end + window)]
                    candidates.append({
                        "Model": model, "Variant_class": label, "Gene": gene, "Window_bp": window,
                        "Tested": len(subset), "Best_P": subset.p_value.min() if len(subset) else np.nan,
                        "Best_ID": subset.loc[subset.p_value.idxmin(), "id"] if len(subset) else "NA",
                    })
    pd.DataFrame(rows).to_csv(outdir / "ABC_model_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(strata).to_csv(outdir / "ABC_lambda_by_MAF.tsv", sep="\t", index=False)
    pd.DataFrame(candidates).to_csv(outdir / "ABC_candidate_peak_stability.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    sens = args.root / "05_model_sensitivity_20260715"
    args.outdir.mkdir(parents=True, exist_ok=True)
    frames: dict[tuple[str, str], pd.DataFrame] = {}
    for model in MODELS:
        for label in CLASSES:
            kind, paths = assoc_paths(args.root, sens, model, label)
            if not paths or any(not path.is_file() or path.stat().st_size == 0 for path in paths):
                raise RuntimeError(f"Missing association files for Model {model}, {label}: {paths}")
            frames[(model, label)] = read_gmmat(paths, model, label) if kind == "gmmat" else read_plink(paths, model, label)
    sample_counts = {
        "SNP": sum(1 for _ in (args.root / "01_structure/plink/coil144.SNP.analysis.fam").open()),
        "INDEL": sum(1 for _ in (args.root / "01_structure/plink/coil144.INDEL.analysis.fam").open()),
        "SV": sum(1 for _ in (args.root / "03_sv_gwas/plink/coil144.SV.GQ20.analysis.fam").open()),
    }
    write_tables(frames, args.outdir, sample_counts)
    draw_qq(frames, args.outdir)
    draw_manhattan(frames, args.outdir)


if __name__ == "__main__":
    main()
