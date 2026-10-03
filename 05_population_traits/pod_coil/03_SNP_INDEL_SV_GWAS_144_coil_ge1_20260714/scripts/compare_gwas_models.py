#!/usr/bin/env python3
import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def clean_score(path):
    frame = pd.read_csv(path, sep=r"\s+", low_memory=False)
    frame = frame.rename(columns={"CHR": "chrom", "SNP": "id", "POS": "pos", "PVAL": "p_value"})
    frame = frame[pd.to_numeric(frame.p_value, errors="coerce").between(0, 1, inclusive="neither")].copy()
    frame["chrom"] = pd.to_numeric(frame.chrom, errors="coerce")
    frame["pos"] = pd.to_numeric(frame.pos, errors="coerce")
    frame["p_value"] = pd.to_numeric(frame.p_value, errors="coerce")
    frame["maf"] = np.minimum(pd.to_numeric(frame.AF, errors="coerce"), 1 - pd.to_numeric(frame.AF, errors="coerce"))
    frame["effect"] = pd.to_numeric(frame.SCORE, errors="coerce") / pd.to_numeric(frame.VAR, errors="coerce")
    frame["se"] = np.sqrt(1 / pd.to_numeric(frame.VAR, errors="coerce"))
    frame["effect_definition"] = "score_beta_approx"
    return frame[["chrom", "pos", "id", "p_value", "maf", "effect", "se", "effect_definition"]]


def clean_plink(paths, frequency):
    frames = []
    for path in paths:
        frame = pd.read_csv(path, sep=r"\s+", low_memory=False)
        frame = frame[frame.TEST.eq("ADD")].copy()
        frames.append(frame)
    frame = pd.concat(frames, ignore_index=True)
    frame = frame.rename(columns={"CHR": "chrom", "BP": "pos", "SNP": "id", "P": "p_value"})
    frame["chrom"] = pd.to_numeric(frame.chrom, errors="coerce")
    frame["pos"] = pd.to_numeric(frame.pos, errors="coerce")
    frame["p_value"] = pd.to_numeric(frame.p_value, errors="coerce")
    frame["effect"] = np.log(pd.to_numeric(frame.OR, errors="coerce"))
    frame["se"] = pd.to_numeric(frame.SE, errors="coerce")
    freq = pd.read_csv(frequency, sep=r"\s+")[["SNP", "MAF"]].rename(columns={"SNP": "id", "MAF": "maf"})
    frame = frame.merge(freq, on="id", how="left")
    frame["effect_definition"] = "log_OR"
    return frame[frame.p_value.between(0, 1, inclusive="neither")][["chrom", "pos", "id", "p_value", "maf", "effect", "se", "effect_definition"]]


def gmmat_paths(root, model, label):
    if model == "A":
        if label in ("SNP", "INDEL"):
            return [root / "02_snp_indel_gwas" / "gmmat" / f"{label}.coil_ge1.pc5.gmmat.score.tsv"]
        paths = list((root / "03_sv_gwas").rglob("SV*.pc5.gmmat.score.tsv"))
        return [next(path for path in paths if not path.name.endswith("model_audit.tsv"))]
    if model == "C":
        base = root / "05_model_sensitivity_20260715" / "model_C_gmmat_grm_pc2" / label
        if label == "SV":
            return [base / "SV.grm_pc2.gmmat.score.tsv"]
        return [base / f"{label}.chr{chrom}.grm_pc2.gmmat.score.tsv" for chrom in range(1, 9)]
    raise ValueError(model)


def plink_paths(root, label):
    base = root / "05_model_sensitivity_20260715" / "model_B_plink_pc3" / label
    if label == "SV":
        return [base / "SV.pc3.assoc.logistic"]
    return [base / f"{label}.chr{chrom}.pc3.assoc.logistic" for chrom in range(1, 9)]


def write_top20(model, label, frame, out):
    top = frame.nsmallest(20, "p_value").copy()
    top.insert(0, "Model", model)
    top.insert(1, "Variant_class", label)
    top = top.rename(columns={"chrom": "CHR", "pos": "POS", "id": "SNP", "p_value": "P", "maf": "MAF", "effect": "Effect", "se": "SE"})
    top.to_csv(out / "top20" / f"Model_{model}_{label}_top20.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--freqdir", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    root = Path(args.root)
    freqdir = Path(args.freqdir)
    out = Path(args.outdir)
    (out / "top20").mkdir(parents=True, exist_ok=True)

    overlap_rows, correlation_rows = [], []
    for label in ("SNP", "INDEL", "SV"):
        frames = {}
        for model in ("A", "C"):
            parts = [clean_score(path) for path in gmmat_paths(root, model, label)]
            frames[model] = pd.concat(parts, ignore_index=True)
        frames["B"] = clean_plink(plink_paths(root, label), freqdir / f"{label}.frq")
        for model, frame in frames.items():
            if frame.empty:
                raise RuntimeError(f"No valid {label} tests for Model {model}")
            write_top20(model, label, frame, out)
        for a, b in combinations(("A", "B", "C"), 2):
            merged = frames[a][["id", "p_value"]].rename(columns={"p_value": "p_a"}).merge(
                frames[b][["id", "p_value"]].rename(columns={"p_value": "p_b"}), on="id", how="inner"
            )
            correlation = pearsonr(-np.log10(merged.p_a), -np.log10(merged.p_b)).statistic if len(merged) >= 3 else np.nan
            top_a = set(frames[a].nsmallest(20, "p_value").id)
            top_b = set(frames[b].nsmallest(20, "p_value").id)
            shared = len(top_a & top_b)
            overlap_rows.append({"Variant_class": label, "Model_a": a, "Model_b": b,
                                 "top20_shared": shared, "top20_jaccard": shared / len(top_a | top_b)})
            correlation_rows.append({"Variant_class": label, "Model_a": a, "Model_b": b,
                                     "shared_tested_variants": len(merged), "minus_log10P_pearson_r": correlation})
    pd.DataFrame(overlap_rows).to_csv(out / "Model_ABC_top20_overlap.tsv", sep="\t", index=False)
    pd.DataFrame(correlation_rows).to_csv(out / "Model_ABC_pvalue_correlation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
