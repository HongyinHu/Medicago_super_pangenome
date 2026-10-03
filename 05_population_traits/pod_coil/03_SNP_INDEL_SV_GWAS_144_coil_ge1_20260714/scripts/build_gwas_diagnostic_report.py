#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def markdown_table(frame):
    if frame.empty:
        return "No rows available.\n"
    table = frame.copy().fillna("NA").astype(str)
    columns = list(table.columns)
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    lines.extend("| " + " | ".join(row) + " |" for row in table.itertuples(index=False, name=None))
    return "\n".join(lines) + "\n"


def resolve_column(frame, requested):
    """Resolve human-readable and underscore-separated versions of a header."""
    normalize = lambda value: "".join(char for char in str(value).lower() if char.isalnum())
    available = {normalize(column): column for column in frame.columns}
    key = normalize(requested)
    aliases = {
        "snplambdagc": ["lambdagcsnp"],
        "indellambdagc": ["lambdagcindel"],
        "svlambdagc": ["lambdagcsv"],
    }
    candidates = [key] + aliases.get(key, [])
    for candidate in candidates:
        if candidate in available:
            return available[candidate]
    if key not in available:
        raise KeyError(f"Missing column {requested!r}; available columns: {list(frame.columns)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    out = Path(args.outdir)
    maf = pd.read_csv(out / "results" / "Model_A_MAF_sensitivity.tsv", sep="\t")
    pc = pd.read_csv(out / "diagnostics" / "phenotype_pc_association.tsv", sep="\t")
    structure = pd.read_csv(out / "diagnostics" / "structure_diagnostics.tsv", sep="\t")
    top_overlap = pd.read_csv(out / "model_comparison" / "Model_ABC_top20_overlap.tsv", sep="\t")
    correlations = pd.read_csv(out / "model_comparison" / "Model_ABC_pvalue_correlation.tsv", sep="\t")
    pc2 = pc[pc.iloc[:, 0].astype(str).eq("PC2")]
    snp_lambda = resolve_column(maf, "SNP lambdaGC")
    indel_lambda = resolve_column(maf, "INDEL lambdaGC")
    maf_cutoff = resolve_column(maf, "MAF cutoff")
    recommended = maf[
        maf[snp_lambda].between(0.8, 1.2)
        & maf[indel_lambda].between(0.8, 1.2)
    ]
    if len(recommended):
        cutoff_text = f"MAF >= {recommended.iloc[0][maf_cutoff]:.2f}"
    else:
        cutoff_text = "no tested cutoff simultaneously calibrated SNP and INDEL to the predefined 0.8-1.2 range"
    report = [
        "# GWAS diagnostic report",
        "",
        "## Scope",
        "Binary pod-coil GWAS: case = spiral circles >=1 and control = <1. SNP/INDEL use 143 samples (103/40); SV uses 137 samples (100/37). Model A was retained throughout: GMMAT binomial logistic mixed model with PC1-PC5 and the SNP-derived GRM.",
        "",
        "## MAF sensitivity",
        markdown_table(maf),
        "MAF 0.05 is the original validated Model A baseline and was reused rather than recomputed. MAF 0.10 and 0.15 were independently rescanned after PLINK re-filtering and with the same threshold explicitly supplied to GMMAT.",
        "",
        "## Population structure",
        markdown_table(pc),
        "PC2 is the primary confounding concern when its phenotype association is strong; this limits causal interpretation even if a locus is statistically ranked.",
        "",
        "## GRM and relatedness",
        markdown_table(structure),
        "The heatmap, pairwise distribution and top-related-pair table are in `diagnostics/`. The GRM was retained because samples span multiple Medicago lineages and the trait is strongly structured.",
        "",
        "## FarmCPU and BLINK",
        "GAPIT was not installed in the validated R environment at run start. FarmCPU/BLINK were therefore not run, rather than substituting software or changing environments during this diagnostic experiment.",
        "",
        "## Model ranking stability",
        "### Top-20 overlap",
        markdown_table(top_overlap),
        "### Genome-wide -log10(P) correlation",
        markdown_table(correlations),
        "Top-20 marker tables, including CHR, POS, SNP, P, MAF, score-effect approximation/log(OR), and SE, are in `model_comparison/top20/` and each MAF result directory.",
        "",
        "## Interpretation and recommended reporting model",
        f"For common SNP/INDEL variants, the calibrated Model A sensitivity choice is **{cutoff_text}**. SV calibration must be interpreted separately because its lambda can remain conservative after filtering. No diagnostic result alone establishes a true association; candidate signals require replication, local LD/structural validation, and explicit consideration of trait-lineage confounding.",
        "",
        "## Reproducibility",
        "All runs used identical phenotype coding, sample-order checks, genome build, SNP-derived GRM, and GMMAT version. Only the MAF threshold changed. Each threshold has a separate filtered PLINK bfile, GMMAT audit, raw score files, QQ plot, Manhattan plot, and completion marker.",
    ]
    (out / "results" / "GWAS_diagnostic_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
