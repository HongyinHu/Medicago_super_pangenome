#!/usr/bin/env python3
"""Create PCA and SNP-GRM diagnostics for the 144-sample coil GWAS."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pointbiserialr


COLORS = {0: "#4C78A8", 1: "#E45756"}
LABELS = {0: "Non-spiral", 1: "Spiral"}


def save(fig: plt.Figure, outdir: Path, stem: str) -> None:
    fig.savefig(outdir / f"{stem}.png", dpi=400, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    model = pd.read_csv(args.root / "01_structure/model/SNP.model.tsv", sep="\t")
    phenotype = pd.read_csv(args.root / "00_input/coil_ge1_phenotype.tsv", sep="\t")
    metadata = model.merge(
        phenotype[["IID", "Latin_name"]], left_on="id", right_on="IID", how="left", validate="one_to_one"
    )
    if metadata["Latin_name"].isna().any():
        raise RuntimeError("Some PCA samples lack Latin_name metadata")

    eigenvalues = np.loadtxt(args.root / "01_structure/plink/coil144.SNP.pca.eigenval")
    explained = eigenvalues / eigenvalues.sum() * 100.0

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.8))
    for ax, xpc, ypc in ((axes[0], "PC1", "PC2"), (axes[1], "PC1", "PC3")):
        for trait in (0, 1):
            subset = metadata[metadata.y == trait]
            ax.scatter(
                subset[xpc], subset[ypc], s=28, alpha=0.78, color=COLORS[trait],
                edgecolor="white", linewidth=0.35, label=f"{LABELS[trait]} (n={len(subset)})"
            )
        xnum, ynum = int(xpc[2:]) - 1, int(ypc[2:]) - 1
        ax.set_xlabel(f"{xpc} ({explained[xnum]:.1f}%)")
        ax.set_ylabel(f"{ypc} ({explained[ynum]:.1f}%)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    axes[2].bar(np.arange(1, len(explained) + 1), explained, color="#72B7B2")
    axes[2].set_xlabel("Principal component")
    axes[2].set_ylabel("Variance explained (%)")
    axes[2].set_xticks(np.arange(1, len(explained) + 1))
    axes[2].spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"SNP PCA: {len(metadata)} samples from {metadata.Latin_name.nunique()} Medicago taxa", y=1.02)
    fig.tight_layout()
    save(fig, args.outdir, "PCA_trait_structure")

    grm_ids = pd.read_csv(
        args.root / "01_structure/plink/coil144.SNP.grm.rel.id", sep=r"\s+", header=None,
        names=["FID", "IID"], dtype=str
    )
    grm = np.loadtxt(args.root / "01_structure/plink/coil144.SNP.grm.rel")
    if grm.shape != (len(grm_ids), len(grm_ids)):
        raise RuntimeError(f"GRM shape {grm.shape} does not match {len(grm_ids)} IDs")
    lookup = metadata.set_index("id")
    y = np.array([int(lookup.loc[iid, "y"]) for iid in grm_ids.IID])
    pc2 = np.array([float(lookup.loc[iid, "PC2"]) for iid in grm_ids.IID])
    order = np.lexsort((pc2, y))
    ordered = grm[np.ix_(order, order)]
    ordered_y = y[order]

    fig = plt.figure(figsize=(7.5, 6.4))
    grid = fig.add_gridspec(2, 2, width_ratios=[0.18, 1], height_ratios=[0.08, 1], wspace=0.04, hspace=0.04)
    ax_top = fig.add_subplot(grid[0, 1])
    ax_left = fig.add_subplot(grid[1, 0])
    ax = fig.add_subplot(grid[1, 1])
    trait_rgb = np.array([matplotlib.colors.to_rgb(COLORS[value]) for value in ordered_y])
    ax_top.imshow(trait_rgb[np.newaxis, :, :], aspect="auto")
    ax_left.imshow(trait_rgb[:, np.newaxis, :], aspect="auto")
    ax_top.set_axis_off()
    ax_left.set_axis_off()
    image = ax.imshow(ordered, cmap="RdBu_r", aspect="equal", interpolation="nearest")
    ax.set_xlabel("Samples ordered by phenotype and PC2")
    ax.set_ylabel("Samples ordered by phenotype and PC2")
    ax.set_xticks([])
    ax.set_yticks([])
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
    colorbar.set_label("SNP relationship coefficient")
    fig.suptitle("SNP-derived genomic relationship matrix")
    save(fig, args.outdir, "GRM_heatmap")

    tri_i, tri_j = np.triu_indices(len(grm), k=1)
    values = grm[tri_i, tri_j]
    pair_class = np.where(
        (y[tri_i] == 1) & (y[tri_j] == 1), "Spiral-Spiral",
        np.where((y[tri_i] == 0) & (y[tri_j] == 0), "Nonspiral-Nonspiral", "Between phenotypes")
    )
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    palette = {"Spiral-Spiral": "#E45756", "Nonspiral-Nonspiral": "#4C78A8", "Between phenotypes": "#72B7B2"}
    for group in ("Spiral-Spiral", "Nonspiral-Nonspiral", "Between phenotypes"):
        subset = values[pair_class == group]
        ax.hist(subset, bins=60, density=True, histtype="step", linewidth=1.5, color=palette[group], label=f"{group} (n={len(subset):,})")
    ax.set_xlabel("Pairwise SNP relationship coefficient")
    ax.set_ylabel("Density")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, args.outdir, "GRM_pairwise_distribution")

    with (args.outdir / "structure_diagnostics.tsv").open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"samples\t{len(metadata)}\n")
        handle.write(f"spiral_cases\t{int((metadata.y == 1).sum())}\n")
        handle.write(f"nonspiral_controls\t{int((metadata.y == 0).sum())}\n")
        handle.write(f"taxa\t{metadata.Latin_name.nunique()}\n")
        handle.write(f"grm_offdiag_median\t{np.median(values):.8g}\n")
        handle.write(f"grm_offdiag_q025\t{np.quantile(values, 0.025):.8g}\n")
        handle.write(f"grm_offdiag_q975\t{np.quantile(values, 0.975):.8g}\n")
        for pc in range(1, 6):
            name = f"PC{pc}"
            correlation = pointbiserialr(metadata.y, metadata[name])
            cases = metadata.loc[metadata.y == 1, name]
            controls = metadata.loc[metadata.y == 0, name]
            mw = mannwhitneyu(cases, controls, alternative="two-sided")
            handle.write(f"{name}_point_biserial_r\t{correlation.statistic:.8g}\n")
            handle.write(f"{name}_point_biserial_p\t{correlation.pvalue:.8g}\n")
            handle.write(f"{name}_mann_whitney_p\t{mw.pvalue:.8g}\n")


if __name__ == "__main__":
    main()
