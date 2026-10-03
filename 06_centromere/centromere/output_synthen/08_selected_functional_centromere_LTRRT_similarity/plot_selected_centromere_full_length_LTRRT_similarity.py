#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch


TYPE_ORDER = ["Copia_SIRE", "Gypsy_Athila", "Gypsy_CRM", "Gypsy_Ogre", "Gypsy_Tekay"]
TYPE_COLORS = {
    "Copia_SIRE": "#8c564b",
    "Gypsy_Athila": "#2ca02c",
    "Gypsy_CRM": "#17becf",
    "Gypsy_Ogre": "#4daf4a",
    "Gypsy_Tekay": "#d95f5f",
}
SPECIES_ORDER = ["R108", "A17", "Mpo", "Msa", "474"]
SPECIES_COLORS = {
    "R108": "#1f77b4",
    "A17": "#ff7f0e",
    "Mpo": "#2ca02c",
    "Msa": "#9467bd",
    "474": "#8c564b",
}


def read_fasta_ids(path):
    ids = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def build_matrix(ids, blast_path, min_cov):
    n = len(ids)
    idx = {seq_id: i for i, seq_id in enumerate(ids)}
    mat = np.zeros((n, n), dtype=float)
    np.fill_diagonal(mat, 1.0)
    best = {}
    with open(blast_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            q, s = f[0], f[1]
            if q not in idx or s not in idx:
                continue
            pident = float(f[2]) / 100.0
            length = float(f[3])
            qlen = float(f[4])
            slen = float(f[5])
            cov = length / min(qlen, slen) if min(qlen, slen) else 0
            if cov < min_cov:
                continue
            key = (q, s)
            score = pident
            if score > best.get(key, 0):
                best[key] = score
    for (q, s), score in best.items():
        mat[idx[q], idx[s]] = score
    mat = np.maximum(mat, mat.T)
    return mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--blast", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-cov", type=float, default=0.2)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ids = read_fasta_ids(args.fasta)
    meta = pd.read_csv(args.metadata, sep="\t")
    meta = meta.set_index("seq_id").loc[ids].reset_index()
    meta["label"] = meta["label"].astype(str).str.strip()
    meta["selected_type"] = meta["selected_type"].astype(str).str.strip()
    meta["superfamily"] = meta["superfamily"].astype(str).str.strip()
    meta["clade"] = meta["clade"].astype(str).str.strip()
    meta["type_rank"] = meta["selected_type"].map({v: i for i, v in enumerate(TYPE_ORDER)}).fillna(999).astype(int)
    meta["species_rank"] = meta["label"].map({v: i for i, v in enumerate(SPECIES_ORDER)}).fillna(999).astype(int)
    meta["time_sort"] = pd.to_numeric(meta["insertion_time_mya"], errors="coerce").fillna(999)
    meta = meta.sort_values(["type_rank", "species_rank", "time_sort", "seq_id"]).reset_index(drop=True)
    ordered_ids = meta["seq_id"].tolist()
    old_idx = {seq_id: i for i, seq_id in enumerate(ids)}
    order = [old_idx[i] for i in ordered_ids]

    mat = build_matrix(ids, args.blast, args.min_cov)
    mat = mat[np.ix_(order, order)]

    pd.DataFrame(mat, index=ordered_ids, columns=ordered_ids).to_csv(
        outdir / "selected_functional_centromere_full_length_LTRRT.similarity_matrix.minCov20.tsv",
        sep="\t",
        float_format="%.4f",
    )
    meta.to_csv(outdir / "selected_functional_centromere_full_length_LTRRT.plot_order.tsv", sep="\t", index=False)

    n = len(meta)
    fig_h = max(7.5, n / 58)
    fig_w = max(9.0, n / 50)
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=200)
    gs = fig.add_gridspec(
        nrows=4,
        ncols=2,
        width_ratios=[1, 0.04],
        height_ratios=[0.18, 1.0, 0.07, 0.07],
        hspace=0.03,
        wspace=0.03,
    )
    ax_type = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0])
    ax_species = fig.add_subplot(gs[2, 0])
    ax_time = fig.add_subplot(gs[3, 0])
    cax = fig.add_subplot(gs[1, 1])

    cmap = LinearSegmentedColormap.from_list("sim_blue_white_red", ["#2166ac", "#f7f7f7", "#d7191c"])
    mask = np.triu(np.ones_like(mat, dtype=bool), k=1)
    plot_mat = np.ma.array(mat, mask=mask)
    im = ax.imshow(plot_mat, cmap=cmap, vmin=0, vmax=1, interpolation="nearest", aspect="auto")
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("intact LTR-RT similarity", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Functional-centromere selected full-length intact LTR-RT similarity", fontsize=10, pad=6)

    type_rgb = np.array([[plt.matplotlib.colors.to_rgb(TYPE_COLORS.get(t, "#999999")) for t in meta["selected_type"]]])
    ax_type.imshow(type_rgb, aspect="auto")
    ax_type.set_xticks([])
    ax_type.set_yticks([])
    ax_type.set_ylabel("type", fontsize=7, rotation=0, labelpad=16, va="center")

    species_rgb = np.array([[plt.matplotlib.colors.to_rgb(SPECIES_COLORS.get(s, "#999999")) for s in meta["label"]]])
    ax_species.imshow(species_rgb, aspect="auto")
    ax_species.set_xticks([])
    ax_species.set_yticks([])
    ax_species.set_ylabel("species", fontsize=7, rotation=0, labelpad=18, va="center")

    times = pd.to_numeric(meta["insertion_time_mya"], errors="coerce").fillna(0).to_numpy()[None, :]
    time_max = max(2.0, float(np.nanpercentile(times, 95)) if times.size else 2.0)
    time_cmap = LinearSegmentedColormap.from_list("time_yellow_red", ["#fff7bc", "#fec44f", "#d95f0e", "#a50f15"])
    ax_time.imshow(times, aspect="auto", cmap=time_cmap, vmin=0, vmax=time_max)
    ax_time.set_xticks([])
    ax_time.set_yticks([])
    ax_time.set_ylabel("Mya", fontsize=7, rotation=0, labelpad=18, va="center")

    type_handles = [Patch(color=TYPE_COLORS[t], label=f"{t} (n={(meta['selected_type'] == t).sum()})") for t in TYPE_ORDER if (meta["selected_type"] == t).any()]
    species_handles = [Patch(color=SPECIES_COLORS[s], label=f"{s} (n={(meta['label'] == s).sum()})") for s in SPECIES_ORDER if (meta["label"] == s).any()]
    fig.legend(handles=type_handles + species_handles, loc="lower center", ncol=5, fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.01))

    for ext in ["png", "pdf"]:
        fig.savefig(outdir / f"selected_functional_centromere_full_length_LTRRT.similarity_heatmap.minCov20.{ext}", bbox_inches="tight")
    plt.close(fig)

    summary = meta.groupby(["label", "superfamily", "clade", "selected_type"]).size().reset_index(name="count")
    summary.to_csv(outdir / "selected_functional_centromere_full_length_LTRRT.plot_counts.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
