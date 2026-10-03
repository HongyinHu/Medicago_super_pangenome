#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import squareform


SPECIES_COLORS = {
    "A17": "#008751",
    "R108": "#ff6a00",
    "Mpo": "#3b6fb6",
    "Msa": "#8a4fb3",
    "474": "#c83f49",
}
GROUP_RANK = {"Copia": 0, "Gypsy": 1, "LTR_unknown": 2, "unknown": 9}
CLADE_COLORS = {
    "Copia": "#bdbdbd",
    "Ale": "#7fc97f",
    "Alesia": "#4daf4a",
    "Angela": "#80b1d3",
    "Bianca": "#fdb462",
    "Bryco": "#b3de69",
    "Ikeros": "#fccde5",
    "Ivana": "#bebada",
    "SIRE": "#fb8072",
    "TAR": "#8dd3c7",
    "Tork": "#bc80bd",
    "Gypsy": "#d9d9d9",
    "Athila": "#1b9e77",
    "CRM": "#1f78b4",
    "Ogre": "#33a02c",
    "Reina": "#ff7f00",
    "Retand": "#6a3d9a",
    "Tekay": "#e31a1c",
    "LTR_unknown": "#969696",
    "unclassified": "#f0f0f0",
}


def read_fasta_lengths(path):
    lengths = {}
    header = None
    length = 0
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    lengths[header] = length
                header = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
        if header is not None:
            lengths[header] = length
    return lengths


def read_tesorter_cls(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fields = line.rstrip("\n").split(None, 6)
            while len(fields) < 7:
                fields.append("")
            rows.append(fields[:7])
    if not rows:
        return pd.DataFrame(columns=["id", "ts_order", "ts_superfamily", "ts_clade", "complete", "ts_strand", "domains"])
    cls = pd.DataFrame(rows, columns=["raw_id", "ts_order", "ts_superfamily", "ts_clade", "complete", "ts_strand", "domains"])
    cls["id"] = cls["raw_id"].astype(str).str.split("#", n=1).str[0]
    cls["ts_clade"] = cls["ts_clade"].replace({"": "unclassified", "-": "unclassified"}).fillna("unclassified")
    cls["ts_superfamily"] = cls["ts_superfamily"].replace({"": "unknown", "-": "unknown"}).fillna("unknown")
    return cls.drop_duplicates("id")


def assign_group(row):
    classification = str(row.get("classification", "unknown"))
    if row.get("ts_order") == "LTR":
        superfamily = row.get("ts_superfamily")
        clade = row.get("ts_clade")
        if superfamily in {"Copia", "Gypsy"} and clade not in {"", "unknown", "unclassified"}:
            return superfamily, clade
    if classification == "LTR/Copia":
        return "Copia", "Copia"
    if classification == "LTR/Gypsy":
        return "Gypsy", "Gypsy"
    return "LTR_unknown", "LTR_unknown"


def apply_display_grouping(meta, collapse_copia=False, gypsy_top_n=0):
    meta = meta.copy()
    meta["plot_group_display"] = meta["plot_group"]
    meta["fine_type_display"] = meta["fine_type"]
    if collapse_copia:
        is_copia = meta["plot_group"] == "Copia"
        meta.loc[is_copia, "fine_type_display"] = "Copia"
    if gypsy_top_n and gypsy_top_n > 0:
        gypsy_counts = (
            meta[meta["plot_group"] == "Gypsy"]["fine_type"]
            .value_counts()
            .head(gypsy_top_n)
        )
        top_gypsy = set(gypsy_counts.index)
        is_minor_gypsy = (meta["plot_group"] == "Gypsy") & (~meta["fine_type"].isin(top_gypsy))
        meta.loc[is_minor_gypsy, "fine_type_display"] = "Gypsy_other"
    return meta


def build_similarity_from_paf(ids, paf_path, min_cov):
    n = len(ids)
    index = {seq_id: i for i, seq_id in enumerate(ids)}
    sim = np.zeros((n, n), dtype=np.float32)
    np.fill_diagonal(sim, 1.0)
    with open(paf_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            q, t = fields[0], fields[5]
            if q not in index or t not in index:
                continue
            qlen = int(fields[1])
            tlen = int(fields[6])
            matches = int(fields[9])
            block_len = int(fields[10])
            if block_len <= 0:
                continue
            cov = block_len / max(1, min(qlen, tlen))
            if cov < min_cov:
                continue
            value = matches / block_len
            i, j = index[q], index[t]
            if value > sim[i, j]:
                sim[i, j] = value
                sim[j, i] = value
    return sim


def build_similarity_from_blastn(ids, blastn_path, min_cov=0.0, min_aln_len=0):
    n = len(ids)
    index = {seq_id: i for i, seq_id in enumerate(ids)}
    sim = np.zeros((n, n), dtype=np.float32)
    np.fill_diagonal(sim, 1.0)
    with open(blastn_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 14:
                continue
            q, t = fields[0], fields[1]
            if q not in index or t not in index:
                continue
            pident = float(fields[2]) / 100.0
            aln_len = int(fields[3])
            qlen = int(fields[12])
            tlen = int(fields[13])
            if aln_len < min_aln_len:
                continue
            cov = aln_len / max(1, min(qlen, tlen))
            if cov < min_cov:
                continue
            i, j = index[q], index[t]
            if pident > sim[i, j]:
                sim[i, j] = pident
                sim[j, i] = pident
    return sim


def color_strip(labels, color_map):
    colors = [mpl.colors.to_rgb(color_map.get(label, "#cccccc")) for label in labels]
    return np.array(colors, dtype=float).reshape(1, len(labels), 3)


def plot(meta, sim, out_prefix):
    dist = 1.0 - sim
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)
    condensed = squareform(dist, checks=False)
    link = linkage(condensed, method="average")
    order = leaves_list(link)

    meta_o = meta.iloc[order].reset_index(drop=True)
    sim_o = sim[np.ix_(order, order)]
    n = len(meta_o)
    mask = np.triu(np.ones_like(sim_o, dtype=bool), k=1)
    sim_masked = np.ma.array(sim_o, mask=mask)

    fig = plt.figure(figsize=(max(10, n * 0.04), 10.5))
    gs = fig.add_gridspec(
        nrows=7,
        ncols=2,
        width_ratios=[28, 1.2],
        height_ratios=[26, 0.55, 0.55, 0.55, 4.5, 1.3, 1.3],
        hspace=0.04,
        wspace=0.04,
    )

    ax_heat = fig.add_subplot(gs[0, 0])
    cmap = mpl.cm.get_cmap("RdBu_r").copy()
    cmap.set_bad("white")
    im = ax_heat.imshow(sim_masked, cmap=cmap, vmin=0, vmax=1, interpolation="nearest", rasterized=True)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    ax_cbar = fig.add_subplot(gs[0, 1])
    cb = fig.colorbar(im, cax=ax_cbar)
    cb.set_label("intact LTR-RT similarity", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    species_labels = meta_o["label"].tolist()
    clade_labels = meta_o["fine_type_display"].tolist()
    insert = meta_o["insertion_time_mya"].astype(float).fillna(0).to_numpy().reshape(1, n)

    ax_species = fig.add_subplot(gs[1, 0])
    ax_species.imshow(color_strip(species_labels, SPECIES_COLORS), aspect="auto")
    ax_species.set_xticks([])
    ax_species.set_yticks([])
    ax_species.set_ylabel("species", rotation=0, ha="right", va="center", fontsize=8)

    ax_clade = fig.add_subplot(gs[2, 0])
    ax_clade.imshow(color_strip(clade_labels, CLADE_COLORS), aspect="auto")
    ax_clade.set_xticks([])
    ax_clade.set_yticks([])
    ax_clade.set_ylabel("clade", rotation=0, ha="right", va="center", fontsize=8)

    ax_time = fig.add_subplot(gs[3, 0])
    time_im = ax_time.imshow(insert, aspect="auto", cmap="YlOrRd", vmin=0, vmax=max(2.0, np.nanpercentile(insert, 99)))
    ax_time.set_xticks([])
    ax_time.set_yticks([])
    ax_time.set_ylabel("time", rotation=0, ha="right", va="center", fontsize=8)
    ax_time_cbar = fig.add_subplot(gs[3, 1])
    time_cb = fig.colorbar(time_im, cax=ax_time_cbar)
    time_cb.set_label("LTR insert time (Mya)", fontsize=8)
    time_cb.ax.tick_params(labelsize=7)

    ax_den = fig.add_subplot(gs[4, 0])
    dendrogram(link, ax=ax_den, no_labels=True, color_threshold=0, above_threshold_color="black")
    ax_den.set_xticks([])
    ax_den.set_yticks([])
    for spine in ax_den.spines.values():
        spine.set_visible(False)

    ax_legend = fig.add_subplot(gs[5:, 0])
    ax_legend.axis("off")
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)
    species_handles = [
        mpl.patches.Patch(color=SPECIES_COLORS[label], label=f"{label} (n={(meta_o['label'] == label).sum()})")
        for label in SPECIES_COLORS
        if label in set(species_labels)
    ]
    clade_order = (
        meta_o[["plot_group_display", "fine_type_display"]]
        .drop_duplicates()
        .assign(rank=lambda d: d["plot_group_display"].map(GROUP_RANK).fillna(8))
        .sort_values(["rank", "plot_group_display", "fine_type_display"])
    )["fine_type_display"].tolist()
    clade_handles = [mpl.patches.Patch(color=CLADE_COLORS.get(label, "#cccccc"), label=label) for label in clade_order]
    leg1 = ax_legend.legend(
        handles=species_handles,
        loc="upper left",
        bbox_to_anchor=(0, 0.98),
        ncol=min(5, len(species_handles)),
        frameon=False,
        fontsize=8,
    )
    ax_legend.add_artist(leg1)
    ax_legend.legend(
        handles=clade_handles,
        loc="upper left",
        bbox_to_anchor=(0, 0.55),
        ncol=min(7, max(1, len(clade_handles))),
        frameon=False,
        fontsize=7,
    )

    fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    ordered_ids = meta_o[
        [
            "id",
            "genome",
            "label",
            "classification",
            "plot_group",
            "fine_type",
            "plot_group_display",
            "fine_type_display",
            "insertion_time_mya",
        ]
    ]
    ordered_ids.to_csv(f"{out_prefix}.ordered_metadata.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--paf")
    parser.add_argument("--blastn")
    parser.add_argument("--tesorter_cls", required=True)
    parser.add_argument("--out_prefix", required=True)
    parser.add_argument("--min_cov", type=float, default=0.5)
    parser.add_argument("--blast_min_cov", type=float, default=0.0)
    parser.add_argument("--blast_min_aln_len", type=int, default=0)
    parser.add_argument("--collapse_copia", action="store_true")
    parser.add_argument("--gypsy_top_n", type=int, default=0)
    args = parser.parse_args()

    meta = pd.read_csv(args.metadata, sep="\t")
    cls = read_tesorter_cls(args.tesorter_cls)
    meta = meta.merge(cls[["id", "ts_order", "ts_superfamily", "ts_clade"]], on="id", how="left")
    labels = meta.apply(assign_group, axis=1)
    meta["plot_group"] = [x[0] for x in labels]
    meta["fine_type"] = [x[1] for x in labels]
    meta = apply_display_grouping(meta, collapse_copia=args.collapse_copia, gypsy_top_n=args.gypsy_top_n)
    meta["insertion_time_mya"] = pd.to_numeric(meta["insertion_time_mya"], errors="coerce")

    ids = meta["id"].tolist()
    if args.blastn:
        sim = build_similarity_from_blastn(
            ids,
            args.blastn,
            min_cov=args.blast_min_cov,
            min_aln_len=args.blast_min_aln_len,
        )
    elif args.paf:
        sim = build_similarity_from_paf(ids, args.paf, args.min_cov)
    else:
        raise SystemExit("Provide either --paf or --blastn")
    pd.DataFrame(sim, index=ids, columns=ids).to_csv(f"{args.out_prefix}.similarity_matrix.tsv", sep="\t")
    plot(meta, sim, args.out_prefix)
    print(f"Plotted {len(ids)} functional-centromere intact LTR-RTs")


if __name__ == "__main__":
    main()
