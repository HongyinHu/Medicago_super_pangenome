#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


SPECIES_ORDER = ["R108", "A17", "Mpo", "Msa", "474"]
SPECIES_COLORS = {
    "R108": "#2c6db2",
    "A17": "#d07a2d",
    "Mpo": "#5b9c48",
    "Msa": "#8e5aa6",
    "474": "#6f6f6f",
}

GROUP_ORDER = [
    "Copia",
    "Gypsy_CRM",
    "Gypsy_Athila",
    "Gypsy_Ogre",
    "Gypsy_Reina",
    "Gypsy_Retand",
    "Gypsy_Tekay",
]
GROUP_COLORS = {
    "Copia": "#9c9c9c",
    "Gypsy_CRM": "#45b8d2",
    "Gypsy_Athila": "#97c95c",
    "Gypsy_Ogre": "#2f9c95",
    "Gypsy_Reina": "#6b79c9",
    "Gypsy_Retand": "#b879c9",
    "Gypsy_Tekay": "#f0aaa0",
    "Other": "#d8d8d8",
}


def major_group(row):
    superfamily = str(row.get("superfamily", ""))
    clade = str(row.get("clade", ""))
    if superfamily == "Copia":
        return "Copia"
    if superfamily == "Gypsy":
        group = "Gypsy_" + clade
        return group if group in GROUP_COLORS else "Other"
    return "Other"


def cluster_order_for_ids(matrix, ids):
    if len(ids) <= 2:
        return list(ids)
    from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
    from scipy.spatial.distance import squareform

    sub = matrix.loc[ids, ids].to_numpy(dtype=float)
    dist = 1.0 - sub
    dist[dist < 0] = 0
    np.fill_diagonal(dist, 0)
    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    try:
        z = optimal_leaf_ordering(z, condensed)
    except Exception:
        pass
    return [ids[i] for i in leaves_list(z)]


def clade_guided_order(matrix, meta):
    ordered = []
    for group in GROUP_ORDER + ["Other"]:
        ids = meta.loc[meta["major_group"] == group, "seq_id"].tolist()
        ordered.extend(cluster_order_for_ids(matrix, ids))
    return ordered


def global_cluster_order(matrix):
    return cluster_order_for_ids(matrix, list(matrix.index))


def tree_leaves(treefile):
    text = Path(treefile).read_text().strip()
    return re.findall(r"(?<=[(,])([^():,;]+):", text)


def current_coord_map(ids):
    coord_map = {}
    pattern = re.compile(r"^([^_]+)__([^_]+)__LTR([0-9]+)_([0-9]+)__")
    for seq_id in ids:
        match = pattern.match(seq_id)
        if match:
            label, chrom, start, end = match.groups()
            coord_map[(label, chrom, start, end)] = seq_id
    return coord_map


def mapped_tree_order(matrix, treefile, tree_id_mode):
    matrix_ids = list(matrix.index)
    matrix_id_set = set(matrix_ids)
    leaves = tree_leaves(treefile)
    order = []
    seen = set()
    coord_map = current_coord_map(matrix_ids)
    old_pattern = re.compile(r"^([^_]+)__.*__(Chr[^_]+)__([0-9]+)_([0-9]+)__")

    for leaf in leaves:
        seq_id = leaf
        if tree_id_mode == "coord":
            match = old_pattern.match(leaf)
            if not match:
                continue
            seq_id = coord_map.get(match.groups())
        if seq_id in matrix_id_set and seq_id not in seen:
            order.append(seq_id)
            seen.add(seq_id)
    order.extend([seq_id for seq_id in matrix_ids if seq_id not in seen])
    return order


def move_group_first(order, meta, group):
    group_by_id = meta.set_index("seq_id")["major_group"].to_dict()
    left = [seq_id for seq_id in order if group_by_id.get(seq_id) == group]
    right = [seq_id for seq_id in order if group_by_id.get(seq_id) != group]
    return left + right


def color_row(values, color_map):
    return np.array([[to_rgb(color_map.get(str(v), "#dddddd")) for v in values]], dtype=float)


def draw_bar(ax, colors, label):
    ax.imshow(colors, aspect="auto", interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8, labelpad=18)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_heatmap(matrix, meta, out_prefix, title, age_max, similarity_label, similarity_white):
    n = matrix.shape[0]
    values = matrix.to_numpy(dtype=float)
    masked = np.ma.array(values, mask=np.triu(np.ones_like(values, dtype=bool), k=1))
    grid = np.arange(n + 1, dtype=float)
    rr, cc = np.meshgrid(grid, grid, indexing="ij")
    x = (rr + cc) / 2.0 - 0.5
    y = (rr - cc) / 2.0

    fig_w = max(10.5, n / 34.0)
    fig_h = fig_w * 0.74
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)
    gs = fig.add_gridspec(
        nrows=5,
        ncols=1,
        height_ratios=[6.2, 0.17, 0.17, 0.17, 1.6],
        hspace=0.045,
    )
    ax_heat = fig.add_subplot(gs[0, 0])
    ax_type = fig.add_subplot(gs[1, 0], sharex=ax_heat)
    ax_species = fig.add_subplot(gs[2, 0], sharex=ax_heat)
    ax_age = fig.add_subplot(gs[3, 0], sharex=ax_heat)
    ax_leg = fig.add_subplot(gs[4, 0])

    similarity_white = min(max(float(similarity_white), 0.01), 0.99)
    sim_cmap = LinearSegmentedColormap.from_list(
        "similarity",
        [
            (0.0, "#6D89BA"),
            (similarity_white * 0.65, "#a8cce4"),
            (similarity_white, "#f7f7f7"),
            (similarity_white + (1.0 - similarity_white) * 0.45, "#f4a6a6"),
            (1.0, "#d7191c"),
        ],
    )
    im = ax_heat.pcolormesh(x, y, masked, cmap=sim_cmap, vmin=0, vmax=1, shading="flat", rasterized=True)
    ax_heat.set_aspect("equal")
    ax_heat.set_xlim(-0.5, n - 0.5)
    ax_heat.set_ylim(-0.5, n / 2 + 2)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    ax_heat.set_title(title, fontsize=13, pad=8)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    sim_cax = inset_axes(ax_heat, width="23%", height="5%", loc="upper left", borderpad=1.0)
    cb = fig.colorbar(im, cax=sim_cax, orientation="horizontal")
    cb.set_label(similarity_label, fontsize=8, labelpad=2)
    cb.ax.tick_params(labelsize=7, length=2)

    draw_bar(ax_type, color_row(meta["major_group"], GROUP_COLORS), "TE type")
    draw_bar(ax_species, color_row(meta["label"], SPECIES_COLORS), "species")

    ages = pd.to_numeric(meta["insertion_time_mya"], errors="coerce").fillna(0.0).to_numpy()[None, :]
    vmax = age_max if age_max else max(1.0, float(np.nanpercentile(ages, 98)))
    age_cmap = LinearSegmentedColormap.from_list("ltr_age", ["#fff7df", "#f6c16b", "#de5b3d", "#9d1d20"])
    ax_age.imshow(ages, aspect="auto", interpolation="nearest", cmap=age_cmap, norm=Normalize(0, vmax))
    ax_age.set_xticks([])
    ax_age.set_yticks([])
    ax_age.set_ylabel("LTR age", rotation=0, ha="right", va="center", fontsize=8, labelpad=18)
    for spine in ax_age.spines.values():
        spine.set_visible(False)

    age_cax = inset_axes(
        ax_heat, width="23%", height="5%", loc="upper left", borderpad=1.0,
        bbox_to_anchor=(0.0, -0.12, 1, 1), bbox_transform=ax_heat.transAxes
    )
    sm = plt.cm.ScalarMappable(norm=Normalize(0, vmax), cmap=age_cmap)
    age_cb = fig.colorbar(sm, cax=age_cax, orientation="horizontal")
    age_cb.set_label("LTR insertion time (Mya)", fontsize=8, labelpad=2)
    age_cb.ax.tick_params(labelsize=7, length=2)

    ax_leg.axis("off")
    type_handles = [
        Patch(color=GROUP_COLORS[g], label="{} (n={})".format(g, int((meta["major_group"] == g).sum())))
        for g in GROUP_ORDER + ["Other"]
        if (meta["major_group"] == g).any()
    ]
    species_handles = [
        Patch(color=SPECIES_COLORS[s], label="{} (n={})".format(s, int((meta["label"] == s).sum())))
        for s in SPECIES_ORDER
        if (meta["label"] == s).any()
    ]
    leg1 = ax_leg.legend(handles=type_handles, loc="upper center", ncol=min(4, len(type_handles)),
                         fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.92), title="TE clade")
    leg1.get_title().set_fontsize(8)
    leg2 = ax_leg.legend(handles=species_handles, loc="lower center", ncol=min(5, len(species_handles)),
                         fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.06), title="Species")
    leg2.get_title().set_fontsize(8)
    ax_leg.add_artist(leg1)

    out_prefix = Path(out_prefix)
    fig.savefig(str(out_prefix) + ".png", bbox_inches="tight")
    fig.savefig(str(out_prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(out_prefix) + ".svg", bbox_inches="tight")
    plt.close(fig)


def write_outputs(matrix, meta, order, out_prefix, title, age_max, similarity_label, similarity_white):
    ordered_matrix = matrix.loc[order, order]
    ordered_meta = meta.set_index("seq_id").loc[order].reset_index()
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    ordered_matrix.to_csv(str(out_prefix) + ".matrix_0to1.tsv", sep="\t", float_format="%.6f")
    (ordered_matrix * 100).to_csv(str(out_prefix) + ".matrix_0to100_TBtools.tsv", sep="\t", float_format="%.2f")
    ordered_meta.to_csv(str(out_prefix) + ".metadata.tsv", sep="\t", index=False)
    ordered_meta[["seq_id", "major_group", "label", "insertion_time_mya"]].to_csv(
        str(out_prefix) + ".annotation_for_TBtools.tsv", sep="\t", index=False
    )
    draw_heatmap(ordered_matrix, ordered_meta, out_prefix, title, age_max, similarity_label, similarity_white)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--order-mode", choices=["cluster", "clade_guided", "tree"], default="clade_guided")
    parser.add_argument("--treefile")
    parser.add_argument("--tree-id-mode", choices=["exact", "coord"], default="exact")
    parser.add_argument("--first-group", help="Move this major_group to the left after ordering, for example Copia.")
    parser.add_argument("--title", default="Functional-centromere RT-domain intact LTR-RT full-length similarity")
    parser.add_argument("--similarity-label", default="intact LTR-RT similarity")
    parser.add_argument("--similarity-white", type=float, default=0.5)
    parser.add_argument("--age-max", type=float, default=1.0)
    args = parser.parse_args()

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    meta = pd.read_csv(args.metadata, sep="\t")
    meta["seq_id"] = meta["seq_id"].astype(str)
    meta["label"] = meta["label"].astype(str).str.strip()
    meta["superfamily"] = meta["superfamily"].astype(str).str.strip()
    meta["clade"] = meta["clade"].astype(str).str.strip()
    meta["major_group"] = meta.apply(major_group, axis=1)
    common = [seq_id for seq_id in matrix.index if seq_id in set(meta["seq_id"])]
    matrix = matrix.loc[common, common]
    meta = meta.set_index("seq_id").loc[common].reset_index()

    if args.order_mode == "cluster":
        order = global_cluster_order(matrix)
    elif args.order_mode == "tree":
        if not args.treefile:
            raise SystemExit("--treefile is required when --order-mode tree")
        order = mapped_tree_order(matrix, args.treefile, args.tree_id_mode)
    else:
        order = clade_guided_order(matrix, meta)
    if args.first_group:
        order = move_group_first(order, meta, args.first_group)
    write_outputs(
        matrix, meta, order, args.out_prefix, args.title, args.age_max,
        args.similarity_label, args.similarity_white
    )
    print("mode={}".format(args.order_mode))
    print("n={}".format(len(order)))
    if args.order_mode == "tree":
        print("treefile={}".format(args.treefile))
    print(meta.groupby(["label", "major_group"]).size().reset_index(name="count").to_string(index=False))


if __name__ == "__main__":
    main()
