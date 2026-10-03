#!/usr/bin/env python3
"""Plot full-length intact LTR-RT similarity as a tree-ordered triangular heatmap."""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


TYPE_ORDER = ["Copia_SIRE", "Gypsy_Athila", "Gypsy_CRM", "Gypsy_Ogre", "Gypsy_Tekay"]
TYPE_COLORS = {
    "Copia_SIRE": "#9c9c9c",
    "Gypsy_Athila": "#97c95c",
    "Gypsy_CRM": "#45b8d2",
    "Gypsy_Ogre": "#2f9c95",
    "Gypsy_Tekay": "#f0aaa0",
}

SPECIES_ORDER = ["R108", "A17", "Mpo", "Msa", "474"]
SPECIES_COLORS = {
    "R108": "#2c6db2",
    "A17": "#d07a2d",
    "Mpo": "#5b9c48",
    "Msa": "#8e5aa6",
    "474": "#6f6f6f",
}


class Node:
    def __init__(self, name="", length=0.0, children=None):
        self.name = name
        self.length = length
        self.children = children or []
        self.x = 0.0
        self.depth = 0.0

    @property
    def is_leaf(self):
        return not self.children


def _read_label(text, i):
    start = i
    while i < len(text) and text[i] not in ":,();":
        i += 1
    return text[start:i].strip().strip("'\""), i


def _read_length(text, i):
    if i >= len(text) or text[i] != ":":
        return 0.0, i
    i += 1
    start = i
    while i < len(text) and text[i] not in ",();":
        i += 1
    raw = text[start:i].strip()
    try:
        return float(raw), i
    except ValueError:
        return 0.0, i


def parse_newick(text):
    text = re.sub(r"\s+", "", text.strip())

    def parse_subtree(i):
        if text[i] == "(":
            i += 1
            children = []
            while True:
                child, i = parse_subtree(i)
                children.append(child)
                if text[i] == ",":
                    i += 1
                    continue
                if text[i] == ")":
                    i += 1
                    break
            name, i = _read_label(text, i)
            length, i = _read_length(text, i)
            return Node(name=name, length=length, children=children), i
        name, i = _read_label(text, i)
        length, i = _read_length(text, i)
        return Node(name=name, length=length), i

    root, _ = parse_subtree(0)
    return root


def assign_tree_coordinates(root, order):
    x_lookup = {name: idx for idx, name in enumerate(order)}
    max_depth = 0.0

    def set_depth(node, depth):
        nonlocal max_depth
        node.depth = depth
        max_depth = max(max_depth, depth)
        for child in node.children:
            set_depth(child, depth + max(child.length, 0.0))

    def set_x(node):
        if node.is_leaf:
            node.x = x_lookup.get(node.name, np.nan)
            return
        for child in node.children:
            set_x(child)
        xs = [child.x for child in node.children if np.isfinite(child.x)]
        node.x = float(np.mean(xs)) if xs else np.nan

    set_depth(root, 0.0)
    set_x(root)
    return max_depth or 1.0


def leaves_in_order(root):
    leaves = []

    def visit(node):
        if node.is_leaf:
            leaves.append(node.name)
            return
        for child in node.children:
            visit(child)

    visit(root)
    return leaves


def leaf_set(root):
    return set(leaves_in_order(root))


def build_adjacency(root):
    adjacency = {}

    def add_edge(a, b, length):
        adjacency.setdefault(id(a), {"node": a, "edges": []})
        adjacency.setdefault(id(b), {"node": b, "edges": []})
        adjacency[id(a)]["edges"].append((id(b), max(float(length), 0.0)))
        adjacency[id(b)]["edges"].append((id(a), max(float(length), 0.0)))

    def visit(node):
        adjacency.setdefault(id(node), {"node": node, "edges": []})
        for child in node.children:
            add_edge(node, child, child.length)
            visit(child)

    visit(root)
    return adjacency


def component_leaves(adjacency, start_id, blocked_id):
    stack = [start_id]
    seen = set([blocked_id])
    leaves = set()
    while stack:
        node_id = stack.pop()
        if node_id in seen:
            continue
        seen.add(node_id)
        node = adjacency[node_id]["node"]
        neighbor_ids = [neighbor_id for neighbor_id, _ in adjacency[node_id]["edges"] if neighbor_id != blocked_id]
        if node.is_leaf:
            leaves.add(node.name)
        for neighbor_id, _ in adjacency[node_id]["edges"]:
            if neighbor_id not in seen:
                stack.append(neighbor_id)
    return leaves


def orient_from_edge(adjacency, node_id, parent_id, incoming_length):
    original = adjacency[node_id]["node"]
    children = []
    for neighbor_id, length in adjacency[node_id]["edges"]:
        if neighbor_id == parent_id:
            continue
        children.append(orient_from_edge(adjacency, neighbor_id, node_id, length))
    return Node(name=original.name, length=incoming_length, children=children)


def reroot_by_leaf_set(root, target_leaves):
    target_leaves = set(target_leaves)
    all_leaves = leaf_set(root)
    if not target_leaves:
        raise ValueError("No leaves matched the requested root clade.")
    if target_leaves == all_leaves:
        raise ValueError("The requested root clade includes all leaves.")

    adjacency = build_adjacency(root)
    best = None
    best_score = (-1, -1, 0)
    checked = set()
    for node_id, entry in adjacency.items():
        for neighbor_id, length in entry["edges"]:
            edge_key = tuple(sorted((node_id, neighbor_id)))
            if edge_key in checked:
                continue
            checked.add(edge_key)
            side = component_leaves(adjacency, node_id, neighbor_id)
            other = all_leaves - side
            for candidate_side, candidate_node, candidate_parent in [
                (side, node_id, neighbor_id),
                (other, neighbor_id, node_id),
            ]:
                if candidate_side == target_leaves:
                    best = (candidate_node, candidate_parent, length, True, len(candidate_side), len(all_leaves) - len(candidate_side))
                    break
                overlap = len(candidate_side & target_leaves)
                extra = len(candidate_side - target_leaves)
                score = (overlap, -extra, -abs(len(candidate_side) - len(target_leaves)))
                if score > best_score:
                    best_score = score
                    best = (candidate_node, candidate_parent, length, False, len(candidate_side), len(all_leaves) - len(candidate_side))
            if best and best[3]:
                break
        if best and best[3]:
            break

    if best is None:
        raise ValueError("Could not find a suitable root edge.")

    copia_node, other_node, edge_length, exact, copia_side_n, other_side_n = best
    half = edge_length / 2.0
    new_root = Node(name="rooted_on_requested_clade", length=0.0, children=[
        orient_from_edge(adjacency, copia_node, other_node, half),
        orient_from_edge(adjacency, other_node, copia_node, half),
    ])
    return new_root, {
        "exact": exact,
        "target_n": len(target_leaves),
        "root_side_n": copia_side_n,
        "other_side_n": other_side_n,
    }


def rotate_target_clade_first(root, target_leaves):
    target_leaves = set(target_leaves)

    def score_and_rotate(node):
        if node.is_leaf:
            return 1 if node.name in target_leaves else 0, 1
        scored_children = []
        target_count = 0
        leaf_count = 0
        for child in node.children:
            child_target_count, child_leaf_count = score_and_rotate(child)
            target_count += child_target_count
            leaf_count += child_leaf_count
            proportion = child_target_count / float(child_leaf_count) if child_leaf_count else 0.0
            scored_children.append((child, child_target_count, child_leaf_count, proportion))
        scored_children.sort(key=lambda item: (item[3], item[1], -item[2]), reverse=True)
        node.children = [item[0] for item in scored_children]
        return target_count, leaf_count

    score_and_rotate(root)


def draw_tree(ax, root, max_depth, n):
    def y(node):
        return max_depth - node.depth

    def draw_node(node):
        if node.is_leaf:
            return
        valid_children = [child for child in node.children if np.isfinite(child.x)]
        if not valid_children:
            return
        xs = [child.x for child in valid_children]
        ax.plot([min(xs), max(xs)], [y(node), y(node)], color="#2b201d", lw=0.45)
        for child in valid_children:
            ax.plot([child.x, child.x], [y(node), y(child)], color="#2b201d", lw=0.45)
            draw_node(child)

    draw_node(root)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(-0.02 * max_depth, max_depth * 1.03)
    ax.axis("off")


def as_color_row(values, color_map):
    return np.array([[to_rgb(color_map.get(str(v), "#dddddd")) for v in values]], dtype=float)


def draw_annotation_bar(ax, colors, label):
    ax.imshow(colors, aspect="auto", interpolation="nearest")
    ax.set_xlim(-0.5, colors.shape[1] - 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8, labelpad=18)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_figure(
    matrix,
    metadata,
    tree,
    out_prefix,
    title,
    age_max,
):
    n = matrix.shape[0]
    values = matrix.to_numpy(dtype=float)
    masked = np.ma.array(values, mask=np.triu(np.ones_like(values, dtype=bool), k=1))

    grid = np.arange(n + 1, dtype=float)
    rr, cc = np.meshgrid(grid, grid, indexing="ij")
    x = (rr + cc) / 2.0 - 0.5
    y = (rr - cc) / 2.0

    fig_w = max(10.8, n / 34.0)
    fig_h = fig_w * 0.88
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)
    gs = fig.add_gridspec(
        nrows=6,
        ncols=1,
        height_ratios=[6.2, 0.16, 0.16, 0.16, 1.6, 1.15],
        hspace=0.045,
    )

    ax_heat = fig.add_subplot(gs[0, 0])
    ax_type = fig.add_subplot(gs[1, 0], sharex=ax_heat)
    ax_species = fig.add_subplot(gs[2, 0], sharex=ax_heat)
    ax_age = fig.add_subplot(gs[3, 0], sharex=ax_heat)
    ax_tree = fig.add_subplot(gs[4, 0], sharex=ax_heat)
    ax_leg = fig.add_subplot(gs[5, 0])

    sim_cmap = LinearSegmentedColormap.from_list(
        "intact_ltrrt_similarity",
        ["#2b83ba", "#a8cce4", "#f7f7f7", "#f4a6a6", "#d7191c"],
    )
    heat = ax_heat.pcolormesh(x, y, masked, cmap=sim_cmap, vmin=0, vmax=1, shading="flat", rasterized=True)
    ax_heat.set_aspect("equal")
    ax_heat.set_xlim(-0.5, n - 0.5)
    ax_heat.set_ylim(-0.5, n / 2 + 2)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    ax_heat.set_title(title, fontsize=13, pad=8)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    sim_cax = inset_axes(ax_heat, width="22%", height="4.5%", loc="upper left", borderpad=1.0)
    sim_cb = fig.colorbar(heat, cax=sim_cax, orientation="horizontal")
    sim_cb.set_label("intact LTR-RT similarity", fontsize=8, labelpad=2)
    sim_cb.ax.tick_params(labelsize=7, length=2)

    draw_annotation_bar(ax_type, as_color_row(metadata["selected_type"], TYPE_COLORS), "TE type")
    draw_annotation_bar(ax_species, as_color_row(metadata["label"], SPECIES_COLORS), "species")

    ages = pd.to_numeric(metadata["insertion_time_mya"], errors="coerce").fillna(0.0).to_numpy()[None, :]
    vmax = age_max if age_max is not None else max(1.0, float(np.nanmax(ages)))
    age_cmap = LinearSegmentedColormap.from_list("ltr_insert_time", ["#fff7df", "#f6c16b", "#de5b3d", "#9d1d20"])
    ax_age.imshow(ages, aspect="auto", interpolation="nearest", cmap=age_cmap, norm=Normalize(0, vmax))
    ax_age.set_xlim(-0.5, n - 0.5)
    ax_age.set_xticks([])
    ax_age.set_yticks([])
    ax_age.set_ylabel("LTR age", rotation=0, ha="right", va="center", fontsize=8, labelpad=18)
    for spine in ax_age.spines.values():
        spine.set_visible(False)

    age_cax = inset_axes(ax_heat, width="22%", height="4.5%", loc="upper left", borderpad=1.0,
                         bbox_to_anchor=(0.0, -0.11, 1, 1), bbox_transform=ax_heat.transAxes)
    age_sm = plt.cm.ScalarMappable(norm=Normalize(0, vmax), cmap=age_cmap)
    age_cb = fig.colorbar(age_sm, cax=age_cax, orientation="horizontal")
    age_cb.set_label("LTR insertion time (Mya)", fontsize=8, labelpad=2)
    age_cb.ax.tick_params(labelsize=7, length=2)

    max_depth = assign_tree_coordinates(tree, list(matrix.index))
    draw_tree(ax_tree, tree, max_depth, n)

    ax_leg.axis("off")
    type_handles = [
        Patch(color=TYPE_COLORS[t], label=f"{t} (n={(metadata['selected_type'] == t).sum()})")
        for t in TYPE_ORDER
        if (metadata["selected_type"] == t).any()
    ]
    species_handles = [
        Patch(color=SPECIES_COLORS[s], label=f"{s} (n={(metadata['label'] == s).sum()})")
        for s in SPECIES_ORDER
        if (metadata["label"] == s).any()
    ]
    ax_leg.text(0.5, 0.88, "TE clade", ha="center", va="center", fontsize=8)
    leg1 = ax_leg.legend(handles=type_handles, loc="upper center", ncol=min(5, len(type_handles)),
                         fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.80))
    ax_leg.text(0.5, 0.34, "Species", ha="center", va="center", fontsize=8)
    leg2 = ax_leg.legend(handles=species_handles, loc="upper center", ncol=min(5, len(species_handles)),
                         fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.26))
    ax_leg.add_artist(leg1)

    fig.savefig(Path(str(out_prefix) + ".png"), bbox_inches="tight")
    fig.savefig(Path(str(out_prefix) + ".pdf"), bbox_inches="tight")
    fig.savefig(Path(str(out_prefix) + ".svg"), bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--treefile", required=True, help="IQ-TREE Newick treefile.")
    parser.add_argument("--matrix", required=True, help="Square similarity matrix TSV.")
    parser.add_argument("--metadata", required=True, help="Metadata TSV containing seq_id, label, selected_type, insertion_time_mya.")
    parser.add_argument("--out-prefix", required=True, help="Output path prefix without extension.")
    parser.add_argument("--title", default="Full-length intact LTR-RT similarity ordered by IQ-TREE")
    parser.add_argument("--age-max", type=float, default=None, help="Maximum LTR insertion time for the age color bar.")
    parser.add_argument("--root-type-prefix", default=None,
                        help="Optional selected_type prefix used to reroot/order the tree, e.g. Copia.")
    parser.add_argument("--root-type-first", action="store_true",
                        help="After rerooting, rotate internal nodes so the requested root type is placed first where possible.")
    args = parser.parse_args()

    tree = parse_newick(Path(args.treefile).read_text())
    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    metadata = pd.read_csv(args.metadata, sep="\t")
    metadata["seq_id"] = metadata["seq_id"].astype(str)
    metadata["label"] = metadata["label"].astype(str)
    metadata["selected_type"] = metadata["selected_type"].astype(str)

    root_info = None
    if args.root_type_prefix:
        target = set(metadata.loc[metadata["selected_type"].str.startswith(args.root_type_prefix), "seq_id"])
        tree, root_info = reroot_by_leaf_set(tree, target)
        if args.root_type_first:
            rotate_target_clade_first(tree, target)

    tree_order = leaves_in_order(tree)

    order = [leaf for leaf in tree_order if leaf in matrix.index]
    if not order:
        raise SystemExit("No tree leaves matched the similarity matrix IDs.")
    missing = sorted(set(matrix.index) - set(order))

    matrix = matrix.loc[order, order]
    metadata = metadata.set_index("seq_id").loc[order].reset_index()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    matrix.to_csv(Path(str(out_prefix) + ".matrix.tsv"), sep="\t", float_format="%.4f")
    metadata.to_csv(Path(str(out_prefix) + ".metadata.tsv"), sep="\t", index=False)
    pd.Series(missing, name="matrix_ids_missing_from_tree").to_csv(Path(str(out_prefix) + ".missing.tsv"), sep="\t", index=False)

    draw_figure(matrix, metadata, tree, out_prefix, args.title, args.age_max)

    print(f"tree_leaves={len(tree_order)}")
    print(f"ordered={len(order)}")
    print(f"missing_from_tree={len(missing)}")
    if root_info:
        print(
            "root_type_prefix={prefix} exact={exact} target_n={target_n} root_side_n={root_side_n} other_side_n={other_side_n}".format(
                prefix=args.root_type_prefix,
                **root_info,
            )
        )
    print(f"output_prefix={out_prefix}")
    print(metadata.groupby(["label", "selected_type"]).size().reset_index(name="count").to_string(index=False))


if __name__ == "__main__":
    main()
