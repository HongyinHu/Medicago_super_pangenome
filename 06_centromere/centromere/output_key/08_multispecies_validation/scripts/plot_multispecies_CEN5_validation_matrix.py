#!/usr/bin/env python3
"""Plot a compact multispecies CEN5 validation matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compact_support(text: str) -> str:
    if not text or "blocks" not in text:
        return text
    parts = text.split("/")
    try:
        bp = int(parts[0].replace("bp", ""))
        if bp >= 1_000_000:
            left = f"{bp / 1_000_000:.2f} Mb"
        elif bp >= 1000:
            left = f"{bp / 1000:.1f} kb"
        else:
            left = f"{bp} bp"
        block_part = parts[1].replace("blocks", " blocks")
        return f"{left} / {block_part}"
    except Exception:
        return text.replace("bp/", " bp / ").replace("blocks", " blocks")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", default="path/to/project/10.centromere_analysis/output_key")
    args = parser.parse_args()

    run_root = Path(args.run_root)
    out_dir = run_root / "08_multispecies_validation"
    species_rows = read_tsv(out_dir / "multispecies_CEN5_species_summary.tsv")

    species_order = ["R108", "A17", "474", "Msa", "Mpo"]
    col_labels = [
        "x=8 active\nCEN5 retained",
        "R108 CEN5\nprojection support",
        "Mpo primary\nactive counterpart",
        "CEN5-associated\nremnant in Mpo",
    ]
    state_colors = {
        "yes": "#2E8B57",
        "moderate": "#8FBC8F",
        "no": "#C43B3B",
        "not_applicable": "#D9D9D9",
        "remnant": "#D89B2B",
    }
    state_labels = {
        "yes": "supported",
        "moderate": "moderate / indirect",
        "no": "not supported",
        "not_applicable": "not applicable",
        "remnant": "local remnant",
    }

    row_map = {row["species"]: row for row in species_rows}
    matrix = []
    annotations = []
    for sp in species_order:
        row = row_map[sp]
        if sp == "Mpo":
            states = ["no", "not_applicable", "not_applicable", "remnant"]
            notes = ["7 active cores", "", "", "Chr4:18.32-18.43 Mb"]
        elif sp == "A17":
            states = ["yes", "moderate", "no", "not_applicable"]
            notes = ["x=8", "5.1 kb / 3 blocks", "absent in Mpo", ""]
        elif sp == "Msa":
            states = ["yes", "moderate", "no", "not_applicable"]
            notes = ["x=8", "pm3Mb support", "absent in Mpo", ""]
        elif sp == "R108":
            states = ["yes", "yes", "no", "not_applicable"]
            notes = ["x=8", "self CEN5", "absent in Mpo", ""]
        else:
            states = ["yes", "yes", "no", "not_applicable"]
            notes = ["x=8", compact_support(row["CEN5_evidence"]), "absent in Mpo", ""]
        matrix.append(states)
        annotations.append(notes)

    fig, ax = plt.subplots(figsize=(8.2, 3.8))
    for i, states in enumerate(matrix):
        for j, state in enumerate(states):
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=state_colors[state], edgecolor="white", linewidth=1.5))
            note = annotations[i][j]
            text = state_labels[state] if not note else note
            color = "white" if state in {"yes", "no", "remnant"} else "#222222"
            ax.text(j + 0.5, i + 0.5, text, ha="center", va="center", fontsize=7.5, color=color, wrap=True)

    ax.set_xlim(0, len(col_labels))
    ax.set_ylim(0, len(species_order))
    ax.invert_yaxis()
    ax.set_xticks([x + 0.5 for x in range(len(col_labels))])
    ax.set_xticklabels(col_labels, fontsize=9)
    ax.set_yticks([y + 0.5 for y in range(len(species_order))])
    ax.set_yticklabels(species_order, fontsize=10, fontstyle="italic")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Multispecies validation of ancestral CEN5 fate in Mpo", fontsize=11, pad=12)

    legend = [
        Patch(facecolor=state_colors["yes"], label="supported"),
        Patch(facecolor=state_colors["moderate"], label="moderate / indirect"),
        Patch(facecolor=state_colors["no"], label="not a primary active Mpo counterpart"),
        Patch(facecolor=state_colors["remnant"], label="local remnant"),
    ]
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, frameon=False, fontsize=8)
    fig.tight_layout()

    png = out_dir / "multispecies_CEN5_validation_matrix.png"
    pdf = out_dir / "multispecies_CEN5_validation_matrix.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
