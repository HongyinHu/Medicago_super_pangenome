#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SKIP_FEATURES = {
    "repeat_region",
    "target_site_duplication",
    "long_terminal_repeat",
    "terminal_inverted_repeat",
    "direct_repeat",
}


def parse_attrs(value):
    attrs = {}
    for item in str(value).split(";"):
        if "=" not in item:
            continue
        key, val = item.split("=", 1)
        attrs[key] = val
    return attrs


def read_bed(path, chrom_regex):
    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end"],
        usecols=[0, 1, 2],
        comment="#",
    )
    if chrom_regex:
        bed = bed[bed["chr"].astype(str).str.match(chrom_regex)].copy()
    return bed


def read_teanno_gff(path, chrom_regex):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr_text = fields
            if chrom_regex and not re.match(chrom_regex, chrom):
                continue
            if feature in SKIP_FEATURES:
                continue
            attrs = parse_attrs(attr_text)
            classification = attrs.get("classification", "unknown")
            name = attrs.get("Name", attrs.get("ID", "NA"))
            rows.append(
                {
                    "chr": chrom,
                    "start": int(start) - 1,
                    "end": int(end),
                    "source": source,
                    "feature": feature,
                    "score": score,
                    "strand": strand,
                    "id": attrs.get("ID", "NA"),
                    "name": name,
                    "classification": classification,
                    "method": attrs.get("method", "NA"),
                }
            )
    return pd.DataFrame(rows)


def split_classification(value):
    if not isinstance(value, str) or value == "":
        return "unknown", "unknown", "unknown"
    if "/" in value:
        order, superfamily = value.split("/", 1)
    else:
        order, superfamily = value, value
    if order == "LTR":
        label = superfamily
    elif order in {"LINE", "SINE"}:
        label = superfamily if superfamily != "unknown" else order
    elif order == "DNA":
        label = superfamily
    else:
        label = value
    return order, superfamily, label


def intersect_te(te, cen_bed):
    hits = []
    for chrom, regions in cen_bed.groupby("chr"):
        sub = te[te["chr"] == chrom]
        if sub.empty:
            continue
        for _, region in regions.iterrows():
            start = int(region["start"])
            end = int(region["end"])
            ov = sub[(sub["start"] < end) & (sub["end"] > start)].copy()
            if ov.empty:
                continue
            ov["cen_start"] = start
            ov["cen_end"] = end
            ov["overlap_start"] = np.maximum(ov["start"].to_numpy(), start)
            ov["overlap_end"] = np.minimum(ov["end"].to_numpy(), end)
            ov["overlap_bp"] = ov["overlap_end"] - ov["overlap_start"]
            hits.append(ov)
    if not hits:
        return pd.DataFrame(columns=list(te.columns) + ["cen_start", "cen_end", "overlap_start", "overlap_end", "overlap_bp"])
    out = pd.concat(hits, ignore_index=True)
    # A TE can overlap only one functional centromere interval per chromosome in these data,
    # but de-duplicate defensively by ID/coordinates and keep summed overlap.
    group_cols = [c for c in te.columns]
    out = out.groupby(group_cols, dropna=False, as_index=False)["overlap_bp"].sum()
    return out


def summarize(hits, genome):
    if hits.empty:
        return pd.DataFrame()
    parts = hits["classification"].apply(split_classification)
    hits = hits.copy()
    hits["te_order"] = [p[0] for p in parts]
    hits["te_superfamily"] = [p[1] for p in parts]
    hits["plot_label"] = [p[2] for p in parts]
    hits["te_length"] = hits["end"] - hits["start"]
    summary = (
        hits.groupby(["te_order", "te_superfamily", "classification", "plot_label"], dropna=False)
        .agg(copy_number=("id", "count"), overlap_bp=("overlap_bp", "sum"), mean_te_length=("te_length", "mean"))
        .reset_index()
    )
    summary["genome"] = genome
    return summary[["genome", "te_order", "te_superfamily", "classification", "plot_label", "copy_number", "overlap_bp", "mean_te_length"]]


def plot_grouped(summary, out_prefix):
    if summary.empty:
        return
    preferred_genome_order = ["genome_A17", "genome_R108", "genome_Mpo", "genome_Msa", "genome_474"]
    present_genomes = set(summary["genome"].dropna().unique())
    genome_order = [genome for genome in preferred_genome_order if genome in present_genomes]
    label_map = {
        "genome_A17": "A17",
        "genome_R108": "R108",
        "genome_Mpo": "Mpo",
        "genome_Msa": "Msa",
        "genome_474": "474",
    }
    colors = {
        "genome_A17": "#008751",
        "genome_R108": "#ff6a00",
        "genome_Mpo": "#3b6fb6",
        "genome_Msa": "#8a4fb3",
        "genome_474": "#c83f49",
    }
    order_rank = {"LTR": 0, "LINE": 1, "SINE": 2, "DNA": 3, "TIR": 3, "Helitron": 3, "unknown": 9}
    labels = (
        summary[["te_order", "te_superfamily", "classification"]]
        .drop_duplicates()
        .assign(rank=lambda d: d["te_order"].map(order_rank).fillna(5))
        .sort_values(["rank", "te_order", "te_superfamily", "classification"])
    )
    label_order = labels["classification"].tolist()
    tick_labels = [value.replace("/", "\n") for value in label_order]
    x = np.arange(len(label_order))
    width = min(0.7 / max(len(genome_order), 1), 0.15)
    offsets = (np.arange(len(genome_order)) - (len(genome_order) - 1) / 2) * width

    fig, ax = plt.subplots(figsize=(max(8, len(label_order) * 0.65), 4.8))
    for i, genome in enumerate(genome_order):
        sub = summary[summary["genome"] == genome].set_index("classification")
        vals = [sub.loc[label, "copy_number"] if label in sub.index else 0 for label in label_order]
        vals = []
        for label in label_order:
            vals.append(int(summary[(summary["genome"] == genome) & (summary["classification"] == label)]["copy_number"].sum()))
        ax.bar(x + offsets[i], vals, width=width, color=colors[genome], label=label_map[genome], edgecolor="white", linewidth=0.2)

    ax.set_ylabel("Copy number")
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, rotation=90)
    ax.text(0.01, 0.98, "Functional centromere", ha="left", va="top", transform=ax.transAxes, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, ncol=min(5, len(genome_order)), loc="upper center", bbox_to_anchor=(0.5, -0.22))

    # Add broad group labels below the axis.
    trans = ax.get_xaxis_transform()
    for order, sublabels in labels.groupby("te_order", sort=False):
        idxs = [label_order.index(v) for v in sublabels["classification"]]
        if not idxs:
            continue
        mid = (min(idxs) + max(idxs)) / 2
        ax.text(mid, -0.30, order, ha="center", va="top", transform=trans, fontsize=10)
        ax.plot([min(idxs) - 0.45, max(idxs) + 0.45], [-0.24, -0.24], color="black", linewidth=0.6, transform=trans, clip_on=False)

    fig.subplots_adjust(bottom=0.36, left=0.08, right=0.99, top=0.96)
    fig.savefig(f"{out_prefix}.classification_copy_number.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.classification_copy_number.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="TSV with columns: genome,teanno_gff,centromere_bed")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--chrom_regex", default=r"^Chr[0-9]+$")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    config = pd.read_csv(args.config, sep="\t")

    all_hits = []
    all_summary = []
    for _, row in config.iterrows():
        genome = row["genome"]
        te = read_teanno_gff(row["teanno_gff"], args.chrom_regex)
        cen = read_bed(row["centromere_bed"], args.chrom_regex)
        hits = intersect_te(te, cen)
        hits["genome"] = genome
        out_hits = outdir / f"{genome}.functional_centromere.TE_hits.tsv"
        hits.to_csv(out_hits, sep="\t", index=False)
        summary = summarize(hits, genome)
        summary.to_csv(outdir / f"{genome}.functional_centromere.TE_type_summary.tsv", sep="\t", index=False)
        all_hits.append(hits)
        all_summary.append(summary)
        print(f"{genome}: {len(hits)} TE records overlap functional centromeres")

    combined_hits = pd.concat(all_hits, ignore_index=True) if all_hits else pd.DataFrame()
    combined_summary = pd.concat(all_summary, ignore_index=True) if all_summary else pd.DataFrame()
    combined_hits.to_csv(outdir / "all_species.functional_centromere.TE_hits.tsv", sep="\t", index=False)
    combined_summary.to_csv(outdir / "all_species.functional_centromere.TE_type_summary.tsv", sep="\t", index=False)
    plot_grouped(combined_summary, str(outdir / "all_species.functional_centromere.TE_type"))
    print(f"Output directory: {outdir}")


if __name__ == "__main__":
    main()
