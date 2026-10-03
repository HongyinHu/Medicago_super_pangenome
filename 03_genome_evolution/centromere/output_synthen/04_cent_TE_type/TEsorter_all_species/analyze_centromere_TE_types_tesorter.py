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
                    "name": attrs.get("Name", attrs.get("ID", "NA")),
                    "classification": attrs.get("classification", "unknown"),
                    "method": attrs.get("method", "NA"),
                }
            )
    return pd.DataFrame(rows)


def read_tesorter_cls(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fields = line.rstrip("\n").split(None, 6)
            if len(fields) < 4:
                continue
            while len(fields) < 7:
                fields.append("")
            rows.append(fields[:7])
    if not rows:
        raise ValueError(f"No TEsorter classification rows parsed: {path}")
    cls = pd.DataFrame(rows, columns=["te_raw", "ts_order", "ts_superfamily", "ts_clade", "complete", "ts_strand", "domains"])
    cls["te_name"] = cls["te_raw"].astype(str).str.split("#", n=1).str[0]
    cls["ts_clade"] = cls["ts_clade"].fillna("unclassified").replace({"": "unclassified", "-": "unclassified"})
    cls["ts_superfamily"] = cls["ts_superfamily"].fillna("unknown").replace({"": "unknown", "-": "unknown"})
    cls = cls.drop_duplicates("te_name", keep="first")
    return cls.set_index("te_name")


def split_classification(value):
    if not isinstance(value, str) or value == "":
        return "unknown", "unknown"
    if "/" in value:
        order, superfamily = value.split("/", 1)
    else:
        order, superfamily = value, value
    return order, superfamily


def fine_label(row):
    order, superfamily = split_classification(row["classification"])
    clade = row.get("ts_clade", np.nan)
    ts_superfamily = row.get("ts_superfamily", np.nan)
    if order == "LTR":
        if isinstance(clade, str) and clade not in {"", "-", "unclassified", "unknown", "nan"}:
            return ts_superfamily if isinstance(ts_superfamily, str) and ts_superfamily else superfamily, clade
        if superfamily in {"Copia", "Gypsy"}:
            return superfamily, f"{superfamily}_unclassified"
        return "LTR", "LTR_unclassified"
    if order == "DNA":
        return "DNA", superfamily
    if order in {"LINE", "SINE", "MITE"}:
        return order, superfamily
    return order, row["classification"]


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
    group_cols = [c for c in te.columns]
    out = out.groupby(group_cols, dropna=False, as_index=False)["overlap_bp"].sum()
    return out


def attach_tesorter(te, cls):
    te = te.copy()
    te["tesorter_key"] = te["name"]
    exact = te["tesorter_key"].isin(cls.index)
    ltr_alt = te["name"].astype(str).str.replace("_LTR$", "_INT", regex=True)
    use_alt = (~exact) & ltr_alt.isin(cls.index)
    te.loc[use_alt, "tesorter_key"] = ltr_alt[use_alt]
    merged = te.join(cls[["ts_order", "ts_superfamily", "ts_clade", "complete", "domains"]], on="tesorter_key")
    labels = merged.apply(fine_label, axis=1)
    merged["plot_group"] = [x[0] for x in labels]
    merged["fine_type"] = [x[1] for x in labels]
    return merged


def summarize(hits, genome):
    if hits.empty:
        return pd.DataFrame()
    hits = hits.copy()
    hits["te_length"] = hits["end"] - hits["start"]
    meta = (
        hits.groupby(["plot_group", "fine_type"], dropna=False)
        .agg(
            edta_classifications=("classification", lambda x: ";".join(sorted(set(str(v) for v in x if str(v) != "nan")))),
            tesorter_superfamilies=("ts_superfamily", lambda x: ";".join(sorted(set(str(v) for v in x if str(v) not in {"nan", ""})))),
            tesorter_clades=("ts_clade", lambda x: ";".join(sorted(set(str(v) for v in x if str(v) not in {"nan", ""})))),
        )
        .reset_index()
    )
    summary = (
        hits.groupby(["plot_group", "fine_type"], dropna=False)
        .agg(copy_number=("id", "count"), overlap_bp=("overlap_bp", "sum"), mean_te_length=("te_length", "mean"))
        .reset_index()
    )
    summary = summary.merge(meta, on=["plot_group", "fine_type"], how="left")
    summary["genome"] = genome
    return summary[
        [
            "genome",
            "plot_group",
            "fine_type",
            "copy_number",
            "overlap_bp",
            "mean_te_length",
            "edta_classifications",
            "tesorter_superfamilies",
            "tesorter_clades",
        ]
    ]


def plot_summary(summary, out_prefix, title="Functional centromere"):
    if summary.empty:
        return
    group_rank = {"Copia": 0, "Gypsy": 1, "LTR": 2, "LINE": 3, "SINE": 4, "DNA": 5, "MITE": 6, "unknown": 9}
    labels = (
        summary[["plot_group", "fine_type"]]
        .drop_duplicates()
        .assign(rank=lambda d: d["plot_group"].map(group_rank).fillna(8))
        .sort_values(["rank", "plot_group", "fine_type"])
    )
    label_order = labels["fine_type"].tolist()
    display_labels = [
        label.replace("_unclassified", "\nunc.").replace("_", "\n")
        for label in label_order
    ]
    x = np.arange(len(label_order))
    vals = [int(summary[summary["fine_type"] == label]["copy_number"].sum()) for label in label_order]

    fig, ax = plt.subplots(figsize=(max(8, len(label_order) * 0.58), 4.8))
    ax.bar(x, vals, width=0.72, color="#ff6a00", edgecolor="white", linewidth=0.2)
    ax.set_ylabel("Copy number")
    ax.set_xticks(x)
    ax.set_xticklabels(display_labels, rotation=90)
    ax.text(0.01, 0.98, title, ha="left", va="top", transform=ax.transAxes, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    trans = ax.get_xaxis_transform()
    for group, sublabels in labels.groupby("plot_group", sort=False):
        idxs = [label_order.index(v) for v in sublabels["fine_type"]]
        mid = (min(idxs) + max(idxs)) / 2
        ax.text(mid, -0.36, group, ha="center", va="top", transform=trans, fontsize=10)
        ax.plot([min(idxs) - 0.45, max(idxs) + 0.45], [-0.29, -0.29], color="black", linewidth=0.6, transform=trans, clip_on=False)

    fig.subplots_adjust(bottom=0.43, left=0.08, right=0.99, top=0.96)
    fig.savefig(f"{out_prefix}.TEsorter_clade_copy_number.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.TEsorter_clade_copy_number.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--teanno_gff", required=True)
    parser.add_argument("--centromere_bed", required=True)
    parser.add_argument("--tesorter_cls", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--chrom_regex", default=r"^Chr[0-9]+$")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    te = read_teanno_gff(args.teanno_gff, args.chrom_regex)
    cls = read_tesorter_cls(args.tesorter_cls)
    te = attach_tesorter(te, cls)
    cen = read_bed(args.centromere_bed, args.chrom_regex)
    hits = intersect_te(te, cen)
    summary = summarize(hits, args.genome)

    hits.to_csv(outdir / f"{args.genome}.functional_centromere.TE_hits.TEsorter_clade.tsv", sep="\t", index=False)
    summary.to_csv(outdir / f"{args.genome}.functional_centromere.TE_type_summary.TEsorter_clade.tsv", sep="\t", index=False)
    plot_summary(summary, str(outdir / f"{args.genome}.functional_centromere.TE_type"))
    print(f"{args.genome}: {len(hits)} TE records overlap functional centromeres")
    print(f"{args.genome}: {hits['ts_clade'].notna().sum()} records have direct TEsorter clade mapping")
    print(f"Output directory: {outdir}")


if __name__ == "__main__":
    main()
