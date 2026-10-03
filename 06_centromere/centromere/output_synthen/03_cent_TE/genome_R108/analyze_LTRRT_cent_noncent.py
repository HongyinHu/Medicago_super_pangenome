#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


LOC_RE = re.compile(r"^([^:]+):(\d+)\.\.(\d+)$")


def parse_ltr_loc(value):
    match = LOC_RE.match(str(value))
    if not match:
        raise ValueError(f"Cannot parse LTR_loc: {value}")
    chrom = match.group(1)
    a = int(match.group(2))
    b = int(match.group(3))
    start_1based = min(a, b)
    end_1based = max(a, b)
    return chrom, start_1based - 1, end_1based


def read_bed(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr", "start", "end"],
        usecols=[0, 1, 2],
        comment="#",
    )


def read_pass_list(path, mutation_rate):
    rows = []
    with open(path) as handle:
        header = None
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split()
                continue
            fields = line.split()
            if header is None:
                raise ValueError("No header found in pass.list")
            record = dict(zip(header, fields))
            chrom, start0, end = parse_ltr_loc(record["LTR_loc"])
            identity = float(record["Identity"])
            insertion_years = record.get("Insertion_Time", "NA")
            if insertion_years in {"NA", ".", ""}:
                insertion_years = (1 - identity) / (2 * mutation_rate)
            else:
                insertion_years = float(insertion_years)
            rows.append(
                {
                    "id": record["LTR_loc"],
                    "chr": chrom,
                    "start": start0,
                    "end": end,
                    "length": end - start0,
                    "category": record.get("Category", ""),
                    "motif": record.get("Motif", ""),
                    "identity": identity,
                    "identity_pct": identity * 100,
                    "strand": record.get("Strand", ""),
                    "superfamily": record.get("SuperFamily", ""),
                    "te_type": record.get("TE_type", ""),
                    "insertion_years": insertion_years,
                    "insertion_mya": insertion_years / 1_000_000,
                }
            )
    return pd.DataFrame(rows)


def add_region(ltrs, cen_bed):
    ltrs = ltrs.copy()
    ltrs["region"] = "Non-cen"
    for chrom, regions in cen_bed.groupby("chr"):
        idx = ltrs["chr"] == chrom
        if not idx.any():
            continue
        overlap = pd.Series(False, index=ltrs.index[idx])
        for _, region in regions.iterrows():
            overlap = overlap | (
                (ltrs.loc[idx, "start"] < int(region["end"]))
                & (ltrs.loc[idx, "end"] > int(region["start"]))
            )
        ltrs.loc[overlap.index[overlap], "region"] = "Cen"
    return ltrs


def make_frequency_table(ltrs, bin_width):
    bins = np.arange(80, 100 + bin_width, bin_width)
    labels = []
    rows = []
    for left, right in zip(bins[:-1], bins[1:]):
        labels.append((left, right))
    for region in ["Non-cen", "Cen"]:
        vals = ltrs.loc[ltrs["region"] == region, "identity_pct"].dropna()
        counts, edges = np.histogram(vals, bins=bins)
        for i, count in enumerate(counts):
            rows.append(
                {
                    "region": region,
                    "bin_start": edges[i],
                    "bin_end": edges[i + 1],
                    "bin_mid": (edges[i] + edges[i + 1]) / 2,
                    "count": int(count),
                }
            )
    return pd.DataFrame(rows)


def plot_identity_hist(ltrs, out_prefix, species, bin_width):
    colors = {"Non-cen": "#ff6a00", "Cen": "#ff6a00"}
    bins = np.arange(90, 100 + bin_width, bin_width)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharex=True)
    for ax, region in zip(axes, ["Non-cen", "Cen"]):
        vals = ltrs.loc[ltrs["region"] == region, "identity_pct"].dropna()
        ax.hist(vals, bins=bins, color=colors[region], edgecolor="none")
        ax.axvline(99.5, color="red", linestyle=":", linewidth=1)
        ax.set_title(f"{species}_{region}_LTR-RTs", fontsize=11)
        ax.text(0.08, 0.78, f"n={len(vals)}", transform=ax.transAxes, fontsize=11)
        ax.set_xlabel("LTR % Identity")
        ax.set_ylabel("FREQ")
        ax.set_xlim(90, 100)
        ax.tick_params(labelsize=9)
    plt.tight_layout()
    fig.savefig(f"{out_prefix}.identity_hist.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.identity_hist.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_insertion_box(ltrs, out_prefix, species):
    data = [
        ltrs.loc[ltrs["region"] == "Non-cen", "insertion_mya"].dropna(),
        ltrs.loc[ltrs["region"] == "Cen", "insertion_mya"].dropna(),
    ]
    colors = ["#ff6a00", "#ff6a00"]
    fig, ax = plt.subplots(figsize=(3.0, 3.2))
    box = ax.boxplot(
        data,
        labels=["Non-cen", "Cen"],
        patch_artist=True,
        showfliers=False,
        widths=0.55,
    )
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor("white")
        patch.set_edgecolor(color)
        patch.set_linewidth(1.2)
    for key in ["whiskers", "caps", "medians"]:
        for item in box[key]:
            item.set_color("#ff6a00")
            item.set_linewidth(1.2)
    ax.set_ylabel("Insertion time\n(million years ago)")
    ax.set_title(f"{species}_LTR-RTs")
    ax.set_ylim(0, max(0.2, np.nanpercentile(ltrs["insertion_mya"], 98) * 1.15))
    ax.tick_params(axis="x", rotation=35)
    plt.tight_layout()
    fig.savefig(f"{out_prefix}.insertion_time_boxplot.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.insertion_time_boxplot.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def write_summary(ltrs, out_path):
    rows = []
    for region, sub in ltrs.groupby("region"):
        rows.append(
            {
                "region": region,
                "n": len(sub),
                "identity_pct_median": sub["identity_pct"].median(),
                "identity_pct_mean": sub["identity_pct"].mean(),
                "insertion_mya_median": sub["insertion_mya"].median(),
                "insertion_mya_mean": sub["insertion_mya"].mean(),
                "insertion_mya_min": sub["insertion_mya"].min(),
                "insertion_mya_max": sub["insertion_mya"].max(),
            }
        )
    pd.DataFrame(rows).sort_values("region").to_csv(out_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pass_list", required=True)
    parser.add_argument("--cent_bed", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--species", default="R108")
    parser.add_argument("--mutation_rate", type=float, default=1.3e-8)
    parser.add_argument("--bin_width", type=float, default=0.5)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_prefix = outdir / f"genome_{args.species}.LTRRT"

    ltrs = read_pass_list(args.pass_list, args.mutation_rate)
    cen_bed = read_bed(args.cent_bed)
    ltrs = add_region(ltrs, cen_bed)

    ltrs.to_csv(f"{out_prefix}.cent_noncent.tsv", sep="\t", index=False)
    freq = make_frequency_table(ltrs, args.bin_width)
    freq.to_csv(f"{out_prefix}.identity_frequency.tsv", sep="\t", index=False)
    write_summary(ltrs, f"{out_prefix}.summary.tsv")
    plot_identity_hist(ltrs, str(out_prefix), args.species, args.bin_width)
    plot_insertion_box(ltrs, str(out_prefix), args.species)

    print(f"Wrote {len(ltrs)} intact LTR-RT records")
    print(ltrs["region"].value_counts().to_string())
    print(f"Output prefix: {out_prefix}")


if __name__ == "__main__":
    main()
