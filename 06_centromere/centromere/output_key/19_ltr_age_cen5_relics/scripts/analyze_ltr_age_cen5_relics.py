#!/usr/bin/env python3
import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})


REGIONS = [
    {
        "region_id": "Chr3_CEN5_relic",
        "group": "Chr3 relic",
        "chrom": "Chr3",
        "start": 20591033,
        "end": 22176458,
        "class": "CEN5_relic",
    },
    {
        "region_id": "Chr3_CEN5_to_CEN6_transition",
        "group": "Chr3 transition",
        "chrom": "Chr3",
        "start": 22119602,
        "end": 22565328,
        "class": "CEN5_to_CEN6_transition",
    },
    {
        "region_id": "Chr5_CEN5_relic",
        "group": "Chr5 relic",
        "chrom": "Chr5",
        "start": 25293760,
        "end": 25645386,
        "class": "CEN5_relic",
    },
    {
        "region_id": "Chr5_CEN5_to_CEN6_transition",
        "group": "Chr5 transition",
        "chrom": "Chr5",
        "start": 25489881,
        "end": 26760313,
        "class": "CEN5_to_CEN6_transition",
    },
    {
        "region_id": "Chr3_relic_plus_transition",
        "group": "Chr3 relic+transition",
        "chrom": "Chr3",
        "start": 20591033,
        "end": 22565328,
        "class": "merged_relic_transition",
    },
    {
        "region_id": "Chr5_relic_plus_transition",
        "group": "Chr5 relic+transition",
        "chrom": "Chr5",
        "start": 25293760,
        "end": 26760313,
        "class": "merged_relic_transition",
    },
]


COLORS = {
    "Chr3 relic": "#2f8f83",
    "Chr3 transition": "#7bb9a9",
    "Chr5 relic": "#d99058",
    "Chr5 transition": "#e5b07d",
    "Chr3 relic+transition": "#28786f",
    "Chr5 relic+transition": "#c9895b",
}


def parse_attrs(attr):
    out = {}
    for item in attr.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def jc_distance(d):
    # Jukes-Cantor correction. d is raw divergence between two LTRs.
    if d <= 0:
        return 0.0
    if d >= 0.74:
        return float("nan")
    return -0.75 * math.log(1 - 4 * d / 3)


def age_mya_from_identity(identity, rate):
    d = 1.0 - identity
    k = jc_distance(d)
    if math.isnan(k):
        return float("nan")
    return k / (2 * rate) / 1_000_000


def read_ltrs(gff, keep_chroms, substitution_rate):
    ltrs = []
    with Path(gff).open(errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in keep_chroms:
                continue
            attrs = parse_attrs(p[8])
            if "ltr_identity" not in attrs:
                continue
            # EDTA writes subfeatures (LTRs, TSDs and internal retrotransposon
            # pieces) with the same ltr_identity. Count only the parent intact
            # LTR record to avoid inflating n by subfeature number.
            if p[2] != "repeat_region":
                continue
            try:
                ltr_identity = float(attrs["ltr_identity"])
            except ValueError:
                continue
            if not (0 < ltr_identity <= 1):
                continue
            chrom = p[0]
            start = int(p[3]) - 1
            end = int(p[4])
            classification = attrs.get("classification", p[2])
            name = attrs.get("Name", attrs.get("ID", "unknown"))
            ltrs.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "feature": p[2],
                "name": name,
                "classification": classification,
                "ltr_identity": ltr_identity,
                "raw_divergence": 1.0 - ltr_identity,
                "jc_distance": jc_distance(1.0 - ltr_identity),
                "age_mya": age_mya_from_identity(ltr_identity, substitution_rate),
                "method": attrs.get("method", "NA"),
            })
    return ltrs


def collect_region_ltrs(ltrs):
    rows = []
    for region in REGIONS:
        for ltr in ltrs:
            if ltr["chrom"] != region["chrom"]:
                continue
            ov = overlap(region["start"], region["end"], ltr["start"], ltr["end"])
            if ov <= 0:
                continue
            row = dict(region)
            row.update({
                "region_length_bp": region["end"] - region["start"],
                "ltr_chrom": ltr["chrom"],
                "ltr_start": ltr["start"],
                "ltr_end": ltr["end"],
                "ltr_length_bp": ltr["end"] - ltr["start"],
                "overlap_bp": ov,
                "overlap_fraction_of_ltr": ov / max(1, ltr["end"] - ltr["start"]),
                "name": ltr["name"],
                "classification": ltr["classification"],
                "feature": ltr["feature"],
                "method": ltr["method"],
                "ltr_identity": ltr["ltr_identity"],
                "raw_divergence": ltr["raw_divergence"],
                "jc_distance": ltr["jc_distance"],
                "age_mya": ltr["age_mya"],
            })
            rows.append(row)
    return rows


def summarize(rows):
    by_region = defaultdict(list)
    for row in rows:
        if row["overlap_fraction_of_ltr"] < 0.10:
            continue
        if math.isnan(row["age_mya"]):
            continue
        by_region[row["region_id"]].append(row)

    summary = []
    for region in REGIONS:
        vals = by_region.get(region["region_id"], [])
        ages = [r["age_mya"] for r in vals]
        identities = [r["ltr_identity"] for r in vals]
        fams = Counter(r["name"] for r in vals)
        classes = Counter(r["classification"] for r in vals)
        if ages:
            row = {
                "region_id": region["region_id"],
                "group": region["group"],
                "class": region["class"],
                "chrom": region["chrom"],
                "start": region["start"],
                "end": region["end"],
                "length_bp": region["end"] - region["start"],
                "intact_LTR_count": len(ages),
                "median_age_mya": statistics.median(ages),
                "mean_age_mya": statistics.mean(ages),
                "min_age_mya": min(ages),
                "max_age_mya": max(ages),
                "median_ltr_identity": statistics.median(identities),
                "top_families": ",".join(f"{k}:{v}" for k, v in fams.most_common(5)),
                "top_classes": ",".join(f"{k}:{v}" for k, v in classes.most_common(5)),
            }
        else:
            row = {
                "region_id": region["region_id"],
                "group": region["group"],
                "class": region["class"],
                "chrom": region["chrom"],
                "start": region["start"],
                "end": region["end"],
                "length_bp": region["end"] - region["start"],
                "intact_LTR_count": 0,
                "median_age_mya": "NA",
                "mean_age_mya": "NA",
                "min_age_mya": "NA",
                "max_age_mya": "NA",
                "median_ltr_identity": "NA",
                "top_families": "NA",
                "top_classes": "NA",
            }
        summary.append(row)
    return summary


def compare_chr3_chr5(summary):
    by_id = {r["region_id"]: r for r in summary}
    pairs = [
        ("relic", "Chr3_CEN5_relic", "Chr5_CEN5_relic"),
        ("transition", "Chr3_CEN5_to_CEN6_transition", "Chr5_CEN5_to_CEN6_transition"),
        ("relic_plus_transition", "Chr3_relic_plus_transition", "Chr5_relic_plus_transition"),
    ]
    out = []
    for label, a, b in pairs:
        ra, rb = by_id[a], by_id[b]
        def num(x):
            return None if x == "NA" else float(x)
        ma, mb = num(ra["median_age_mya"]), num(rb["median_age_mya"])
        out.append({
            "comparison": label,
            "chr3_region": a,
            "chr5_region": b,
            "chr3_LTR_count": ra["intact_LTR_count"],
            "chr5_LTR_count": rb["intact_LTR_count"],
            "chr3_median_age_mya": ra["median_age_mya"],
            "chr5_median_age_mya": rb["median_age_mya"],
            "median_age_difference_mya": "NA" if ma is None or mb is None else abs(ma - mb),
            "interpretation": (
                "comparable_age_distribution_candidate"
                if ma is not None and mb is not None and abs(ma - mb) <= 0.5
                else "different_or_insufficient"
            ),
        })
    return out


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_ages(rows, summary, outdir):
    filtered = [r for r in rows if r["overlap_fraction_of_ltr"] >= 0.10 and not math.isnan(r["age_mya"])]
    groups = [r["group"] for r in REGIONS]
    ages_by_group = {g: [r["age_mya"] for r in filtered if r["group"] == g] for g in groups}

    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.2), constrained_layout=True)
    ax = axes[0]
    x = np.arange(len(groups))
    data = [ages_by_group[g] for g in groups]
    parts = ax.violinplot([d if d else [np.nan] for d in data], positions=x, widths=0.72,
                          showmeans=False, showmedians=False, showextrema=False)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(COLORS.get(groups[i], "#999999"))
        body.set_edgecolor("none")
        body.set_alpha(0.45 if data[i] else 0.08)
    for i, g in enumerate(groups):
        vals = ages_by_group[g]
        if vals:
            jitter = np.linspace(-0.16, 0.16, len(vals)) if len(vals) > 1 else [0]
            ax.scatter(np.array([i] * len(vals)) + np.array(jitter), vals,
                       s=14, color=COLORS.get(g, "#666666"), edgecolor="white", linewidth=0.3, zorder=3)
            ax.plot([i - 0.22, i + 0.22], [statistics.median(vals)] * 2,
                    color="#222222", lw=1.2, zorder=4)
            ax.text(i, max(vals) + 0.08, f"n={len(vals)}", ha="center", va="bottom", fontsize=6)
        else:
            ax.text(i, 0.08, "n=0", ha="center", va="bottom", fontsize=6, color="#777777")
    label_map = {
        "Chr3 relic": "Chr3\nrelic",
        "Chr3 transition": "Chr3\ntransition",
        "Chr5 relic": "Chr5\nrelic",
        "Chr5 transition": "Chr5\ntransition",
        "Chr3 relic+transition": "Chr3\nrelic+trans.",
        "Chr5 relic+transition": "Chr5\nrelic+trans.",
    }
    ax.set_xticks(x)
    ax.set_xticklabels([label_map.get(g, g) for g in groups], rotation=0)
    ax.set_ylabel("LTR insertion age (Mya)")
    ax.set_title("LTR age around CEN5 relic and transition intervals", loc="left", fontweight="bold")

    ax2 = axes[1]
    merged_groups = ["Chr3 relic+transition", "Chr5 relic+transition"]
    bins = np.linspace(0, max([r["age_mya"] for r in filtered] + [1]) + 0.5, 12)
    for g in merged_groups:
        vals = ages_by_group[g]
        if vals:
            ax2.hist(vals, bins=bins, alpha=0.55, label=f"{g} (n={len(vals)})",
                     color=COLORS.get(g, "#666666"))
    ax2.set_xlabel("LTR insertion age (Mya)")
    ax2.set_ylabel("Count")
    ax2.set_title("Merged relic+transition windows", loc="left", fontweight="bold")
    ax2.legend(fontsize=6)

    prefix = outdir / "Mpo_CEN5_relic_transition_LTR_insertion_age"
    fig.savefig(str(prefix) + ".png", dpi=450, bbox_inches="tight")
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")


def write_report(path, summary, comparison, substitution_rate):
    by_id = {r["region_id"]: r for r in summary}
    lines = [
        "# Mpo CEN5 relic / transition LTR insertion-age analysis",
        "",
        f"Substitution rate used for approximate absolute age: `{substitution_rate:g}` substitutions/site/year. Ages are best interpreted comparatively because the exact Medicago LTR substitution rate may vary.",
        "",
        "## Region-level summary",
        "",
    ]
    for r in summary:
        lines.append(
            f"- {r['group']} ({r['chrom']}:{r['start']}-{r['end']}): "
            f"intact LTR n={r['intact_LTR_count']}, median age={r['median_age_mya']}, "
            f"top classes={r['top_classes']}."
        )
    lines.extend(["", "## Chr3 versus Chr5 comparison", ""])
    for c in comparison:
        lines.append(
            f"- {c['comparison']}: Chr3 n={c['chr3_LTR_count']}, median={c['chr3_median_age_mya']}; "
            f"Chr5 n={c['chr5_LTR_count']}, median={c['chr5_median_age_mya']}; "
            f"difference={c['median_age_difference_mya']}; interpretation={c['interpretation']}."
        )
    lines.extend([
        "",
        "## Writing boundary",
        "",
        "Similar LTR age distributions would support temporally comparable TE activity or repeat turnover around the two CEN5 relic/transition regions. This alone does not prove that TE insertion initiated the rearrangement or that both relics were generated in the same single event.",
    ])
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--rate", type=float, default=1.3e-8)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    keep_chroms = {r["chrom"] for r in REGIONS}
    ltrs = read_ltrs(args.gff, keep_chroms, args.rate)
    rows = collect_region_ltrs(ltrs)
    summary = summarize(rows)
    comparison = compare_chr3_chr5(summary)

    detail_fields = [
        "region_id", "group", "class", "chrom", "start", "end", "region_length_bp",
        "ltr_chrom", "ltr_start", "ltr_end", "ltr_length_bp", "overlap_bp",
        "overlap_fraction_of_ltr", "name", "classification", "feature", "method",
        "ltr_identity", "raw_divergence", "jc_distance", "age_mya",
    ]
    summary_fields = [
        "region_id", "group", "class", "chrom", "start", "end", "length_bp",
        "intact_LTR_count", "median_age_mya", "mean_age_mya", "min_age_mya",
        "max_age_mya", "median_ltr_identity", "top_families", "top_classes",
    ]
    write_tsv(outdir / "Mpo_CEN5_relic_transition_intact_LTR_detail.tsv", rows, detail_fields)
    write_tsv(outdir / "Mpo_CEN5_relic_transition_LTR_age_summary.tsv", summary, summary_fields)
    write_tsv(outdir / "Mpo_CEN5_relic_transition_LTR_age_chr3_chr5_comparison.tsv", comparison)
    plot_ages(rows, summary, outdir)
    write_report(outdir / "Mpo_CEN5_relic_transition_LTR_age_summary_CN.md", summary, comparison, args.rate)


if __name__ == "__main__":
    main()
