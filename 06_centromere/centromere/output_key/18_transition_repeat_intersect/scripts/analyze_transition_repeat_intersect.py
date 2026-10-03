#!/usr/bin/env python3
import argparse
import csv
import random
import statistics
from collections import defaultdict
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


PAL = {
    "edta": "#4d83b5",
    "trash": "#40a486",
    "trf": "#c9895b",
    "combined": "#d95f5f",
    "random": "#b8b8b8",
    "text": "#222222",
}


TRANSITIONS = [
    {
        "interval_id": "Mpo_Chr3_CEN5_to_CEN6_transition",
        "chrom": "Chr3",
        "start": 22119602,
        "end": 22565328,
        "label": "Chr3 CEN5->CEN6",
    },
    {
        "interval_id": "Mpo_Chr5_CEN5_to_CEN6_transition",
        "chrom": "Chr5",
        "start": 25489881,
        "end": 26760313,
        "label": "Chr5 CEN5->CEN6",
    },
]


def read_fai(path):
    lengths = {}
    with Path(path).open() as handle:
        for line in handle:
            p = line.split("\t")
            if len(p) >= 2:
                lengths[p[0]] = int(p[1])
    return lengths


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def merge_intervals(intervals, gap=0):
    if not intervals:
        return []
    intervals = sorted((int(a), int(b)) for a, b in intervals if int(b) > int(a))
    merged = [intervals[0]]
    for a, b in intervals[1:]:
        la, lb = merged[-1]
        if a <= lb + gap:
            merged[-1] = (la, max(lb, b))
        else:
            merged.append((a, b))
    return merged


def covered_bp(intervals):
    return sum(b - a for a, b in merge_intervals(intervals))


def parse_attr(attr):
    out = {}
    for part in attr.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def major_te_class(classification, feature_type):
    cls = classification or feature_type or "unknown"
    if "/" in cls:
        return cls.split("/", 1)[0]
    if "LTR" in cls:
        return "LTR"
    if "TIR" in cls:
        return "TIR"
    if "LINE" in cls:
        return "LINE"
    if "Helitron" in cls:
        return "Helitron"
    return cls


def load_edta(path, keep_chroms):
    by_chrom = defaultdict(list)
    with Path(path).open(errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in keep_chroms:
                continue
            chrom = p[0]
            start = int(p[3]) - 1
            end = int(p[4])
            attr = parse_attr(p[8])
            classification = attr.get("classification", p[2])
            name = attr.get("Name", attr.get("ID", "unknown"))
            by_chrom[chrom].append({
                "start": start,
                "end": end,
                "feature": p[2],
                "classification": classification,
                "major": major_te_class(classification, p[2]),
                "name": name,
            })
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x["start"])
    return by_chrom


def load_trash_bed(path, keep_chroms):
    by_chrom = defaultdict(list)
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] not in keep_chroms:
                continue
            by_chrom[p[0]].append({
                "start": int(p[1]),
                "end": int(p[2]),
                "name": p[3] if len(p) > 3 else "TRASH_array",
                "period": p[4] if len(p) > 4 else "NA",
            })
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x["start"])
    return by_chrom


def parse_trf_dat(path, id_to_interval):
    trf = defaultdict(list)
    current = None
    if not Path(path).exists():
        return trf
    with Path(path).open(errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("Sequence:"):
                current = line.split("Sequence:", 1)[1].strip().split()[0]
                continue
            if not current or not line or not line[0].isdigit():
                continue
            p = line.split()
            if len(p) < 14:
                continue
            try:
                rel_start = int(p[0]) - 1
                rel_end = int(p[1])
            except ValueError:
                continue
            if current not in id_to_interval:
                continue
            iv = id_to_interval[current]
            trf[iv["interval_id"]].append({
                "start": iv["start"] + rel_start,
                "end": iv["start"] + rel_end,
                "period": p[2],
                "copies": p[3],
                "consensus_size": p[13],
            })
    return trf


def intersect_features(features, chrom, start, end):
    ivals = []
    class_ivals = defaultdict(list)
    family_bp = defaultdict(int)
    family_class = {}
    n = 0
    for f in features.get(chrom, []):
        if f["end"] <= start:
            continue
        if f["start"] >= end:
            break
        ov = overlap(start, end, f["start"], f["end"])
        if ov <= 0:
            continue
        a = max(start, f["start"])
        b = min(end, f["end"])
        ivals.append((a, b))
        if "major" in f:
            class_ivals[f["major"]].append((a, b))
            family_bp[f["name"]] += ov
            family_class[f["name"]] = f["classification"]
        n += 1
    return ivals, class_ivals, family_bp, family_class, n


def interval_stats(interval, edta, trash, trf_by_interval):
    chrom, start, end = interval["chrom"], interval["start"], interval["end"]
    length = end - start
    edta_ivals, class_ivals, family_bp, family_class, edta_n = intersect_features(edta, chrom, start, end)
    trash_ivals, _, trash_family_bp, _, trash_n = intersect_features(trash, chrom, start, end)
    trf_ivals = [(x["start"], x["end"]) for x in trf_by_interval.get(interval["interval_id"], [])]
    combined_edta_trash_ivals = edta_ivals + trash_ivals
    combined_ivals = edta_ivals + trash_ivals + trf_ivals
    row = {
        "interval_id": interval["interval_id"],
        "label": interval["label"],
        "chrom": chrom,
        "start": start,
        "end": end,
        "length_bp": length,
        "edta_feature_count": edta_n,
        "edta_bp": covered_bp(edta_ivals),
        "edta_pct": covered_bp(edta_ivals) / length * 100,
        "trash_array_count": trash_n,
        "trash_bp": covered_bp(trash_ivals),
        "trash_pct": covered_bp(trash_ivals) / length * 100,
        "trf_array_count": len(trf_ivals),
        "trf_bp": covered_bp(trf_ivals),
        "trf_pct": covered_bp(trf_ivals) / length * 100,
        "combined_edta_trash_bp": covered_bp(combined_edta_trash_ivals),
        "combined_edta_trash_pct": covered_bp(combined_edta_trash_ivals) / length * 100,
        "combined_repeat_bp": covered_bp(combined_ivals),
        "combined_repeat_pct": covered_bp(combined_ivals) / length * 100,
    }
    class_rows = []
    for cls, ivals in sorted(class_ivals.items()):
        class_rows.append({
            "interval_id": interval["interval_id"],
            "repeat_class": cls,
            "bp": covered_bp(ivals),
            "pct": covered_bp(ivals) / length * 100,
        })
    family_rows = []
    for fam, bp in sorted(family_bp.items(), key=lambda kv: -kv[1])[:20]:
        family_rows.append({
            "interval_id": interval["interval_id"],
            "family": fam,
            "classification": family_class.get(fam, "NA"),
            "overlap_bp_unmerged": bp,
            "pct_of_interval_unmerged": bp / length * 100,
        })
    trash_rows = []
    for fam, bp in sorted(trash_family_bp.items(), key=lambda kv: -kv[1])[:20]:
        trash_rows.append({
            "interval_id": interval["interval_id"],
            "trash_array": fam,
            "overlap_bp_unmerged": bp,
            "pct_of_interval_unmerged": bp / length * 100,
        })
    return row, class_rows, family_rows, trash_rows


def random_stats(interval, chrom_lengths, edta, trash, n=1000, seed=20260616):
    rng = random.Random(seed + sum(ord(c) for c in interval["interval_id"]))
    chrom = interval["chrom"]
    length = interval["end"] - interval["start"]
    max_start = chrom_lengths[chrom] - length
    out = []
    for i in range(n):
        s = rng.randint(0, max_start)
        e = s + length
        edta_ivals, _, _, _, _ = intersect_features(edta, chrom, s, e)
        trash_ivals, _, _, _, _ = intersect_features(trash, chrom, s, e)
        out.append({
            "edta_pct": covered_bp(edta_ivals) / length * 100,
            "trash_pct": covered_bp(trash_ivals) / length * 100,
            "combined_edta_trash_pct": covered_bp(edta_ivals + trash_ivals) / length * 100,
        })
    return out


def summarize_random(values, observed, key):
    vals = [v[key] for v in values]
    mean = statistics.mean(vals)
    median = statistics.median(vals)
    sd = statistics.pstdev(vals)
    p95 = np.percentile(vals, 95)
    p_emp = (sum(v >= observed for v in vals) + 1) / (len(vals) + 1)
    z = (observed - mean) / sd if sd > 0 else 0.0
    return {
        f"{key}_random_mean": mean,
        f"{key}_random_median": median,
        f"{key}_random_p95": p95,
        f"{key}_empirical_p_ge_observed": p_emp,
        f"{key}_zscore": z,
    }


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_summary(outdir, rows):
    labels = [r["label"].replace(" ", "\n") for r in rows]
    x = np.arange(len(rows))
    w = 0.18
    fig, ax = plt.subplots(figsize=(5.2, 3.2), constrained_layout=True)
    ax.bar(x - 1.5 * w, [r["edta_pct"] for r in rows], width=w, color=PAL["edta"], label="EDTA TE")
    ax.bar(x - 0.5 * w, [r["trash_pct"] for r in rows], width=w, color=PAL["trash"], label="TRASH array")
    ax.bar(x + 0.5 * w, [r["trf_pct"] for r in rows], width=w, color=PAL["trf"], label="TRF tandem")
    ax.bar(x + 1.5 * w, [r["combined_repeat_pct"] for r in rows], width=w, color=PAL["combined"], label="Combined")
    for i, r in enumerate(rows):
        ax.scatter([i - 1.5 * w], [r["edta_pct_random_mean"]], color=PAL["random"], s=16, zorder=5)
        ax.scatter([i + 1.5 * w], [r["combined_edta_trash_pct_random_mean"]], color=PAL["random"], s=16, zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Interval coverage (%)")
    ax.set_title("Repeat content of Mpo Chr5-to-Chr6 transition intervals", loc="left", fontweight="bold")
    ax.legend(ncol=2, fontsize=6)
    ax.set_ylim(0, min(105, max([r["combined_repeat_pct"] for r in rows] + [80]) + 10))
    for i, r in enumerate(rows):
        ax.text(i + 1.5 * w, r["combined_repeat_pct"] + 2,
                f"{r['combined_repeat_pct']:.1f}%", ha="center", va="bottom", fontsize=6)
    prefix = outdir / "Mpo_Chr5_to_Chr6_transition_repeat_content"
    fig.savefig(str(prefix) + ".png", dpi=450, bbox_inches="tight")
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")


def write_report(path, rows):
    lines = [
        "# Mpo Chr5-to-Chr6 transition repeat intersection",
        "",
        "## Summary",
        "",
        "This analysis intersects the two refined Mpo Chr5-to-Chr6 transition intervals with genome-wide Mpo EDTA TE annotation and TRASH satellite-array annotation, and adds targeted TRF tandem-repeat calls on the extracted transition sequences.",
        "",
    ]
    for r in rows:
        lines.append(
            f"- {r['label']} ({r['chrom']}:{r['start']}-{r['end']}, {r['length_bp'] / 1_000_000:.2f} Mb): "
            f"EDTA TE {r['edta_pct']:.2f}% "
            f"(same-chromosome random mean {r['edta_pct_random_mean']:.2f}%, empirical P>=obs {r['edta_pct_empirical_p_ge_observed']:.3f}); "
            f"TRASH {r['trash_pct']:.2f}%; TRF {r['trf_pct']:.2f}%; "
            f"combined repeat {r['combined_repeat_pct']:.2f}%."
        )
    lines.extend([
        "",
        "## Writing note",
        "",
        "If the observed TE/repeat coverage is higher than the same-chromosome random background, it is reasonable to write that the transition intervals are repeat-rich. This still supports a permissive or associated role for repeats, not a direct proof that TE/repeat was the primary cause of the rearrangement.",
    ])
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--trf-dat", required=True)
    parser.add_argument("--random-n", type=int, default=1000)
    args = parser.parse_args()

    run_root = Path(args.run_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chroms = {x["chrom"] for x in TRANSITIONS}
    chrom_lengths = read_fai(run_root / "data/00_genome/genome_Mpo.fa.fai")
    edta = load_edta(run_root / "data/03_repeat/genome_Mpo.EDTA/genome_Mpo.fa.mod.EDTA.TEanno.gff3", chroms)
    trash = load_trash_bed(run_root / "data/03_repeat/genome_Mpo.TRASH/genome_Mpo.TRASH_arrays.sorted.bed", chroms)
    id_to_interval = {x["interval_id"]: x for x in TRANSITIONS}
    trf = parse_trf_dat(args.trf_dat, id_to_interval)

    summary_rows = []
    class_rows = []
    family_rows = []
    trash_rows = []
    random_rows = []
    for interval in TRANSITIONS:
        row, cls, fam, tr = interval_stats(interval, edta, trash, trf)
        rnd = random_stats(interval, chrom_lengths, edta, trash, n=args.random_n)
        for i, r in enumerate(rnd):
            rr = {"interval_id": interval["interval_id"], "random_index": i}
            rr.update(r)
            random_rows.append(rr)
        row.update(summarize_random(rnd, row["edta_pct"], "edta_pct"))
        row.update(summarize_random(rnd, row["trash_pct"], "trash_pct"))
        row.update(summarize_random(rnd, row["combined_edta_trash_pct"], "combined_edta_trash_pct"))
        summary_rows.append(row)
        class_rows.extend(cls)
        family_rows.extend(fam)
        trash_rows.extend(tr)

    fields = [
        "interval_id", "label", "chrom", "start", "end", "length_bp",
        "edta_feature_count", "edta_bp", "edta_pct",
        "edta_pct_random_mean", "edta_pct_random_median", "edta_pct_random_p95",
        "edta_pct_empirical_p_ge_observed", "edta_pct_zscore",
        "trash_array_count", "trash_bp", "trash_pct",
        "trash_pct_random_mean", "trash_pct_random_median", "trash_pct_random_p95",
        "trash_pct_empirical_p_ge_observed", "trash_pct_zscore",
        "trf_array_count", "trf_bp", "trf_pct",
        "combined_edta_trash_bp", "combined_edta_trash_pct",
        "combined_repeat_bp", "combined_repeat_pct",
        "combined_edta_trash_pct_random_mean", "combined_edta_trash_pct_random_median",
        "combined_edta_trash_pct_random_p95", "combined_edta_trash_pct_empirical_p_ge_observed",
        "combined_edta_trash_pct_zscore",
    ]
    write_tsv(outdir / "transition_repeat_intersect_summary.tsv", summary_rows, fields)
    write_tsv(outdir / "transition_EDTA_class_breakdown.tsv", class_rows,
              ["interval_id", "repeat_class", "bp", "pct"])
    write_tsv(outdir / "transition_EDTA_top_families.tsv", family_rows,
              ["interval_id", "family", "classification", "overlap_bp_unmerged", "pct_of_interval_unmerged"])
    write_tsv(outdir / "transition_TRASH_top_arrays.tsv", trash_rows,
              ["interval_id", "trash_array", "overlap_bp_unmerged", "pct_of_interval_unmerged"])
    write_tsv(outdir / "transition_random_background.tsv", random_rows,
              ["interval_id", "random_index", "edta_pct", "trash_pct", "combined_edta_trash_pct"])
    plot_summary(outdir, summary_rows)
    write_report(outdir / "transition_repeat_intersect_summary_CN.md", summary_rows)


if __name__ == "__main__":
    main()
