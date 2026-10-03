#!/usr/bin/env python3
import argparse
import gzip
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2


GENES = {
    "Chr14120": (3, 53_196_386, 53_203_316),
    "Chr24229": (4, 91_173_307, 91_178_034),
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assoc", action="append", required=True, help="CLASS=path")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", required=True)
    return parser.parse_args()


def load_assoc(label, path):
    usecols = ["chr", "rs", "ps", "n_miss", "af", "beta", "se", "p_wald"]
    frame = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    frame = frame.rename(columns={"chr": "chrom", "rs": "id", "ps": "pos"})
    frame["class"] = label
    frame["chrom"] = frame["chrom"].astype(str).str.replace("Chr", "", regex=False)
    frame["chrom"] = pd.to_numeric(frame["chrom"], errors="coerce")
    frame["pos"] = pd.to_numeric(frame["pos"], errors="coerce")
    frame["p_wald"] = pd.to_numeric(frame["p_wald"], errors="coerce")
    frame = frame[
        frame["chrom"].between(1, 8)
        & frame["pos"].notna()
        & frame["p_wald"].gt(0)
        & frame["p_wald"].le(1)
    ].copy()
    frame["chrom"] = frame["chrom"].astype(int)
    frame["pos"] = frame["pos"].astype(int)
    frame["minus_log10_p"] = -np.log10(frame["p_wald"].clip(lower=np.nextafter(0, 1)))
    return frame


def lambda_gc(pvalues):
    values = np.asarray(pvalues, dtype=float)
    values = values[(values > 0) & (values < 1)]
    if len(values) == 0:
        return float("nan")
    return float(np.median(chi2.isf(values, 1)) / 0.4549364)


def interval_distance(pos, start, end):
    pos = np.asarray(pos)
    return np.where(pos < start, start - pos, np.where(pos > end, pos - end, 0))


def candidate_summary(frames, out):
    rows = []
    for label, frame in frames.items():
        n = len(frame)
        bonf = 0.05 / n
        suggestive = 1.0 / n
        for gene, (chrom, start, end) in GENES.items():
            sub = frame[frame["chrom"] == chrom].copy()
            sub["distance_to_gene_bp"] = interval_distance(sub["pos"], start, end)
            record = {
                "class": label,
                "gene_id": gene,
                "chrom": chrom,
                "gene_start": start,
                "gene_end": end,
                "tested_markers": n,
                "bonferroni_0.05": bonf,
                "suggestive_1_over_n": suggestive,
            }
            for window in (0, 50_000, 500_000, 1_000_000):
                local = sub[sub["distance_to_gene_bp"] <= window]
                key = "gene_body" if window == 0 else f"plusminus_{window // 1000}kb"
                record[f"{key}_markers"] = len(local)
                if len(local):
                    best = local.loc[local["p_wald"].idxmin()]
                    record[f"{key}_best_id"] = best["id"]
                    record[f"{key}_best_pos"] = int(best["pos"])
                    record[f"{key}_best_distance_bp"] = int(best["distance_to_gene_bp"])
                    record[f"{key}_best_p"] = float(best["p_wald"])
                    record[f"{key}_best_minus_log10_p"] = float(best["minus_log10_p"])
                    if best["p_wald"] <= bonf:
                        status = "genome_wide_significant"
                    elif best["p_wald"] <= suggestive:
                        status = "suggestive"
                    else:
                        status = "not_significant"
                    record[f"{key}_status"] = status
                else:
                    for suffix in ("best_id", "best_pos", "best_distance_bp", "best_p", "best_minus_log10_p"):
                        record[f"{key}_{suffix}"] = "NA"
                    record[f"{key}_status"] = "no_marker"
            rows.append(record)
    pd.DataFrame(rows).to_csv(out / "candidate_gene_peak_summary.tsv", sep="\t", index=False)


def add_cumulative(frame):
    maxima = frame.groupby("chrom")["pos"].max().reindex(range(1, 9)).fillna(0).astype(int)
    offsets = {}
    total = 0
    for chrom in range(1, 9):
        offsets[chrom] = total
        total += int(maxima.loc[chrom]) + 5_000_000
    result = frame.copy()
    result["cum_pos"] = result.apply(lambda row: row["pos"] + offsets[int(row["chrom"])], axis=1)
    centers = {
        chrom: offsets[chrom] + int(maxima.loc[chrom]) / 2 for chrom in range(1, 9)
    }
    return result, offsets, centers


def manhattan(frames, out, prefix):
    colors = ["#75BFE8", "#E83E8C"]
    fig, axes = plt.subplots(len(frames), 1, figsize=(12, 3.4 * len(frames)), squeeze=False)
    for ax, (label, frame) in zip(axes[:, 0], frames.items()):
        data, offsets, centers = add_cumulative(frame)
        for chrom in range(1, 9):
            part = data[data["chrom"] == chrom]
            ax.scatter(
                part["cum_pos"], part["minus_log10_p"], s=2.2, alpha=0.65,
                color=colors[(chrom - 1) % 2], edgecolors="none", rasterized=True,
            )
        n = len(data)
        ax.axhline(-math.log10(0.05 / n), color="#D73027", linestyle="--", linewidth=0.9, label="Bonferroni 0.05")
        ax.axhline(-math.log10(1 / n), color="#6B6B6B", linestyle=":", linewidth=0.9, label="Suggestive 1/N")
        ymax = max(4.0, float(data["minus_log10_p"].quantile(0.9999)) * 1.15)
        top = float(data["minus_log10_p"].max())
        ymax = max(ymax, min(top * 1.08, 40))
        for gene, (chrom, start, end) in GENES.items():
            x = offsets[chrom] + (start + end) / 2
            ax.axvline(x, color="#333333", linewidth=0.7, alpha=0.7)
            ax.text(x, ymax * 0.98, gene, ha="center", va="top", fontsize=8, rotation=90)
        ax.set_ylim(0, ymax)
        ax.set_ylabel(r"$-\log_{10}(P)$")
        ax.set_title(f"{label} GWAS (n={n:,})", loc="left", fontweight="bold")
        ax.set_xticks([centers[c] for c in range(1, 9)])
        ax.set_xticklabels([str(c) for c in range(1, 9)])
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1, 0].set_xlabel("Chromosome")
    fig.tight_layout()
    fig.savefig(out / f"{prefix}.manhattan.png", dpi=400)
    fig.savefig(out / f"{prefix}.manhattan.pdf", dpi=300)
    plt.close(fig)


def qqplot(frames, out, prefix):
    fig, axes = plt.subplots(1, len(frames), figsize=(4.4 * len(frames), 4), squeeze=False)
    for ax, (label, frame) in zip(axes[0], frames.items()):
        obs = np.sort(frame["p_wald"].to_numpy())
        n = len(obs)
        exp = (np.arange(1, n + 1) - 0.5) / n
        if n > 100_000:
            idx = np.unique(np.linspace(0, n - 1, 100_000).astype(int))
            obs, exp = obs[idx], exp[idx]
        x = -np.log10(exp)
        y = -np.log10(np.clip(obs, np.nextafter(0, 1), 1))
        limit = max(float(x.max()), float(y.max())) * 1.02
        ax.plot([0, limit], [0, limit], color="#777777", linestyle="--", linewidth=0.8)
        ax.scatter(x, y, s=4, color="#3182BD", alpha=0.55, edgecolors="none", rasterized=True)
        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
        ax.set_xlabel(r"Expected $-\log_{10}(P)$")
        ax.set_ylabel(r"Observed $-\log_{10}(P)$")
        ax.set_title(f"{label}, $\lambda_{{GC}}$={lambda_gc(frame['p_wald']):.3f}")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"{prefix}.qq.png", dpi=400)
    fig.savefig(out / f"{prefix}.qq.pdf", dpi=300)
    plt.close(fig)


def regional(frames, out, prefix):
    fig, axes = plt.subplots(len(frames), len(GENES), figsize=(12, 3.5 * len(frames)), squeeze=False)
    for row_index, (label, frame) in enumerate(frames.items()):
        n = len(frame)
        for col_index, (gene, (chrom, start, end)) in enumerate(GENES.items()):
            ax = axes[row_index, col_index]
            lo, hi = start - 1_000_000, end + 1_000_000
            part = frame[(frame["chrom"] == chrom) & frame["pos"].between(lo, hi)]
            ax.scatter(part["pos"] / 1e6, part["minus_log10_p"], s=8, color="#3182BD", alpha=0.7, edgecolors="none", rasterized=True)
            ax.axvspan(start / 1e6, end / 1e6, color="#D73027", alpha=0.16)
            ax.axhline(-math.log10(0.05 / n), color="#D73027", linestyle="--", linewidth=0.8)
            ax.axhline(-math.log10(1 / n), color="#777777", linestyle=":", linewidth=0.8)
            ax.set_title(f"{label}: {gene} (Chr{chrom})")
            ax.set_xlabel("Position (Mb)")
            ax.set_ylabel(r"$-\log_{10}(P)$")
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"{prefix}.candidate_regions.png", dpi=400)
    fig.savefig(out / f"{prefix}.candidate_regions.pdf", dpi=300)
    plt.close(fig)


def main():
    args = parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    frames = {}
    for item in args.assoc:
        label, path = item.split("=", 1)
        frames[label] = load_assoc(label, path)
    combined = pd.concat(frames.values(), ignore_index=True)
    combined = combined.sort_values("p_wald")
    combined.to_csv(out / f"{args.prefix}.assoc.tsv.gz", sep="\t", index=False, compression="gzip")
    combined.head(500).to_csv(out / f"{args.prefix}.top500.tsv", sep="\t", index=False)

    with open(out / f"{args.prefix}.summary.tsv", "w", encoding="utf-8") as handle:
        handle.write("class\ttested\tbonferroni_0.05\tsuggestive_1_over_n\tlambda_gc\tbest_id\tbest_chr\tbest_pos\tbest_p\n")
        for label, frame in frames.items():
            best = frame.loc[frame["p_wald"].idxmin()]
            handle.write(
                f"{label}\t{len(frame)}\t{0.05/len(frame):.12g}\t{1/len(frame):.12g}\t"
                f"{lambda_gc(frame['p_wald']):.6g}\t{best['id']}\t{int(best['chrom'])}\t"
                f"{int(best['pos'])}\t{best['p_wald']:.12g}\n"
            )

    candidate_summary(frames, out)
    manhattan(frames, out, args.prefix)
    qqplot(frames, out, args.prefix)
    regional(frames, out, args.prefix)
    print(f"Summarized {len(combined):,} tests across {len(frames)} variant classes")


if __name__ == "__main__":
    main()
