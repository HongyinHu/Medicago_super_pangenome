#!/usr/bin/env python3
import argparse
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


def choose(columns, *candidates):
    lookup = {column.upper(): column for column in columns}
    for candidate in candidates:
        if candidate.upper() in lookup:
            return lookup[candidate.upper()]
    raise RuntimeError(f"None of {candidates} found in columns {list(columns)}")


def load_assoc(label, path):
    frame = pd.read_csv(path, sep=r"\s+", low_memory=False)
    chrom = choose(frame.columns, "CHR", "CHROM")
    marker = choose(frame.columns, "SNP", "ID", "MARKER")
    pos = choose(frame.columns, "POS", "BP")
    pval = choose(frame.columns, "PVAL", "P", "PVALUE")
    rename = {chrom: "chrom", marker: "id", pos: "pos", pval: "p_value"}
    frame = frame.rename(columns=rename)
    frame["class"] = label
    frame["chrom"] = frame["chrom"].astype(str).str.replace("Chr", "", regex=False)
    frame["chrom"] = pd.to_numeric(frame["chrom"], errors="coerce")
    frame["pos"] = pd.to_numeric(frame["pos"], errors="coerce")
    frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce")
    frame = frame[
        frame["chrom"].between(1, 8)
        & frame["pos"].notna()
        & frame["p_value"].gt(0)
        & frame["p_value"].le(1)
    ].copy()
    frame["chrom"] = frame["chrom"].astype(int)
    frame["pos"] = frame["pos"].astype(int)
    frame["minus_log10_p"] = -np.log10(frame["p_value"].clip(lower=np.nextafter(0, 1)))
    if "SCORE" in frame.columns and "VAR" in frame.columns:
        var = pd.to_numeric(frame["VAR"], errors="coerce")
        score = pd.to_numeric(frame["SCORE"], errors="coerce")
        frame["beta_score_approx"] = score / var
        frame["se_score_approx"] = np.sqrt(1 / var)
    return frame


def lambda_gc(values):
    p = np.asarray(values, dtype=float)
    p = p[(p > 0) & (p < 1)]
    return float(np.median(chi2.isf(p, 1)) / 0.4549364) if len(p) else float("nan")


def cumulative(frame):
    maxima = frame.groupby("chrom")["pos"].max().reindex(range(1, 9)).fillna(0)
    offsets, centers, total = {}, {}, 0
    for chrom in range(1, 9):
        offsets[chrom] = total
        centers[chrom] = total + maxima.loc[chrom] / 2
        total += int(maxima.loc[chrom]) + 5_000_000
    data = frame.copy()
    data["cum_pos"] = [pos + offsets[chrom] for chrom, pos in zip(data.chrom, data.pos)]
    return data, offsets, centers


def draw_manhattan(frames, out, prefix):
    fig, axes = plt.subplots(len(frames), 1, figsize=(12, 3.4 * len(frames)), squeeze=False)
    colors = ("#75BFE8", "#E83E8C")
    for ax, (label, frame) in zip(axes[:, 0], frames.items()):
        data, offsets, centers = cumulative(frame)
        for chrom in range(1, 9):
            part = data[data.chrom == chrom]
            ax.scatter(part.cum_pos, part.minus_log10_p, s=2, alpha=0.68,
                       color=colors[(chrom - 1) % 2], edgecolors="none", rasterized=True)
        n = len(data)
        ax.axhline(-math.log10(0.05 / n), color="#D73027", ls="--", lw=0.9)
        ax.axhline(-math.log10(1 / n), color="#666666", ls=":", lw=0.9)
        for gene, (chrom, start, end) in GENES.items():
            x = offsets[chrom] + (start + end) / 2
            ax.axvline(x, color="#333333", lw=0.65)
            ax.text(x, ax.get_ylim()[1] * 0.97, gene, rotation=90, ha="center", va="top", fontsize=8)
        ax.set_ylabel(r"$-\log_{10}(P)$")
        ax.set_title(f"{label} logistic mixed-model GWAS (n markers={n:,})", loc="left")
        ax.set_xticks([centers[c] for c in range(1, 9)])
        ax.set_xticklabels(range(1, 9))
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1, 0].set_xlabel("Chromosome")
    fig.tight_layout()
    fig.savefig(out / f"{prefix}.manhattan.png", dpi=400)
    fig.savefig(out / f"{prefix}.manhattan.pdf", dpi=300)
    plt.close(fig)


def draw_qq(frames, out, prefix):
    fig, axes = plt.subplots(1, len(frames), figsize=(4.5 * len(frames), 4), squeeze=False)
    for ax, (label, frame) in zip(axes[0], frames.items()):
        obs = np.sort(frame.p_value.to_numpy())
        exp = (np.arange(1, len(obs) + 1) - 0.5) / len(obs)
        if len(obs) > 100_000:
            idx = np.unique(np.linspace(0, len(obs) - 1, 100_000).astype(int))
            obs, exp = obs[idx], exp[idx]
        x = -np.log10(exp)
        y = -np.log10(np.clip(obs, np.nextafter(0, 1), 1))
        limit = max(x.max(), y.max()) * 1.02
        ax.plot([0, limit], [0, limit], color="#777777", ls="--", lw=0.8)
        ax.scatter(x, y, s=4, color="#3182BD", alpha=0.55, edgecolors="none", rasterized=True)
        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"Expected $-\log_{10}(P)$")
        ax.set_ylabel(r"Observed $-\log_{10}(P)$")
        ax.set_title(f"{label}, $\lambda_{{GC}}$={lambda_gc(frame.p_value):.3f}")
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"{prefix}.qq.png", dpi=400)
    fig.savefig(out / f"{prefix}.qq.pdf", dpi=300)
    plt.close(fig)


def candidate_summary(frames, out):
    rows = []
    for label, frame in frames.items():
        n = len(frame)
        for gene, (chrom, start, end) in GENES.items():
            same_chr = frame[frame.chrom == chrom].copy()
            same_chr["distance_bp"] = np.where(
                same_chr.pos < start,
                start - same_chr.pos,
                np.where(same_chr.pos > end, same_chr.pos - end, 0),
            )
            for window in (0, 50_000, 500_000, 1_000_000):
                local = same_chr[same_chr.distance_bp <= window]
                row = {
                    "class": label,
                    "gene_id": gene,
                    "chrom": chrom,
                    "gene_start": start,
                    "gene_end": end,
                    "window_bp": window,
                    "tested_markers": n,
                    "markers_in_window": len(local),
                }
                if len(local):
                    best = local.loc[local.p_value.idxmin()]
                    row.update(best_id=best.id, best_pos=int(best.pos),
                               best_distance_bp=int(best.distance_bp), best_p=float(best.p_value),
                               best_minus_log10_p=float(best.minus_log10_p))
                rows.append(row)
    pd.DataFrame(rows).to_csv(out / "candidate_gene_peak_summary.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assoc", action="append", required=True, help="CLASS=path")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    frames = {}
    for item in args.assoc:
        label, path = item.split("=", 1)
        frames[label] = load_assoc(label, path)
    combined = pd.concat(frames.values(), ignore_index=True).sort_values("p_value")
    combined.to_csv(out / f"{args.prefix}.assoc.tsv.gz", sep="\t", index=False, compression="gzip")
    combined.head(500).to_csv(out / f"{args.prefix}.top500.tsv", sep="\t", index=False)
    with open(out / f"{args.prefix}.summary.tsv", "w", encoding="utf-8") as handle:
        handle.write("class\ttested\tbonferroni_0.05\tsuggestive_1_over_n\tlambda_gc\tbest_id\tbest_chr\tbest_pos\tbest_p\n")
        for label, frame in frames.items():
            best = frame.loc[frame.p_value.idxmin()]
            handle.write(f"{label}\t{len(frame)}\t{0.05/len(frame):.12g}\t{1/len(frame):.12g}\t"
                         f"{lambda_gc(frame.p_value):.6g}\t{best.id}\t{int(best.chrom)}\t"
                         f"{int(best.pos)}\t{best.p_value:.12g}\n")
    candidate_summary(frames, out)
    draw_manhattan(frames, out, args.prefix)
    draw_qq(frames, out, args.prefix)
    print(f"Summarized {len(combined):,} GMMAT tests across {len(frames)} classes")


if __name__ == "__main__":
    main()
