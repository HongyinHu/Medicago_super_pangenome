#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


CHR_RE = re.compile(r"^Chr(\d+)$")


def natural_chr_key(chrom: str) -> tuple[int, str]:
    m = CHR_RE.match(chrom)
    if m:
        return (int(m.group(1)), chrom)
    return (10**9, chrom)


def read_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip():
                continue
            chrom, size = raw.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(float(size))
    return sizes


def main_chromosomes(sizes: dict[str, int]) -> list[str]:
    chrs = [c for c in sizes if CHR_RE.match(c)]
    return sorted(chrs, key=natural_chr_key)


def read_bedgraph(path: Path, keep_chroms: set[str]) -> list[dict]:
    rows = []
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith(("track", "#")):
                continue
            chrom, start, end, value = raw.rstrip("\n").split()[:4]
            if chrom not in keep_chroms:
                continue
            try:
                signal = float(value)
            except ValueError:
                continue
            if not math.isfinite(signal):
                continue
            rows.append(
                {
                    "chrom": chrom,
                    "start": int(float(start)),
                    "end": int(float(end)),
                    "signal": signal,
                }
            )
    rows.sort(key=lambda r: (natural_chr_key(r["chrom"]), r["start"], r["end"]))
    return rows


def weighted_mean(rows: list[dict]) -> float:
    total = 0.0
    weight = 0
    for r in rows:
        bp = max(0, r["end"] - r["start"])
        total += r["signal"] * bp
        weight += bp
    return total / weight if weight else 0.0


def cluster_high_bins(rows: list[dict], threshold: float, max_gap: int) -> list[dict]:
    high = [r for r in rows if r["signal"] >= threshold]
    clusters: list[list[dict]] = []
    current: list[dict] = []
    for r in high:
        if not current:
            current = [r]
            continue
        last = current[-1]
        if r["chrom"] == last["chrom"] and r["start"] - last["end"] <= max_gap:
            current.append(r)
        else:
            clusters.append(current)
            current = [r]
    if current:
        clusters.append(current)

    out = []
    counters: dict[str, int] = defaultdict(int)
    for xs in clusters:
        chrom = xs[0]["chrom"]
        counters[chrom] += 1
        start = min(x["start"] for x in xs)
        end = max(x["end"] for x in xs)
        covered_bp = sum(max(0, x["end"] - x["start"]) for x in xs)
        score = sum(max(0.0, x["signal"]) * max(0, x["end"] - x["start"]) for x in xs)
        out.append(
            {
                "cluster_id": f"{chrom}_cluster_{counters[chrom]}",
                "chrom": chrom,
                "start": start,
                "end": end,
                "span_bp": end - start,
                "covered_high_bp": covered_bp,
                "bin_count": len(xs),
                "mean_signal": weighted_mean(xs),
                "max_signal": max(x["signal"] for x in xs),
                "score": score,
                "high_fraction_in_span": covered_bp / (end - start) if end > start else 0.0,
            }
        )
    return out


def merge_clusters_to_domains(clusters: list[dict], max_gap: int, min_seed_score: float) -> list[dict]:
    clusters = [c for c in clusters if float(c["score"]) >= min_seed_score]
    clusters.sort(key=lambda r: (natural_chr_key(r["chrom"]), r["start"], r["end"]))
    grouped: list[list[dict]] = []
    current: list[dict] = []
    for c in clusters:
        if not current:
            current = [c]
            continue
        last = current[-1]
        if c["chrom"] == last["chrom"] and c["start"] - last["end"] <= max_gap:
            current.append(c)
        else:
            grouped.append(current)
            current = [c]
    if current:
        grouped.append(current)

    domains = []
    counters: dict[str, int] = defaultdict(int)
    for xs in grouped:
        chrom = xs[0]["chrom"]
        counters[chrom] += 1
        start = min(int(x["start"]) for x in xs)
        end = max(int(x["end"]) for x in xs)
        covered = sum(int(x["covered_high_bp"]) for x in xs)
        score = sum(float(x["score"]) for x in xs)
        max_signal = max(float(x["max_signal"]) for x in xs)
        mean_signal = (
            sum(float(x["mean_signal"]) * int(x["covered_high_bp"]) for x in xs) / covered
            if covered
            else 0.0
        )
        domains.append(
            {
                "domain_id": f"{chrom}_domain_{counters[chrom]}",
                "chrom": chrom,
                "start": start,
                "end": end,
                "span_bp": end - start,
                "covered_high_bp": covered,
                "subcluster_count": len(xs),
                "mean_signal": mean_signal,
                "max_signal": max_signal,
                "score": score,
                "high_fraction_in_span": covered / (end - start) if end > start else 0.0,
                "subclusters": ",".join(x["cluster_id"] for x in xs),
            }
        )
    return domains


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_bed(path: Path) -> list[dict]:
    if not path.exists():
        return []
    rows = []
    with path.open() as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            rows.append(
                {
                    "chrom": parts[0],
                    "start": int(float(parts[1])),
                    "end": int(float(parts[2])),
                    "name": parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}",
                }
            )
    return rows


def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--genome-sizes", required=True)
    parser.add_argument("--old-core-bed", default="")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--cluster-gap", type=int, default=50000)
    parser.add_argument("--domain-gap", type=int, default=1200000)
    parser.add_argument("--min-seed-score", type=float, default=50000)
    parser.add_argument("bedgraphs", nargs="+", help="mode=path.bedGraph")
    args = parser.parse_args()

    species = args.species
    outdir = Path(args.outdir)
    sizes = read_sizes(Path(args.genome_sizes))
    chroms = main_chromosomes(sizes)
    keep = set(chroms)
    old = read_bed(Path(args.old_core_bed)) if args.old_core_bed else []

    all_primary = []
    all_domains = []
    all_clusters = []

    cluster_fields = [
        "mode",
        "cluster_id",
        "chrom",
        "start",
        "end",
        "span_bp",
        "covered_high_bp",
        "bin_count",
        "mean_signal",
        "max_signal",
        "score",
        "high_fraction_in_span",
    ]
    domain_fields = [
        "mode",
        "domain_id",
        "chrom",
        "start",
        "end",
        "span_bp",
        "covered_high_bp",
        "subcluster_count",
        "mean_signal",
        "max_signal",
        "score",
        "high_fraction_in_span",
        "subclusters",
    ]

    for item in args.bedgraphs:
        mode, path_s = item.split("=", 1)
        rows = read_bedgraph(Path(path_s), keep)
        clusters = cluster_high_bins(rows, args.threshold, args.cluster_gap)
        domains = merge_clusters_to_domains(clusters, args.domain_gap, args.min_seed_score)
        for c in clusters:
            c["mode"] = mode
        for d in domains:
            d["mode"] = mode
        all_clusters.extend(clusters)
        all_domains.extend(domains)

        primary_by_chr = {}
        for chrom in chroms:
            candidates = [d for d in domains if d["chrom"] == chrom]
            if candidates:
                candidates.sort(
                    key=lambda r: (
                        float(r["score"]),
                        int(r["covered_high_bp"]),
                        int(r["span_bp"]),
                    ),
                    reverse=True,
                )
                hit = dict(candidates[0])
                hit["mode"] = mode
                hit["primary_rank_basis"] = "max_score_per_chromosome"
                primary_by_chr[chrom] = hit
                all_primary.append(hit)

        write_tsv(outdir / f"{mode}.high_bin_clusters.tsv", clusters, cluster_fields)
        write_tsv(outdir / f"{mode}.broad_domains.tsv", domains, domain_fields)
        write_tsv(
            outdir / f"{mode}.primary_domains_per_chr.tsv",
            list(primary_by_chr.values()),
            domain_fields + ["primary_rank_basis"],
        )

    primary_fields = domain_fields + ["primary_rank_basis"]
    write_tsv(outdir / "all_modes.primary_domains_per_chr.tsv", all_primary, primary_fields)
    write_tsv(outdir / "all_modes.broad_domains.tsv", all_domains, domain_fields)
    write_tsv(outdir / "all_modes.high_bin_clusters.tsv", all_clusters, cluster_fields)

    by_chr_mode = defaultdict(dict)
    for p in all_primary:
        by_chr_mode[p["chrom"]][p["mode"]] = p

    consensus = []
    for chrom in chroms:
        modes = by_chr_mode.get(chrom, {})
        starts = [int(v["start"]) for v in modes.values()]
        ends = [int(v["end"]) for v in modes.values()]
        if not starts:
            continue
        raw = modes.get("raw")
        relaxed = modes.get("relaxed")
        q20 = modes.get("q20")
        old_hits = [o for o in old if o["chrom"] == chrom]
        old_overlap = 0
        old_names = []
        for o in old_hits:
            ov = overlap_bp(min(starts), max(ends), o["start"], o["end"])
            if ov:
                old_overlap += ov
                old_names.append(o["name"])
        consensus.append(
            {
                "chrom": chrom,
                "consensus_start": min(starts),
                "consensus_end": max(ends),
                "consensus_span_bp": max(ends) - min(starts),
                "support_modes": ",".join(sorted(modes)),
                "raw_domain": f"{raw['start']}-{raw['end']}" if raw else "",
                "relaxed_domain": f"{relaxed['start']}-{relaxed['end']}" if relaxed else "",
                "q20_domain": f"{q20['start']}-{q20['end']}" if q20 else "",
                "raw_relaxed_overlap_bp": overlap_bp(int(raw["start"]), int(raw["end"]), int(relaxed["start"]), int(relaxed["end"])) if raw and relaxed else 0,
                "raw_q20_overlap_bp": overlap_bp(int(raw["start"]), int(raw["end"]), int(q20["start"]), int(q20["end"])) if raw and q20 else 0,
                "relaxed_q20_overlap_bp": overlap_bp(int(relaxed["start"]), int(relaxed["end"]), int(q20["start"]), int(q20["end"])) if relaxed and q20 else 0,
                "overlap_old_core_bp": old_overlap,
                "old_core_names": ",".join(old_names),
            }
        )

    consensus_fields = [
        "chrom",
        "consensus_start",
        "consensus_end",
        "consensus_span_bp",
        "support_modes",
        "raw_domain",
        "relaxed_domain",
        "q20_domain",
        "raw_relaxed_overlap_bp",
        "raw_q20_overlap_bp",
        "relaxed_q20_overlap_bp",
        "overlap_old_core_bp",
        "old_core_names",
    ]
    write_tsv(outdir / f"{species}_CENH3_recall_consensus.tsv", consensus, consensus_fields)

    final_bed = outdir / f"{species}.CENH3.functional_centromere.recall.consensus.bed"
    with final_bed.open("w") as fh:
        for r in consensus:
            fh.write(
                "\t".join(
                    [
                        r["chrom"],
                        str(r["consensus_start"]),
                        str(r["consensus_end"]),
                        f"{r['chrom']}_CENH3_recall_consensus",
                        r["support_modes"],
                    ]
                )
                + "\n"
            )

    md = [
        f"# {species} CENH3 Recall Summary",
        "",
        "Input signal: CENH3/Input log2 enrichment recalculated from raw, relaxed and q20 BAMs with 10 kb bins and 50 kb smoothing.",
        "",
        f"Main chromosomes detected from idxstats: {', '.join(chroms)}.",
        "",
        f"Final BED: `{final_bed.name}`.",
        "",
        "Consensus domains:",
        "",
        "| Chromosome | Consensus interval | Span (Mb) | Support modes | Old core overlap (bp) |",
        "|---|---:|---:|---|---:|",
    ]
    for r in consensus:
        md.append(
            f"| {r['chrom']} | {r['consensus_start']}-{r['consensus_end']} | "
            f"{int(r['consensus_span_bp']) / 1e6:.2f} | {r['support_modes']} | {r['overlap_old_core_bp']} |"
        )
    (outdir / f"{species}_CENH3_recall_summary.md").write_text("\n".join(md) + "\n")


if __name__ == "__main__":
    main()
