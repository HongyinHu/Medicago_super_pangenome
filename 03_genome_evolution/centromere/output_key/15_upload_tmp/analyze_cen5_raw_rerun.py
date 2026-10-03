#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


PRIMARY_CHRS = {f"Chr{i}" for i in range(1, 9)}
FUSION_CHRS = {"Chr3", "Chr5", "Chr6"}
MPO_TARGET_CHRS = {"Chr3", "Chr5"}


def read_counts(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end, count, *_ = line.rstrip("\n").split("\t")
            rows.append({"chrom": chrom, "start": int(start), "end": int(end), "count": int(count)})
    return rows


def rolling_mean(values: list[float], radius: int = 2) -> list[float]:
    out = []
    for i in range(len(values)):
        lo = max(0, i - radius)
        hi = min(len(values), i + radius + 1)
        out.append(sum(values[lo:hi]) / max(1, hi - lo))
    return out


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def call_cluster(chrom_rows: list[dict[str, object]], min_score: float = 0.25) -> dict[str, object] | None:
    if not chrom_rows:
        return None
    scores = [float(r["smooth_log2"]) for r in chrom_rows]
    max_i = max(range(len(scores)), key=scores.__getitem__)
    max_score = scores[max_i]
    if max_score < min_score:
        return None
    threshold = max(min_score, max_score * 0.35)
    left = max_i
    right = max_i
    while left > 0 and scores[left - 1] >= threshold:
        left -= 1
    while right < len(scores) - 1 and scores[right + 1] >= threshold:
        right += 1
    cluster_rows = chrom_rows[left : right + 1]
    return {
        "chrom": cluster_rows[0]["chrom"],
        "start": cluster_rows[0]["start"],
        "end": cluster_rows[-1]["end"],
        "max_window_start": chrom_rows[max_i]["start"],
        "max_window_end": chrom_rows[max_i]["end"],
        "max_smooth_log2": round(max_score, 4),
        "mean_log2": round(sum(float(r["log2"]) for r in cluster_rows) / len(cluster_rows), 4),
        "window_count": len(cluster_rows),
        "threshold": round(threshold, 4),
    }


def call_signal(args: argparse.Namespace) -> None:
    run = Path(args.run)
    outdir = run / "results"
    signal_dir = outdir / "signal"
    domain_rows = []
    all_signal_rows = []
    for species in ["genome_R108", "genome_Mpo"]:
        cen_path = run / "coverage" / f"{species}.CENH3.q20.50k.counts.tsv"
        inp_path = run / "coverage" / f"{species}.Input.q20.50k.counts.tsv"
        cen = read_counts(cen_path)
        inp = read_counts(inp_path)
        if len(cen) != len(inp):
            raise RuntimeError(f"Window count mismatch for {species}")
        cen_total = sum(r["count"] for r in cen)
        inp_total = sum(r["count"] for r in inp)
        rows = []
        for c, i in zip(cen, inp):
            if c["chrom"] != i["chrom"] or c["start"] != i["start"] or c["end"] != i["end"]:
                raise RuntimeError(f"Window mismatch for {species}")
            cen_rpm = (c["count"] / max(1, cen_total)) * 1_000_000
            inp_rpm = (i["count"] / max(1, inp_total)) * 1_000_000
            log2 = math.log2((cen_rpm + 0.25) / (inp_rpm + 0.25))
            rows.append(
                {
                    "species": species,
                    "chrom": c["chrom"],
                    "start": c["start"],
                    "end": c["end"],
                    "CENH3_count": c["count"],
                    "Input_count": i["count"],
                    "CENH3_rpm": round(cen_rpm, 5),
                    "Input_rpm": round(inp_rpm, 5),
                    "log2": round(log2, 5),
                    "smooth_log2": 0,
                }
            )
        by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in rows:
            by_chrom[str(row["chrom"])].append(row)
        for chrom, chrom_rows in by_chrom.items():
            smoothed = rolling_mean([float(r["log2"]) for r in chrom_rows], radius=2)
            for row, score in zip(chrom_rows, smoothed):
                row["smooth_log2"] = round(score, 5)
        write_tsv(
            signal_dir / f"{species}.q20.50k.log2.tsv",
            rows,
            [
                "species",
                "chrom",
                "start",
                "end",
                "CENH3_count",
                "Input_count",
                "CENH3_rpm",
                "Input_rpm",
                "log2",
                "smooth_log2",
            ],
        )
        all_signal_rows.extend(rows)
        for chrom in sorted([x for x in by_chrom if x in PRIMARY_CHRS]):
            cluster = call_cluster(by_chrom[chrom])
            if cluster:
                cluster["species"] = species
                cluster["domain_id"] = f"{species}_{chrom}_rawQ20_CENH3"
                domain_rows.append(cluster)

    fields = [
        "species",
        "domain_id",
        "chrom",
        "start",
        "end",
        "max_window_start",
        "max_window_end",
        "max_smooth_log2",
        "mean_log2",
        "window_count",
        "threshold",
    ]
    write_tsv(outdir / "raw_q20_called_CENH3_domains.tsv", domain_rows, fields)
    with (outdir / "raw_q20_called_CENH3_domains.bed").open("w", encoding="utf-8") as handle:
        for row in domain_rows:
            handle.write(
                f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['domain_id']}\t{row['max_smooth_log2']}\t{row['species']}\n"
            )

    r108_chr5 = [r for r in domain_rows if r["species"] == "genome_R108" and r["chrom"] == "Chr5"]
    if not r108_chr5:
        raise RuntimeError("Could not call R108 Chr5 CENH3 domain from raw q20 signal")
    cen5 = r108_chr5[0]
    sizes = read_sizes(run / "chrom_sizes" / "genome_R108.sizes")
    chrom_len = sizes["Chr5"]
    core_start = int(cen5["start"])
    core_end = int(cen5["end"])
    regions = {
        "core": (core_start, core_end),
        "pm1Mb": (max(0, core_start - 1_000_000), min(chrom_len, core_end + 1_000_000)),
        "pm3Mb": (max(0, core_start - 3_000_000), min(chrom_len, core_end + 3_000_000)),
    }
    for label, (start, end) in regions.items():
        with (outdir / f"R108_CEN5_raw_q20_{label}.bed").open("w", encoding="utf-8") as handle:
            handle.write(f"Chr5\t{start}\t{end}\tR108_CEN5_raw_q20_{label}\n")


def read_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, size, *_ = line.rstrip("\n").split("\t")
            sizes[chrom] = int(size)
    return sizes


def parse_paf(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq = p[:12]
            rows.append(
                {
                    "qname": qname,
                    "qlen": int(qlen),
                    "qstart": int(qstart),
                    "qend": int(qend),
                    "strand": strand,
                    "tname": tname,
                    "tlen": int(tlen),
                    "tstart": int(tstart),
                    "tend": int(tend),
                    "nmatch": int(nmatch),
                    "alen": int(alen),
                    "mapq": int(mapq),
                    "identity": int(nmatch) / max(1, int(alen)),
                }
            )
    return rows


def interval_overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def read_signal_table(path: Path) -> list[dict[str, object]]:
    with path.open() as handle:
        return [
            {
                "chrom": r["chrom"],
                "start": int(r["start"]),
                "end": int(r["end"]),
                "log2": float(r["log2"]),
                "smooth_log2": float(r["smooth_log2"]),
            }
            for r in csv.DictReader(handle, delimiter="\t")
        ]


def mean_signal(signal_rows: list[dict[str, object]], chrom: str, start: int, end: int) -> tuple[float, float, int]:
    values = []
    for row in signal_rows:
        if row["chrom"] != chrom:
            continue
        ov = interval_overlap(start, end, int(row["start"]), int(row["end"]))
        if ov > 0:
            values.append(float(row["smooth_log2"]))
    if not values:
        return (float("nan"), float("nan"), 0)
    return (sum(values) / len(values), max(values), len(values))


def merge_intervals(rows: list[dict[str, object]], gap: int = 100_000) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["target_chr"])].append(row)
    merged = []
    for chrom, items in grouped.items():
        items = sorted(items, key=lambda r: (int(r["target_start"]), int(r["target_end"])))
        current = None
        for row in items:
            if current is None:
                current = {
                    "target_chr": chrom,
                    "start": int(row["target_start"]),
                    "end": int(row["target_end"]),
                    "hit_count": 1,
                    "aligned_bp": int(row["aligned_bp"]),
                    "max_identity": float(row["identity"]),
                }
                continue
            if int(row["target_start"]) <= int(current["end"]) + gap:
                current["end"] = max(int(current["end"]), int(row["target_end"]))
                current["hit_count"] = int(current["hit_count"]) + 1
                current["aligned_bp"] = int(current["aligned_bp"]) + int(row["aligned_bp"])
                current["max_identity"] = max(float(current["max_identity"]), float(row["identity"]))
            else:
                merged.append(current)
                current = {
                    "target_chr": chrom,
                    "start": int(row["target_start"]),
                    "end": int(row["target_end"]),
                    "hit_count": 1,
                    "aligned_bp": int(row["aligned_bp"]),
                    "max_identity": float(row["identity"]),
                }
        if current is not None:
            merged.append(current)
    return merged


def summarize_paf(args: argparse.Namespace) -> None:
    run = Path(args.run)
    outdir = run / "results"
    paf_dir = run / "paf"
    sizes = read_sizes(run / "chrom_sizes" / "genome_Mpo.sizes")
    whole = parse_paf(paf_dir / "Mpo_to_R108.raw_asm5.paf")
    whole = [
        r
        for r in whole
        if r["qname"] in MPO_TARGET_CHRS
        and r["tname"] in FUSION_CHRS
        and int(r["alen"]) >= 5000
        and float(r["identity"]) >= 0.75
    ]

    window_rows = []
    for chrom in ["Chr3", "Chr5"]:
        length = sizes[chrom]
        for start in range(0, length, 1_000_000):
            end = min(length, start + 1_000_000)
            weights = {c: 0 for c in ["Chr5", "Chr6", "Chr3"]}
            for r in whole:
                if r["qname"] != chrom:
                    continue
                ov = interval_overlap(start, end, int(r["qstart"]), int(r["qend"]))
                if ov > 0:
                    weights[str(r["tname"])] += ov
            total = sum(weights.values())
            if total:
                dominant = max(weights, key=weights.get)
                frac = weights[dominant] / total
            else:
                dominant = "NA"
                frac = 0
            window_rows.append(
                {
                    "target_chr": chrom,
                    "window_start": start,
                    "window_end": end,
                    "R108_Chr5_bp": weights["Chr5"],
                    "R108_Chr6_bp": weights["Chr6"],
                    "R108_Chr3_bp": weights["Chr3"],
                    "total_bp": total,
                    "dominant_ancestry": dominant,
                    "dominant_fraction": round(frac, 4),
                }
            )
    write_tsv(
        outdir / "raw_minimap_Mpo_Chr3_Chr5_ancestry_1Mb.tsv",
        window_rows,
        [
            "target_chr",
            "window_start",
            "window_end",
            "R108_Chr5_bp",
            "R108_Chr6_bp",
            "R108_Chr3_bp",
            "total_bp",
            "dominant_ancestry",
            "dominant_fraction",
        ],
    )

    transitions = []
    for chrom in ["Chr3", "Chr5"]:
        rows = [r for r in window_rows if r["target_chr"] == chrom and r["dominant_ancestry"] != "NA" and r["dominant_fraction"] >= 0.45]
        prev = None
        for row in rows:
            if prev and prev["dominant_ancestry"] != row["dominant_ancestry"]:
                transitions.append(
                    {
                        "target_chr": chrom,
                        "left_ancestry": prev["dominant_ancestry"],
                        "right_ancestry": row["dominant_ancestry"],
                        "left_window": f"{prev['window_start']}-{prev['window_end']}",
                        "right_window": f"{row['window_start']}-{row['window_end']}",
                        "candidate_start": prev["window_end"],
                        "candidate_end": row["window_start"],
                        "midpoint": int((int(prev["window_end"]) + int(row["window_start"])) / 2),
                    }
                )
            prev = row
    write_tsv(
        outdir / "raw_minimap_Mpo_Chr3_Chr5_breakpoint_windows.tsv",
        transitions,
        [
            "target_chr",
            "left_ancestry",
            "right_ancestry",
            "left_window",
            "right_window",
            "candidate_start",
            "candidate_end",
            "midpoint",
        ],
    )

    # CEN5 homology/relic candidates.
    signal = read_signal_table(outdir / "signal" / "genome_Mpo.q20.50k.log2.tsv")
    cen5_hits_paf = parse_paf(paf_dir / "R108_CEN5_pm3Mb_to_Mpo.raw_asm20.paf")
    hit_rows = []
    for r in cen5_hits_paf:
        if int(r["alen"]) < 1000 or float(r["identity"]) < 0.70:
            continue
        mean_log2, max_log2, nwin = mean_signal(signal, str(r["tname"]), int(r["tstart"]), int(r["tend"]))
        hit_rows.append(
            {
                "query": r["qname"],
                "query_start": r["qstart"],
                "query_end": r["qend"],
                "target_chr": r["tname"],
                "target_start": min(int(r["tstart"]), int(r["tend"])),
                "target_end": max(int(r["tstart"]), int(r["tend"])),
                "strand": r["strand"],
                "aligned_bp": r["alen"],
                "identity": round(float(r["identity"]), 4),
                "mapq": r["mapq"],
                "Mpo_q20_CENH3_mean_smooth_log2": round(mean_log2, 4) if not math.isnan(mean_log2) else "NA",
                "Mpo_q20_CENH3_max_smooth_log2": round(max_log2, 4) if not math.isnan(max_log2) else "NA",
                "signal_windows": nwin,
            }
        )
    hit_rows = sorted(hit_rows, key=lambda x: (str(x["target_chr"]), int(x["target_start"]), -int(x["aligned_bp"])))
    write_tsv(
        outdir / "R108_CEN5_pm3Mb_raw_minimap_hits_to_Mpo.tsv",
        hit_rows,
        [
            "query",
            "query_start",
            "query_end",
            "target_chr",
            "target_start",
            "target_end",
            "strand",
            "aligned_bp",
            "identity",
            "mapq",
            "Mpo_q20_CENH3_mean_smooth_log2",
            "Mpo_q20_CENH3_max_smooth_log2",
            "signal_windows",
        ],
    )

    merged = merge_intervals(hit_rows, gap=100_000)
    for row in merged:
        mean_log2, max_log2, nwin = mean_signal(signal, str(row["target_chr"]), int(row["start"]), int(row["end"]))
        row["mean_CENH3_log2"] = round(mean_log2, 4) if not math.isnan(mean_log2) else "NA"
        row["max_CENH3_log2"] = round(max_log2, 4) if not math.isnan(max_log2) else "NA"
        row["signal_windows"] = nwin
        row["length"] = int(row["end"]) - int(row["start"])
    merged = sorted(merged, key=lambda r: (str(r["target_chr"]), int(r["start"])))
    write_tsv(
        outdir / "Mpo_CEN5_relic_candidate_merged_intervals.tsv",
        merged,
        [
            "target_chr",
            "start",
            "end",
            "length",
            "hit_count",
            "aligned_bp",
            "max_identity",
            "mean_CENH3_log2",
            "max_CENH3_log2",
            "signal_windows",
        ],
    )
    with (outdir / "Mpo_CEN5_relic_candidate_merged_intervals.bed").open("w", encoding="utf-8") as handle:
        for i, row in enumerate(merged, start=1):
            handle.write(
                f"{row['target_chr']}\t{row['start']}\t{row['end']}\tMpo_CEN5relic_raw_{i};hits={row['hit_count']};bp={row['aligned_bp']};maxCENH3={row['max_CENH3_log2']}\n"
            )


def read_fasta(path: Path) -> dict[str, str]:
    seqs = {}
    name = None
    parts = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name:
            seqs[name] = "".join(parts).upper()
    return seqs


def count_kmers(seq: str, k: int) -> Counter[str]:
    c = Counter()
    for i in range(0, max(0, len(seq) - k + 1)):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        c[kmer] += 1
    return c


def kmer_decay(args: argparse.Namespace) -> None:
    run = Path(args.run)
    outdir = run / "results"
    r108_core = read_fasta(run / "fasta" / "R108_CEN5_raw_q20_core.fa")
    candidate = read_fasta(run / "fasta" / "Mpo_CEN5_relic_candidates.fa")
    active = read_fasta(run / "fasta" / "Mpo_raw_q20_active_CENH3_domains.fa")
    if not r108_core:
        raise RuntimeError("Missing R108 core FASTA")
    ref_seq = "".join(r108_core.values()).upper()
    k = int(args.k)
    ref_counts = count_kmers(ref_seq, k)
    top = {kmer for kmer, count in ref_counts.most_common(int(args.top)) if count >= int(args.min_count)}
    rows = []
    for group, seqs in [("R108_CEN5_core", r108_core), ("Mpo_CEN5_relic_candidates", candidate), ("Mpo_active_CENH3_domains", active)]:
        for name, seq in seqs.items():
            counts = count_kmers(seq, k)
            total_hits = sum(counts.get(kmer, 0) for kmer in top)
            density = total_hits / max(1, len(seq) - k + 1) * 1000
            rows.append(
                {
                    "group": group,
                    "sequence": name,
                    "length": len(seq),
                    "top_CEN5_kmer_set_size": len(top),
                    "top_CEN5_kmer_hits": total_hits,
                    "top_CEN5_kmer_hits_per_kb": round(density, 5),
                }
            )
    write_tsv(
        outdir / "CEN5_specific_kmer_decay_raw.tsv",
        rows,
        [
            "group",
            "sequence",
            "length",
            "top_CEN5_kmer_set_size",
            "top_CEN5_kmer_hits",
            "top_CEN5_kmer_hits_per_kb",
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("call-signal")
    p1.add_argument("--run", required=True)
    p1.set_defaults(func=call_signal)
    p2 = sub.add_parser("summarize-paf")
    p2.add_argument("--run", required=True)
    p2.set_defaults(func=summarize_paf)
    p3 = sub.add_parser("kmer-decay")
    p3.add_argument("--run", required=True)
    p3.add_argument("-k", default=31)
    p3.add_argument("--top", default=2000)
    p3.add_argument("--min-count", default=5)
    p3.set_defaults(func=kmer_decay)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
