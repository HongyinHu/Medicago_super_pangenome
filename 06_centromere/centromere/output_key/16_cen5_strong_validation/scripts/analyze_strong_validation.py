#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


PRIMARY_CHRS = {"Chr%d" % i for i in range(1, 9)}
TRANSITIONS = [
    {
        "name": "Chr3_CEN5_to_CEN6",
        "mpo_chr": "Chr3",
        "left": "Chr5",
        "right": "Chr6",
        "coarse_start": 22075815,
        "coarse_end": 22569156,
    },
    {
        "name": "Chr3_CEN6_to_CEN3",
        "mpo_chr": "Chr3",
        "left": "Chr6",
        "right": "Chr3",
        "coarse_start": 39528301,
        "coarse_end": 40374862,
    },
    {
        "name": "Chr5_CEN5_to_CEN6",
        "mpo_chr": "Chr5",
        "left": "Chr5",
        "right": "Chr6",
        "coarse_start": 25451336,
        "coarse_end": 26794443,
    },
    {
        "name": "Chr5_CEN6_to_CEN3",
        "mpo_chr": "Chr5",
        "left": "Chr6",
        "right": "Chr3",
        "coarse_start": 55703848,
        "coarse_end": 56883318,
    },
]


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def merge_intervals(intervals, gap=0):
    if not intervals:
        return []
    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(a, b) for a, b in merged]


def parse_region_name(name):
    # samtools faidx names look like Chr3:21075816-23569156.
    m = re.match(r"([^:]+):(\d+)-(\d+)$", name)
    if not m:
        return name, 0, None
    chrom = m.group(1)
    start0 = int(m.group(2)) - 1
    end = int(m.group(3))
    return chrom, start0, end


def parse_paf(path):
    rows = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq = p[:12]
            qchrom, qoff, qregion_end = parse_region_name(qname)
            rows.append(
                {
                    "qname": qname,
                    "qchrom": qchrom,
                    "qregion_start": qoff,
                    "qregion_end": qregion_end,
                    "qstart": int(qstart),
                    "qend": int(qend),
                    "qgenomic_start": qoff + int(qstart),
                    "qgenomic_end": qoff + int(qend),
                    "strand": strand,
                    "tname": tname,
                    "tstart": int(tstart),
                    "tend": int(tend),
                    "alen": int(alen),
                    "nmatch": int(nmatch),
                    "identity": int(nmatch) / max(1, int(alen)),
                    "mapq": int(mapq),
                }
            )
    return rows


def summarize_junction(args):
    run = Path(args.run)
    paf = parse_paf(run / "01_junction" / "Mpo_transition_windows_vs_R108_Chr3_5_6.asm20.paf")
    min_aln = int(args.min_aln)
    min_id = float(args.min_identity)
    blocks = []
    boundaries = []
    for tr in TRANSITIONS:
        q0 = tr["coarse_start"] - 2_000_000
        q1 = tr["coarse_end"] + 2_000_000
        local = [
            r
            for r in paf
            if r["qchrom"] == tr["mpo_chr"]
            and r["tname"] in {tr["left"], tr["right"], "Chr3", "Chr5", "Chr6"}
            and r["alen"] >= min_aln
            and r["identity"] >= min_id
            and overlap(r["qgenomic_start"], r["qgenomic_end"], q0, q1) > 0
        ]
        local = sorted(local, key=lambda r: (r["qgenomic_start"], r["qgenomic_end"], -r["alen"]))
        for r in local:
            blocks.append(
                {
                    "transition": tr["name"],
                    "Mpo_chr": tr["mpo_chr"],
                    "Mpo_start": r["qgenomic_start"],
                    "Mpo_end": r["qgenomic_end"],
                    "R108_chr": r["tname"],
                    "R108_start": min(r["tstart"], r["tend"]),
                    "R108_end": max(r["tstart"], r["tend"]),
                    "strand": r["strand"],
                    "aligned_bp": r["alen"],
                    "identity": round(r["identity"], 4),
                    "mapq": r["mapq"],
                }
            )
        left_blocks = [r for r in local if r["tname"] == tr["left"]]
        right_blocks = [r for r in local if r["tname"] == tr["right"]]
        left_before = [r for r in left_blocks if r["qgenomic_start"] <= tr["coarse_end"]]
        right_after = [r for r in right_blocks if r["qgenomic_end"] >= tr["coarse_start"]]
        last_left = max(left_before, key=lambda r: (r["qgenomic_end"], r["alen"])) if left_before else None
        first_right = min(right_after, key=lambda r: (r["qgenomic_start"], -r["alen"])) if right_after else None
        if last_left and first_right:
            b0 = min(last_left["qgenomic_end"], first_right["qgenomic_start"])
            b1 = max(last_left["qgenomic_end"], first_right["qgenomic_start"])
            width = b1 - b0
            status = "sub100kb_candidate" if width <= 100000 else "repeat_gap_interval"
            boundaries.append(
                {
                    "transition": tr["name"],
                    "Mpo_chr": tr["mpo_chr"],
                    "left_ancestry": tr["left"],
                    "right_ancestry": tr["right"],
                    "refined_interval_start": b0,
                    "refined_interval_end": b1,
                    "refined_interval_bp": width,
                    "last_left_block": "%s:%d-%d" % (tr["mpo_chr"], last_left["qgenomic_start"], last_left["qgenomic_end"]),
                    "first_right_block": "%s:%d-%d" % (tr["mpo_chr"], first_right["qgenomic_start"], first_right["qgenomic_end"]),
                    "last_left_identity": round(last_left["identity"], 4),
                    "first_right_identity": round(first_right["identity"], 4),
                    "status": status,
                }
            )
        else:
            boundaries.append(
                {
                    "transition": tr["name"],
                    "Mpo_chr": tr["mpo_chr"],
                    "left_ancestry": tr["left"],
                    "right_ancestry": tr["right"],
                    "refined_interval_start": "NA",
                    "refined_interval_end": "NA",
                    "refined_interval_bp": "NA",
                    "last_left_block": "NA" if not last_left else "%s:%d-%d" % (tr["mpo_chr"], last_left["qgenomic_start"], last_left["qgenomic_end"]),
                    "first_right_block": "NA" if not first_right else "%s:%d-%d" % (tr["mpo_chr"], first_right["qgenomic_start"], first_right["qgenomic_end"]),
                    "last_left_identity": "NA" if not last_left else round(last_left["identity"], 4),
                    "first_right_identity": "NA" if not first_right else round(first_right["identity"], 4),
                    "status": "missing_flank_alignment",
                }
            )
    write_tsv(
        run / "01_junction" / "junction_local_alignment_blocks.tsv",
        blocks,
        ["transition", "Mpo_chr", "Mpo_start", "Mpo_end", "R108_chr", "R108_start", "R108_end", "strand", "aligned_bp", "identity", "mapq"],
    )
    write_tsv(
        run / "results" / "junction_refined_intervals.tsv",
        boundaries,
        [
            "transition",
            "Mpo_chr",
            "left_ancestry",
            "right_ancestry",
            "refined_interval_start",
            "refined_interval_end",
            "refined_interval_bp",
            "last_left_block",
            "first_right_block",
            "last_left_identity",
            "first_right_identity",
            "status",
        ],
    )


def parse_fasta_lengths(path):
    lengths = {}
    name = None
    n = 0
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].strip().split()[0]
                n = 0
            else:
                n += len(line.strip())
        if name is not None:
            lengths[name] = n
    return lengths


def parse_blast(path, target_lengths):
    intervals = defaultdict(list)
    weighted_ident = defaultdict(float)
    aligned = defaultdict(int)
    qfamilies = defaultdict(set)
    if not Path(path).exists():
        return {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            qseqid, sseqid = p[0], p[1]
            pident = float(p[2])
            length = int(p[3])
            sstart, send = int(p[8]), int(p[9])
            a, b = sorted((sstart - 1, send))
            if length < 80:
                continue
            intervals[sseqid].append((a, b))
            weighted_ident[sseqid] += pident * length
            aligned[sseqid] += length
            qfamilies[sseqid].add(qseqid)
    out = {}
    for target, ivals in intervals.items():
        merged = merge_intervals(ivals, gap=20)
        cov = sum(b - a for a, b in merged)
        frag_lengths = [b - a for a, b in ivals]
        out[target] = {
            "trash_hit_fragments": len(ivals),
            "trash_hit_families": len(qfamilies[target]),
            "trash_merged_bp": cov,
            "trash_coverage_pct": cov / max(1, target_lengths.get(target, 1)) * 100,
            "trash_longest_fragment": max(frag_lengths) if frag_lengths else 0,
            "trash_mean_identity": weighted_ident[target] / max(1, aligned[target]),
        }
    return out


def parse_trf_dat(path, target_lengths):
    stats = defaultdict(lambda: {"trf_arrays": 0, "trf_bp_raw": 0, "trf_longest_array": 0, "ivals": []})
    current = None
    if not Path(path).exists():
        return {}
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
                start = int(p[0]) - 1
                end = int(p[1])
            except ValueError:
                continue
            length = max(0, end - start)
            stats[current]["trf_arrays"] += 1
            stats[current]["trf_bp_raw"] += length
            stats[current]["trf_longest_array"] = max(stats[current]["trf_longest_array"], length)
            stats[current]["ivals"].append((start, end))
    out = {}
    for target, st in stats.items():
        merged = merge_intervals(st["ivals"], gap=20)
        cov = sum(b - a for a, b in merged)
        out[target] = {
            "trf_arrays": st["trf_arrays"],
            "trf_merged_bp": cov,
            "trf_coverage_pct": cov / max(1, target_lengths.get(target, 1)) * 100,
            "trf_longest_array": st["trf_longest_array"],
        }
    return out


def parse_target_manifest(path):
    rows = []
    with Path(path).open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            r["start"] = int(r["start"])
            r["end"] = int(r["end"])
            r["length"] = r["end"] - r["start"]
            rows.append(r)
    return rows


def parse_edta_gff(path, chrom, start, end):
    intervals = []
    families = defaultdict(int)
    if not Path(path).exists():
        return {"edta_bp": 0, "edta_pct": 0.0, "edta_top": "NA"}
    with Path(path).open(errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != chrom:
                continue
            s = int(p[3]) - 1
            e = int(p[4])
            ov = overlap(start, end, s, e)
            if ov <= 0:
                continue
            intervals.append((max(start, s) - start, min(end, e) - start))
            attr = p[8]
            fam = "unknown"
            for key in ["Classification=", "Name=", "ID="]:
                if key in attr:
                    fam = attr.split(key, 1)[1].split(";", 1)[0]
                    break
            families[fam] += ov
    merged = merge_intervals(intervals, gap=0)
    cov = sum(b - a for a, b in merged)
    top = "NA"
    if families:
        top = ",".join("%s:%d" % (k, v) for k, v in sorted(families.items(), key=lambda x: -x[1])[:3])
    return {"edta_bp": cov, "edta_pct": cov / max(1, end - start) * 100, "edta_top": top}


def summarize_repeat(args):
    run = Path(args.run)
    repeat_dir = run / "02_repeat_decay"
    target_lengths = parse_fasta_lengths(repeat_dir / "target_regions.fa")
    target_manifest = parse_target_manifest(repeat_dir / "target_regions.tsv")
    blast_stats = parse_blast(repeat_dir / "blast_R108_TRASH_consensus_vs_targets.tsv", target_lengths)
    trf_candidates = list(repeat_dir.glob("target_regions.fa.*.dat"))
    trf_stats = parse_trf_dat(trf_candidates[0], target_lengths) if trf_candidates else {}
    rows = []
    for r in target_manifest:
        name = r["target_id"]
        b = blast_stats.get(name, {})
        t = trf_stats.get(name, {})
        gff = run.parent / "data" / "03_repeat" / ("%s.EDTA" % r["species"]) / ("%s.fa.mod.EDTA.TEanno.gff3" % r["species"])
        e = parse_edta_gff(gff, r["chrom"], r["start"], r["end"])
        rows.append(
            {
                "target_id": name,
                "species": r["species"],
                "class": r["class"],
                "region": "%s:%d-%d" % (r["chrom"], r["start"], r["end"]),
                "length": r["length"],
                "trash_hit_fragments": b.get("trash_hit_fragments", 0),
                "trash_hit_families": b.get("trash_hit_families", 0),
                "trash_merged_bp": int(b.get("trash_merged_bp", 0)),
                "trash_coverage_pct": round(b.get("trash_coverage_pct", 0.0), 4),
                "trash_longest_fragment": int(b.get("trash_longest_fragment", 0)),
                "trash_mean_identity": round(b.get("trash_mean_identity", 0.0), 4),
                "trf_arrays": t.get("trf_arrays", 0),
                "trf_merged_bp": int(t.get("trf_merged_bp", 0)),
                "trf_coverage_pct": round(t.get("trf_coverage_pct", 0.0), 4),
                "trf_longest_array": int(t.get("trf_longest_array", 0)),
                "edta_bp": int(e["edta_bp"]),
                "edta_pct": round(e["edta_pct"], 4),
                "edta_top": e["edta_top"],
            }
        )
    write_tsv(
        run / "results" / "targeted_repeat_decay_summary.tsv",
        rows,
        [
            "target_id",
            "species",
            "class",
            "region",
            "length",
            "trash_hit_fragments",
            "trash_hit_families",
            "trash_merged_bp",
            "trash_coverage_pct",
            "trash_longest_fragment",
            "trash_mean_identity",
            "trf_arrays",
            "trf_merged_bp",
            "trf_coverage_pct",
            "trf_longest_array",
            "edta_bp",
            "edta_pct",
            "edta_top",
        ],
    )


def read_counts(path):
    rows = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end, count, *_ = line.rstrip("\n").split("\t")
            rows.append({"chrom": chrom, "start": int(start), "end": int(end), "count": int(count)})
    return rows


def rolling_mean(values, radius=2):
    out = []
    for i in range(len(values)):
        lo = max(0, i - radius)
        hi = min(len(values), i + radius + 1)
        out.append(sum(values[lo:hi]) / max(1, hi - lo))
    return out


def call_cluster(chrom_rows, min_score=0.25):
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


def call_x8_domains(args):
    run = Path(args.run)
    species_list = args.species.split(",")
    domain_rows = []
    signal_rows = []
    for species in species_list:
        cen = read_counts(run / "03_x8_projection" / "coverage" / ("%s.CENH3.q20.50k.counts.tsv" % species))
        inp = read_counts(run / "03_x8_projection" / "coverage" / ("%s.Input.q20.50k.counts.tsv" % species))
        cen_total = sum(r["count"] for r in cen)
        inp_total = sum(r["count"] for r in inp)
        rows = []
        for c, i in zip(cen, inp):
            if c["chrom"] != i["chrom"] or c["start"] != i["start"] or c["end"] != i["end"]:
                raise RuntimeError("Window mismatch for %s" % species)
            cen_rpm = c["count"] / max(1, cen_total) * 1_000_000
            inp_rpm = i["count"] / max(1, inp_total) * 1_000_000
            log2 = math.log((cen_rpm + 0.25) / (inp_rpm + 0.25), 2)
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
        by_chrom = defaultdict(list)
        for row in rows:
            by_chrom[row["chrom"]].append(row)
        for chrom, chrom_rows in by_chrom.items():
            smoothed = rolling_mean([float(r["log2"]) for r in chrom_rows], radius=2)
            for row, score in zip(chrom_rows, smoothed):
                row["smooth_log2"] = round(score, 5)
        signal_rows.extend(rows)
        write_tsv(
            run / "03_x8_projection" / "signal" / ("%s.q20.50k.log2.tsv" % species),
            rows,
            ["species", "chrom", "start", "end", "CENH3_count", "Input_count", "CENH3_rpm", "Input_rpm", "log2", "smooth_log2"],
        )
        for chrom in sorted([x for x in by_chrom if x in PRIMARY_CHRS]):
            cluster = call_cluster(by_chrom[chrom])
            if cluster:
                cluster["species"] = species
                cluster["domain_id"] = "%s_%s_rawQ20_CENH3" % (species, chrom)
                domain_rows.append(cluster)
    write_tsv(
        run / "results" / "x8_raw_q20_called_CENH3_domains.tsv",
        domain_rows,
        [
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
        ],
    )
    with (run / "results" / "x8_raw_q20_called_CENH3_domains.bed").open("w", encoding="utf-8") as handle:
        for r in domain_rows:
            handle.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (r["chrom"], r["start"], r["end"], r["domain_id"], r["max_smooth_log2"], r["species"]))


def read_bed(path):
    rows = []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            rows.append({"chrom": p[0], "start": int(p[1]), "end": int(p[2]), "name": p[3] if len(p) > 3 else "."})
    return rows


def summarize_x8_projection(args):
    run = Path(args.run)
    active = read_bed(run / "03_x8_projection" / "Mpo_active_CENH3_domains.bed")
    relics = read_bed(run / "03_x8_projection" / "Mpo_CEN5_relic_focus.bed")
    rows = []
    for paf_path in sorted((run / "03_x8_projection" / "paf").glob("*_CEN5_pm3Mb_to_Mpo.asm20.paf")):
        species = paf_path.name.split("_CEN5_pm3Mb_to_Mpo")[0]
        paf_rows = parse_paf(paf_path)
        hits = []
        for r in paf_rows:
            if r["alen"] < 1000 or r["identity"] < 0.70:
                continue
            hits.append({"chrom": r["tname"], "start": min(r["tstart"], r["tend"]), "end": max(r["tstart"], r["tend"]), "alen": r["alen"], "identity": r["identity"]})
        by_chrom = defaultdict(list)
        for h in hits:
            by_chrom[h["chrom"]].append(h)
        merged_rows = []
        for chrom, hs in by_chrom.items():
            intervals = merge_intervals([(h["start"], h["end"]) for h in hs], gap=100000)
            for start, end in intervals:
                supporting = [h for h in hs if overlap(start, end, h["start"], h["end"]) > 0]
                aligned_bp = sum(h["alen"] for h in supporting)
                max_id = max(h["identity"] for h in supporting)
                active_ov = sum(overlap(start, end, a["start"], a["end"]) for a in active if a["chrom"] == chrom)
                relic_ov = sum(overlap(start, end, a["start"], a["end"]) for a in relics if a["chrom"] == chrom)
                merged_rows.append(
                    {
                        "species": species,
                        "Mpo_chr": chrom,
                        "Mpo_start": start,
                        "Mpo_end": end,
                        "length": end - start,
                        "hit_count": len(supporting),
                        "aligned_bp": aligned_bp,
                        "max_identity": round(max_id, 4),
                        "overlap_Mpo_active_CENH3_bp": active_ov,
                        "overlap_Mpo_CEN5_relic_bp": relic_ov,
                    }
                )
        merged_rows = sorted(merged_rows, key=lambda x: (-x["aligned_bp"], x["Mpo_chr"], x["Mpo_start"]))
        rows.extend(merged_rows[:20])
    write_tsv(
        run / "results" / "x8_CEN5_pm3Mb_projection_to_Mpo_summary.tsv",
        rows,
        [
            "species",
            "Mpo_chr",
            "Mpo_start",
            "Mpo_end",
            "length",
            "hit_count",
            "aligned_bp",
            "max_identity",
            "overlap_Mpo_active_CENH3_bp",
            "overlap_Mpo_CEN5_relic_bp",
        ],
    )


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")
    s = sub.add_parser("summarize-junction")
    s.add_argument("--run", required=True)
    s.add_argument("--min-aln", default=200)
    s.add_argument("--min-identity", default=0.70)
    s.set_defaults(func=summarize_junction)
    s = sub.add_parser("summarize-repeat")
    s.add_argument("--run", required=True)
    s.set_defaults(func=summarize_repeat)
    s = sub.add_parser("call-x8-domains")
    s.add_argument("--run", required=True)
    s.add_argument("--species", required=True)
    s.set_defaults(func=call_x8_domains)
    s = sub.add_parser("summarize-x8-projection")
    s.add_argument("--run", required=True)
    s.set_defaults(func=summarize_x8_projection)
    args = p.parse_args()
    if not hasattr(args, "func"):
        p.print_help()
        raise SystemExit(2)
    args.func(args)


if __name__ == "__main__":
    main()
