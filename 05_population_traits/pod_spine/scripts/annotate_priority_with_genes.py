#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path


def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_attrs(attr):
    out = {}
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
            v = v.strip('"')
        else:
            continue
        out[k] = v
    return out


def load_genes(gff):
    genes = defaultdict(list)
    with open_maybe_gzip(gff) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = parts
            a = parse_attrs(attrs)
            gid = a.get("ID") or a.get("Name") or a.get("gene_id") or f"{chrom}:{start}-{end}"
            name = a.get("Name") or a.get("gene") or a.get("gene_name") or gid
            note = a.get("Note") or a.get("product") or a.get("description") or ""
            genes[chrom].append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "id": gid,
                    "name": name,
                    "note": note,
                }
            )
    index = {}
    for chrom, arr in genes.items():
        arr.sort(key=lambda x: (x["start"], x["end"]))
        starts = [g["start"] for g in arr]
        index[chrom] = (arr, starts)
    return index


def gene_label(g):
    base = f"{g['id']}|{g['chrom']}:{g['start']}-{g['end']}({g['strand']})"
    if g["name"] and g["name"] != g["id"]:
        base += f"|{g['name']}"
    if g["note"]:
        base += f"|{g['note'][:120]}"
    return base


def annotate_interval(index, chrom, start, end, window):
    if not chrom or chrom in {"NA", "."} or not start or not end:
        return {
            "overlap_genes": "NA",
            "overlap_gene_count": "0",
            "nearest_gene": "NA",
            "nearest_gene_distance": "NA",
            "nearest_gene_locus": "NA",
            "gene_context": "no_coordinate",
        }
    try:
        start, end = int(start), int(end)
    except ValueError:
        return {
            "overlap_genes": "NA",
            "overlap_gene_count": "0",
            "nearest_gene": "NA",
            "nearest_gene_distance": "NA",
            "nearest_gene_locus": "NA",
            "gene_context": "bad_coordinate",
        }
    if start > end:
        start, end = end, start
    if chrom not in index:
        return {
            "overlap_genes": "NA",
            "overlap_gene_count": "0",
            "nearest_gene": "NA",
            "nearest_gene_distance": "NA",
            "nearest_gene_locus": "NA",
            "gene_context": "chrom_missing_in_gff",
        }
    arr, starts = index[chrom]
    pos = bisect_right(starts, end)
    scan_start = max(0, pos - 200)
    overlaps = []
    for g in arr[scan_start:pos]:
        if g["end"] >= start and g["start"] <= end:
            overlaps.append(g)
    if overlaps:
        nearest = sorted(overlaps, key=lambda g: (g["start"], g["end"]))[0]
        return {
            "overlap_genes": ",".join(gene_label(g) for g in overlaps[:20]),
            "overlap_gene_count": str(len(overlaps)),
            "nearest_gene": nearest["id"],
            "nearest_gene_distance": "0",
            "nearest_gene_locus": gene_label(nearest),
            "gene_context": "overlap_gene",
        }

    candidates = []
    left_pos = bisect_right(starts, start)
    for i in range(max(0, left_pos - 20), min(len(arr), left_pos + 20)):
        g = arr[i]
        if g["end"] < start:
            dist = start - g["end"]
        elif g["start"] > end:
            dist = g["start"] - end
        else:
            dist = 0
        candidates.append((dist, g))
    if not candidates:
        return {
            "overlap_genes": "NA",
            "overlap_gene_count": "0",
            "nearest_gene": "NA",
            "nearest_gene_distance": "NA",
            "nearest_gene_locus": "NA",
            "gene_context": "no_gene_on_chrom",
        }
    dist, nearest = min(candidates, key=lambda x: (x[0], x[1]["start"]))
    if dist <= 5_000:
        context = "within_5kb"
    elif dist <= 20_000:
        context = "within_20kb"
    elif dist <= window:
        context = f"within_{window//1000}kb"
    else:
        context = "distal_intergenic"
    return {
        "overlap_genes": "NA",
        "overlap_gene_count": "0",
        "nearest_gene": nearest["id"],
        "nearest_gene_distance": str(dist),
        "nearest_gene_locus": gene_label(nearest),
        "gene_context": context,
    }


def as_float(v):
    try:
        return float(v)
    except Exception:
        return 0.0


def context_bonus(context):
    if context == "overlap_gene":
        return 25
    if context == "within_5kb":
        return 18
    if context == "within_20kb":
        return 10
    if context.startswith("within_"):
        return 4
    return 0


def best_context(ref_annos):
    ranked = sorted(
        ref_annos,
        key=lambda x: (
            -context_bonus(x["gene_context"]),
            int(x["nearest_gene_distance"]) if re.match(r"^\d+$", str(x["nearest_gene_distance"])) else 10**12,
        ),
    )
    return ranked[0]


def annotate_table(input_path, output_path, refa_gff, refb_gff, window):
    ref_indexes = {"RefA": load_genes(refa_gff), "RefB": load_genes(refb_gff)}
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        base_fields = reader.fieldnames or []
        rows = []
        for row in reader:
            ref_annos = []
            for ref in ["RefA", "RefB"]:
                anno = annotate_interval(
                    ref_indexes[ref],
                    row.get(f"{ref}_CHROM", ""),
                    row.get(f"{ref}_POS", ""),
                    row.get(f"{ref}_END", ""),
                    window,
                )
                for k, v in anno.items():
                    row[f"{ref}_{k}"] = v
                ref_annos.append({"ref": ref, **anno})
            best = best_context(ref_annos)
            row["best_gene_ref"] = best["ref"]
            row["best_gene_context"] = best["gene_context"]
            row["best_gene_id"] = best["nearest_gene"]
            row["best_gene_distance"] = best["nearest_gene_distance"]
            row["best_gene_locus"] = best["nearest_gene_locus"]
            row["gene_context_bonus"] = context_bonus(best["gene_context"])
            row["gene_validation_score"] = round(as_float(row.get("validation_score")) + row["gene_context_bonus"], 3)
            rows.append(row)

    rows.sort(
        key=lambda r: (
            -as_float(r["gene_validation_score"]),
            -as_float(r.get("validation_score")),
            r.get("best_gene_context", ""),
            r.get("MetaSV_ID", ""),
        )
    )
    for i, row in enumerate(rows, 1):
        row["gene_validation_rank"] = i

    added = [
        "gene_validation_rank",
        "gene_validation_score",
        "gene_context_bonus",
        "best_gene_ref",
        "best_gene_context",
        "best_gene_id",
        "best_gene_distance",
        "best_gene_locus",
        "RefA_gene_context",
        "RefA_overlap_gene_count",
        "RefA_overlap_genes",
        "RefA_nearest_gene",
        "RefA_nearest_gene_distance",
        "RefA_nearest_gene_locus",
        "RefB_gene_context",
        "RefB_overlap_gene_count",
        "RefB_overlap_genes",
        "RefB_nearest_gene",
        "RefB_nearest_gene_distance",
        "RefB_nearest_gene_locus",
    ]
    fields = added + [f for f in base_fields if f not in added]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--refa-gff", required=True)
    parser.add_argument("--refb-gff", required=True)
    parser.add_argument("--window", type=int, default=50_000)
    args = parser.parse_args()
    rows = annotate_table(args.input, args.output, args.refa_gff, args.refb_gff, args.window)
    print(f"wrote\t{args.output}\trows={len(rows)}")


if __name__ == "__main__":
    main()
