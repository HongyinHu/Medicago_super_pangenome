#!/usr/bin/env python3
import argparse
import bisect
import csv
import os
from collections import defaultdict, Counter

UPDOWN_WINDOW = 5000
BIN_SIZE = 500

META_COLS = {
    "SV_ID", "source_methods", "reference", "CHROM", "POS", "END", "SVTYPE", "SVLEN",
    "read_SV_ID", "svgap_ids", "final_confidence"
}


def parse_attrs(attr):
    out = {}
    for part in attr.strip().split(";"):
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


def load_gff(gff_path):
    genes = defaultdict(list)
    exons = defaultdict(list)
    mrna_to_gene = {}

    with open(gff_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, ftype, start, end, _score, strand, _phase, attrs = fields[:9]
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            a = parse_attrs(attrs)
            if ftype == "gene":
                gid = a.get("ID") or a.get("Name") or f"{chrom}:{start}-{end}"
                genes[chrom].append((start, end, strand if strand in ("+", "-") else "+", gid))
            elif ftype in ("mRNA", "transcript"):
                mid = a.get("ID")
                parent = a.get("Parent")
                if mid and parent:
                    mrna_to_gene[mid] = parent.split(",")[0]
            elif ftype == "exon":
                parent = a.get("Parent", "")
                if parent:
                    gene_parent = mrna_to_gene.get(parent.split(",")[0], parent.split(",")[0])
                else:
                    gene_parent = a.get("ID") or f"{chrom}:{start}-{end}"
                exons[chrom].append((start, end, gene_parent))

    if not genes:
        raise SystemExit(f"No gene features found in GFF: {gff_path}")

    gene_index = {}
    exon_index = {}
    for chrom, recs in genes.items():
        recs.sort(key=lambda x: (x[0], x[1]))
        starts = [r[0] for r in recs]
        ends_sorted_recs = sorted((r[1], r[0], r[2], r[3]) for r in recs)
        ends = [r[0] for r in ends_sorted_recs]
        prefix_max_end = []
        m = -1
        for rec in recs:
            m = max(m, rec[1])
            prefix_max_end.append(m)
        gene_index[chrom] = {
            "recs": recs,
            "starts": starts,
            "prefix_max_end": prefix_max_end,
            "ends_sorted_recs": ends_sorted_recs,
            "ends": ends,
        }

    for chrom, recs in exons.items():
        recs.sort(key=lambda x: (x[0], x[1]))
        starts = [r[0] for r in recs]
        prefix_max_end = []
        m = -1
        for rec in recs:
            m = max(m, rec[1])
            prefix_max_end.append(m)
        exon_index[chrom] = {"recs": recs, "starts": starts, "prefix_max_end": prefix_max_end}

    return gene_index, exon_index


def find_overlaps(index_by_chrom, chrom, start, end):
    idx = index_by_chrom.get(chrom)
    if not idx:
        return []
    recs = idx["recs"]
    starts = idx["starts"]
    prefix = idx["prefix_max_end"]
    j = bisect.bisect_right(starts, end) - 1
    out = []
    while j >= 0 and prefix[j] >= start:
        rec = recs[j]
        if rec[1] >= start and rec[0] <= end:
            out.append(rec)
        j -= 1
    return out


def nearest_gene(gene_index, chrom, start, end):
    idx = gene_index.get(chrom)
    if not idx:
        return None
    left = None
    li = bisect.bisect_left(idx["ends"], start) - 1
    if li >= 0:
        gene_end, gene_start, strand, gid = idx["ends_sorted_recs"][li]
        left = (start - gene_end, "left", gene_start, gene_end, strand, gid)
    right = None
    ri = bisect.bisect_right(idx["starts"], end)
    if ri < len(idx["recs"]):
        gene_start, gene_end, strand, gid = idx["recs"][ri]
        right = (gene_start - end, "right", gene_start, gene_end, strand, gid)
    candidates = [x for x in (left, right) if x is not None]
    if not candidates:
        return None
    return min(candidates, key=lambda x: (x[0], 0 if x[1] == "left" else 1))


def presence_value(v):
    if v is None:
        return False
    s = str(v).strip()
    if s in ("", ".", "NA", "NaN", "nan"):
        return False
    try:
        return float(s) > 0
    except ValueError:
        return s not in ("0", "./.", "0/0")


def sv_interval(row):
    chrom = row["CHROM"]
    pos = int(float(row["POS"]))
    end_raw = row.get("END", "")
    try:
        end = int(float(end_raw)) if end_raw not in ("", ".", "NA") else pos
    except ValueError:
        end = pos
    if end < pos:
        pos, end = end, pos
    return chrom, pos, end


def distance_bin(category4, signed_dist):
    if category4 == "Gene body":
        return "Gene", 0
    if signed_dist is None:
        return ">5k", UPDOWN_WINDOW + BIN_SIZE
    if signed_dist < -UPDOWN_WINDOW:
        return "<-5k", -UPDOWN_WINDOW - BIN_SIZE
    if signed_dist > UPDOWN_WINDOW:
        return ">5k", UPDOWN_WINDOW + BIN_SIZE
    if signed_dist < 0:
        left = (signed_dist // BIN_SIZE) * BIN_SIZE
        right = left + BIN_SIZE
    else:
        left = (signed_dist // BIN_SIZE) * BIN_SIZE
        right = left + BIN_SIZE
    if left == 0 and signed_dist < 0:
        left = -BIN_SIZE
        right = 0
    label = f"{left/1000:g}k~{right/1000:g}k"
    center = int((left + right) / 2)
    return label, center


def classify(row, gene_index, exon_index):
    chrom, start, end = sv_interval(row)
    gene_overlaps = find_overlaps(gene_index, chrom, start, end)
    exon_overlaps = find_overlaps(exon_index, chrom, start, end)

    if gene_overlaps:
        best = max(gene_overlaps, key=lambda g: min(end, g[1]) - max(start, g[0]) + 1)
        return {
            "category4": "Gene body",
            "gene_body_detail": "exon_overlap" if exon_overlaps else "intron_only",
            "nearest_gene_id": best[3],
            "nearest_gene_start": best[0],
            "nearest_gene_end": best[1],
            "nearest_gene_strand": best[2],
            "distance_to_gene_bp": 0,
            "signed_distance_bp": 0,
        }

    ng = nearest_gene(gene_index, chrom, start, end)
    if ng is None:
        return {
            "category4": "Intergenic",
            "gene_body_detail": "not_gene_body",
            "nearest_gene_id": ".",
            "nearest_gene_start": ".",
            "nearest_gene_end": ".",
            "nearest_gene_strand": ".",
            "distance_to_gene_bp": ".",
            "signed_distance_bp": ".",
        }

    dist, side, gstart, gend, strand, gid = ng
    if side == "left":
        signed = dist if strand == "+" else -dist
    else:
        signed = -dist if strand == "+" else dist

    if dist <= UPDOWN_WINDOW:
        category4 = "Upstream" if signed < 0 else "Downstream"
    else:
        category4 = "Intergenic"

    return {
        "category4": category4,
        "gene_body_detail": "not_gene_body",
        "nearest_gene_id": gid,
        "nearest_gene_start": gstart,
        "nearest_gene_end": gend,
        "nearest_gene_strand": strand,
        "distance_to_gene_bp": dist,
        "signed_distance_bp": signed,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pav", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--exclude-sample", default="genome_Msa")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    gene_index, exon_index = load_gff(args.gff)

    event_out = args.out_prefix + ".events.tsv"
    category_out = args.out_prefix + ".category_counts.tsv"
    detail_out = args.out_prefix + ".gene_body_detail_counts.tsv"
    dist_out = args.out_prefix + ".distance_bins.tsv"
    svtype_cat_out = args.out_prefix + ".svtype_by_category.tsv"
    summary_out = args.out_prefix + ".summary.tsv"

    category_counts = Counter()
    detail_counts = Counter()
    dist_counts = Counter()
    svtype_cat = Counter()
    total = 0
    skipped_no_remaining_presence = 0

    with open(args.pav, "r", encoding="utf-8", errors="replace", newline="") as inf, \
         open(event_out, "w", encoding="utf-8", newline="") as outf:
        reader = csv.DictReader(inf, delimiter="\t")
        sample_cols = [c for c in reader.fieldnames if c not in META_COLS and c != args.exclude_sample]
        writer_fields = [
            "SV_ID", "CHROM", "POS", "END", "SVTYPE", "SVLEN", "source_methods", "final_confidence",
            "remaining_sample_presence_count", "category4", "gene_body_detail", "nearest_gene_id",
            "nearest_gene_start", "nearest_gene_end", "nearest_gene_strand", "distance_to_gene_bp",
            "signed_distance_bp", "distance_bin", "distance_bin_center_bp"
        ]
        writer = csv.DictWriter(outf, fieldnames=writer_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for row in reader:
            present = sum(1 for c in sample_cols if presence_value(row.get(c)))
            if present == 0:
                skipped_no_remaining_presence += 1
                continue
            total += 1
            info = classify(row, gene_index, exon_index)
            signed = info["signed_distance_bp"]
            signed_for_bin = None if signed == "." else int(signed)
            bin_label, bin_center = distance_bin(info["category4"], signed_for_bin)

            category_counts[info["category4"]] += 1
            detail_counts[info["gene_body_detail"]] += 1
            dist_counts[(bin_center, bin_label)] += 1
            svtype_cat[(row.get("SVTYPE", "."), info["category4"])] += 1

            outrow = {k: row.get(k, ".") for k in ("SV_ID", "CHROM", "POS", "END", "SVTYPE", "SVLEN", "source_methods", "final_confidence")}
            outrow.update(info)
            outrow["remaining_sample_presence_count"] = present
            outrow["distance_bin"] = bin_label
            outrow["distance_bin_center_bp"] = bin_center
            writer.writerow(outrow)

    with open(category_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["category", "count", "percent"])
        for cat in ["Upstream", "Downstream", "Gene body", "Intergenic"]:
            cnt = category_counts[cat]
            w.writerow([cat, cnt, f"{(cnt / total * 100) if total else 0:.4f}"])

    with open(detail_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_body_detail", "count", "percent"])
        for cat in ["exon_overlap", "intron_only", "not_gene_body"]:
            cnt = detail_counts[cat]
            w.writerow([cat, cnt, f"{(cnt / total * 100) if total else 0:.4f}"])

    ordered_bins = [(-UPDOWN_WINDOW - BIN_SIZE, "<-5k")]
    x = -UPDOWN_WINDOW
    while x < 0:
        ordered_bins.append((x + BIN_SIZE // 2, f"{x/1000:g}k~{(x+BIN_SIZE)/1000:g}k"))
        x += BIN_SIZE
    ordered_bins.append((0, "Gene"))
    x = 0
    while x < UPDOWN_WINDOW:
        ordered_bins.append((x + BIN_SIZE // 2, f"{x/1000:g}k~{(x+BIN_SIZE)/1000:g}k"))
        x += BIN_SIZE
    ordered_bins.append((UPDOWN_WINDOW + BIN_SIZE, ">5k"))

    with open(dist_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["bin_label", "bin_center_bp", "event_count"])
        for center, label in ordered_bins:
            w.writerow([label, center, dist_counts.get((center, label), 0)])

    with open(svtype_cat_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["SVTYPE", "category", "count"])
        for svtype in sorted({k[0] for k in svtype_cat}):
            for cat in ["Upstream", "Downstream", "Gene body", "Intergenic"]:
                w.writerow([svtype, cat, svtype_cat.get((svtype, cat), 0)])

    with open(summary_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerow(["input_pav", args.pav])
        w.writerow(["gff", args.gff])
        w.writerow(["excluded_sample", args.exclude_sample])
        w.writerow(["remaining_sample_columns", ",".join(sample_cols)])
        w.writerow(["events_used_remaining_presence_gt0", total])
        w.writerow(["events_skipped_no_remaining_presence", skipped_no_remaining_presence])
        w.writerow(["updown_window_bp", UPDOWN_WINDOW])
        w.writerow(["distance_bin_size_bp", BIN_SIZE])

    print(f"DONE events_used={total} skipped_no_remaining_presence={skipped_no_remaining_presence}")
    print(event_out)
    print(category_out)
    print(dist_out)


if __name__ == "__main__":
    main()
