#!/usr/bin/env python3
import argparse
import csv
import multiprocessing as mp
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


SKIP_FEATURES = {
    "repeat_region",
    "target_site_duplication",
    "long_terminal_repeat",
    "terminal_inverted_repeat",
    "direct_repeat",
}

BROAD_ORDER = [
    "Tandem repeat",
    "LTR/Gypsy",
    "LTR/Copia",
    "LTR unknown",
    "SINE / LINE",
    "DNA transposons",
    "Unclassified TE",
    "Others",
]

BROAD_PRIORITY = {
    "Tandem repeat": 0,
    "LTR/Gypsy": 10,
    "LTR/Copia": 11,
    "LTR unknown": 12,
    "SINE / LINE": 13,
    "DNA transposons": 14,
    "Unclassified TE": 15,
    "Others": 20,
}


def parse_attrs(text):
    attrs = {}
    for item in str(text).split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def read_regions(path, genome, label, row_prefix):
    regions = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            if not re.match(r"^Chr[0-9]+$", chrom):
                continue
            start = int(fields[1])
            end = int(fields[2])
            if end <= start:
                continue
            chr_num = int(re.sub(r"\D+", "", chrom))
            regions.append(
                {
                    "genome": genome,
                    "label": label,
                    "row_prefix": row_prefix,
                    "centromere": f"{row_prefix}.CEN{chr_num:02d}",
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "size_bp": end - start,
                }
            )
    return regions


def read_tesorter_cls(path):
    cls = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fields = line.rstrip("\n").split(None, 6)
            if len(fields) < 4:
                continue
            while len(fields) < 7:
                fields.append("")
            raw, order, superfamily, clade, complete, strand, domains = fields[:7]
            key = raw.split("#", 1)[0]
            rec = {
                "tesorter_id": key,
                "ts_order": order or "unknown",
                "ts_superfamily": superfamily or "unknown",
                "ts_clade": clade or "unclassified",
                "ts_complete": complete or "",
                "ts_domains": domains or "",
            }
            keys = {key}
            if key.endswith("_INT") or key.endswith("_LTR"):
                keys.add(re.sub(r"_(INT|LTR)$", "", key))
            if key.endswith("_INT"):
                keys.add(re.sub(r"_INT$", "_LTR", key))
            if key.endswith("_LTR"):
                keys.add(re.sub(r"_LTR$", "_INT", key))
            for item in keys:
                cls.setdefault(item, rec)
    return cls


def classify_te(feature, classification, name, cls):
    classification = classification or "unknown"
    ts = cls.get(name)
    if ts is None and (name.endswith("_LTR") or name.endswith("_INT")):
        ts = cls.get(re.sub(r"_(LTR|INT)$", "", name))

    if classification.startswith("LTR/"):
        if ts and ts["ts_order"] == "LTR" and ts["ts_superfamily"] in {"Gypsy", "Copia"}:
            broad = f"LTR/{ts['ts_superfamily']}"
            clade = ts["ts_clade"] if ts["ts_clade"] not in {"", "-", "unknown"} else "unclassified"
            fine = f"{ts['ts_superfamily']}/{clade}"
            return broad, fine, ts
        if classification.startswith("LTR/Gypsy"):
            return "LTR/Gypsy", "Gypsy/unclassified", ts
        if classification.startswith("LTR/Copia"):
            return "LTR/Copia", "Copia/unclassified", ts
        return "LTR unknown", "LTR_unknown", ts

    if classification.startswith("SINE") or classification.startswith("LINE"):
        return "SINE / LINE", classification, ts
    if classification.startswith("DNA/") or classification.startswith("MITE/") or feature == "helitron":
        return "DNA transposons", classification, ts
    if classification == "unknown" or feature in {"repeat_fragment", "repeat_region"}:
        return "Unclassified TE", "Unclassified_TE", ts
    return "Others", classification, ts


def add_interval(intervals_by_chr, chrom, start, end, broad, fine, source, name, feature="", classification="", ts=None):
    if end <= start or not re.match(r"^Chr[0-9]+$", chrom):
        return
    intervals_by_chr[chrom].append(
        {
            "chr": chrom,
            "start": start,
            "end": end,
            "broad_category": broad,
            "fine_category": fine,
            "priority": BROAD_PRIORITY.get(broad, 99),
            "source": source,
            "name": name,
            "feature": feature,
            "classification": classification,
            "ts_superfamily": "" if ts is None else ts.get("ts_superfamily", ""),
            "ts_clade": "" if ts is None else ts.get("ts_clade", ""),
        }
    )


def read_trash(trash_bed, intervals_by_chr):
    with open(trash_bed) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(float(fields[1]))
            end = int(float(fields[2]))
            name = fields[3] if len(fields) > 3 else "TRASH"
            monomer = fields[4] if len(fields) > 4 else "NA"
            fine = f"TRASH_{monomer}bp" if monomer != "NA" else "TRASH"
            add_interval(intervals_by_chr, chrom, start, end, "Tandem repeat", fine, "TRASH", name)


def read_teanno_gff(gff_path, cls, intervals_by_chr):
    with open(gff_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr_text = fields
            if feature in SKIP_FEATURES or not re.match(r"^Chr[0-9]+$", chrom):
                continue
            attrs = parse_attrs(attr_text)
            name = attrs.get("Name", attrs.get("ID", "NA"))
            classification = attrs.get("classification", "unknown")
            broad, fine, ts = classify_te(feature, classification, name, cls)
            add_interval(
                intervals_by_chr,
                chrom,
                int(start) - 1,
                int(end),
                broad,
                fine,
                "EDTA_TEanno",
                name,
                feature=feature,
                classification=classification,
                ts=ts,
            )


def summarize_region(task):
    region, intervals = task
    start = region["start"]
    end = region["end"]
    clipped = []
    breakpoints = {start, end}
    for iv in intervals:
        s = max(start, iv["start"])
        e = min(end, iv["end"])
        if e <= s:
            continue
        new_iv = dict(iv)
        new_iv["clip_start"] = s
        new_iv["clip_end"] = e
        clipped.append(new_iv)
        breakpoints.add(s)
        breakpoints.add(e)

    points = sorted(breakpoints)
    broad_bp = defaultdict(int)
    fine_bp = defaultdict(int)
    interval_rows = []

    for left, right in zip(points[:-1], points[1:]):
        if right <= left:
            continue
        candidates = [iv for iv in clipped if iv["clip_start"] < right and iv["clip_end"] > left]
        if not candidates:
            broad = "Others"
            fine = "Unannotated_or_other"
            source = "unassigned"
            name = ""
            feature = ""
            classification = ""
            ts_superfamily = ""
            ts_clade = ""
        else:
            best = min(candidates, key=lambda x: (x["priority"], x["broad_category"], x["fine_category"], x["name"]))
            broad = best["broad_category"]
            fine = best["fine_category"]
            source = best["source"]
            name = best["name"]
            feature = best.get("feature", "")
            classification = best.get("classification", "")
            ts_superfamily = best.get("ts_superfamily", "")
            ts_clade = best.get("ts_clade", "")
        length = right - left
        broad_bp[broad] += length
        fine_bp[(broad, fine)] += length
        interval_rows.append(
            {
                **region,
                "segment_start": left,
                "segment_end": right,
                "segment_length_bp": length,
                "broad_category": broad,
                "fine_category": fine,
                "source": source,
                "name": name,
                "feature": feature,
                "classification": classification,
                "ts_superfamily": ts_superfamily,
                "ts_clade": ts_clade,
            }
        )

    broad_rows = []
    for broad in BROAD_ORDER:
        length = broad_bp.get(broad, 0)
        broad_rows.append(
            {
                **region,
                "category": broad,
                "length_bp": length,
                "ratio_pct": length / region["size_bp"] * 100,
            }
        )
    fine_rows = []
    for (broad, fine), length in sorted(fine_bp.items(), key=lambda x: (BROAD_PRIORITY.get(x[0][0], 99), x[0][0], x[0][1])):
        fine_rows.append(
            {
                **region,
                "broad_category": broad,
                "fine_category": fine,
                "length_bp": length,
                "ratio_pct": length / region["size_bp"] * 100,
            }
        )
    return broad_rows, fine_rows, interval_rows


def write_wide(broad_df, out_path):
    rows = []
    key_cols = ["genome", "label", "row_prefix", "centromere", "chr", "start", "end", "size_bp"]
    for keys, sub in broad_df.groupby(key_cols, sort=False):
        row = dict(zip(key_cols, keys))
        lookup = sub.set_index("category")
        for broad in BROAD_ORDER:
            safe = broad.replace("/", "_").replace(" ", "_")
            row[f"{safe}_length_bp"] = int(lookup.loc[broad, "length_bp"]) if broad in lookup.index else 0
            row[f"{safe}_ratio_pct"] = float(lookup.loc[broad, "ratio_pct"]) if broad in lookup.index else 0.0
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False, float_format="%.4f")


def write_table_s11(broad_df, outdir, genome):
    rows = []
    for _, sub in broad_df.groupby(["centromere", "chr"], sort=False):
        sub = sub.set_index("category")
        row = {"Centromere": sub["centromere"].iloc[0]}
        for category in BROAD_ORDER:
            row[f"{category} length (bp)"] = int(round(float(sub.loc[category, "length_bp"]))) if category in sub.index else 0
            row[f"{category} ratio (%)"] = round(float(sub.loc[category, "ratio_pct"]), 2) if category in sub.index else 0.0
        rows.append(row)
    table = pd.DataFrame(rows)
    avg = {"Centromere": "Average"}
    for category in BROAD_ORDER:
        avg[f"{category} length (bp)"] = int(round(table[f"{category} length (bp)"].mean()))
        avg[f"{category} ratio (%)"] = round(float(table[f"{category} ratio (%)"].mean()), 2)
    table = pd.concat([table, pd.DataFrame([avg])], ignore_index=True)
    flat = outdir / f"{genome}.functional_centromere_sequence_composition.TableS11_style.flat.tsv"
    pretty = outdir / f"{genome}.functional_centromere_sequence_composition.TableS11_style.tsv"
    table.to_csv(flat, sep="\t", index=False)
    first = ["Centromere"]
    second = [""]
    for category in BROAD_ORDER:
        first.extend([category, ""])
        second.extend(["length (bp)", "ratio (%)"])
    with open(pretty, "w") as handle:
        handle.write("\t".join(first) + "\n")
        handle.write("\t".join(second) + "\n")
        for _, row in table.iterrows():
            values = [row["Centromere"]]
            for category in BROAD_ORDER:
                values.append(str(int(row[f"{category} length (bp)"])))
                values.append(f"{float(row[f'{category} ratio (%)']):.2f}")
            handle.write("\t".join(values) + "\n")
    return table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--row_prefix", required=True)
    parser.add_argument("--centromere_bed", required=True)
    parser.add_argument("--teanno_gff", required=True)
    parser.add_argument("--tesorter_cls", required=True)
    parser.add_argument("--trash_bed", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    regions = read_regions(args.centromere_bed, args.genome, args.label, args.row_prefix)
    cls = read_tesorter_cls(args.tesorter_cls)
    intervals_by_chr = defaultdict(list)
    read_trash(args.trash_bed, intervals_by_chr)
    read_teanno_gff(args.teanno_gff, cls, intervals_by_chr)

    tasks = []
    for region in regions:
        chr_intervals = [
            iv
            for iv in intervals_by_chr.get(region["chr"], [])
            if iv["start"] < region["end"] and iv["end"] > region["start"]
        ]
        tasks.append((region, chr_intervals))

    workers = max(1, min(args.threads, len(tasks)))
    if workers == 1:
        results = [summarize_region(task) for task in tasks]
    else:
        with mp.Pool(processes=workers) as pool:
            results = pool.map(summarize_region, tasks)

    broad_rows, fine_rows, interval_rows = [], [], []
    for broad, fine, intervals in results:
        broad_rows.extend(broad)
        fine_rows.extend(fine)
        interval_rows.extend(intervals)

    broad_df = pd.DataFrame(broad_rows)
    fine_df = pd.DataFrame(fine_rows)
    interval_df = pd.DataFrame(interval_rows)

    broad_long = outdir / f"{args.genome}.functional_centromere_sequence_composition.broad.long.tsv"
    broad_wide = outdir / f"{args.genome}.functional_centromere_sequence_composition.broad.wide.tsv"
    fine_path = outdir / f"{args.genome}.functional_centromere_sequence_composition.TEsorter_fine.tsv"
    segment_path = outdir / f"{args.genome}.functional_centromere_sequence_composition.assigned_segments.tsv"

    broad_df.to_csv(broad_long, sep="\t", index=False, float_format="%.4f")
    write_wide(broad_df, broad_wide)
    fine_df.to_csv(fine_path, sep="\t", index=False, float_format="%.4f")
    interval_df.to_csv(segment_path, sep="\t", index=False)

    avg_broad = (
        broad_df.groupby("category", as_index=False)
        .agg(length_bp=("length_bp", "mean"), ratio_pct=("ratio_pct", "mean"))
        .assign(genome=args.genome, label=args.label)
    )
    avg_broad[["genome", "label", "category", "length_bp", "ratio_pct"]].to_csv(
        outdir / f"{args.genome}.functional_centromere_sequence_composition.broad.average.tsv",
        sep="\t",
        index=False,
        float_format="%.4f",
    )
    table = write_table_s11(broad_df, outdir, args.genome)

    print(f"{args.genome}: regions={len(regions)}, workers={workers}, tesorter_entries={len(cls)}")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
