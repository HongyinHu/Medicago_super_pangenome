#!/usr/bin/env python3
import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


CATEGORIES = [
    ("Tandem repeat", 1),
    ("LTR/Gypsy", 2),
    ("LTR/Copia", 3),
    ("LTR unknown", 4),
    ("SINEs", 5),
    ("DNA transposons", 6),
    ("Unclassified TE", 7),
]
CAT_TO_CODE = dict(CATEGORIES)
CODE_TO_CAT = {code: cat for cat, code in CATEGORIES}

SKIP_FEATURES = {
    "repeat_region",
    "target_site_duplication",
    "long_terminal_repeat",
    "terminal_inverted_repeat",
    "direct_repeat",
}


def parse_attrs(value):
    attrs = {}
    for item in str(value).split(";"):
        if "=" not in item:
            continue
        key, val = item.split("=", 1)
        attrs[key] = val
    return attrs


def load_regions(path):
    regions = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row["chr"]
            if not re.match(r"^Chr[0-9]+$", chrom):
                continue
            start = int(row["start"])
            end = int(row["end"])
            if end <= start:
                continue
            label = row["label"]
            chr_num = re.sub(r"\D+", "", chrom)
            cen_id = f"{label}.CEN{int(chr_num):02d}" if chr_num else f"{label}.{chrom}"
            regions.append(
                {
                    "genome": row["genome"],
                    "label": label,
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "centromere": cen_id,
                    "size_bp": end - start,
                }
            )
    return regions


def load_config(path):
    df = pd.read_csv(path, sep="\t")
    return {row["genome"]: row.to_dict() for _, row in df.iterrows()}


def classify_te(feature, classification):
    classification = classification or "unknown"
    if classification.startswith("LTR/Gypsy"):
        return "LTR/Gypsy"
    if classification.startswith("LTR/Copia"):
        return "LTR/Copia"
    if classification.startswith("LTR/"):
        return "LTR unknown"
    if classification.startswith("SINE"):
        return "SINEs"
    if classification.startswith("DNA/") or classification.startswith("MITE/") or feature == "helitron":
        return "DNA transposons"
    if classification == "unknown" or feature in {"repeat_fragment", "repeat_region"}:
        return "Unclassified TE"
    return None


def region_index(regions):
    idx = defaultdict(list)
    for i, region in enumerate(regions):
        idx[(region["genome"], region["chr"])].append((i, region["start"], region["end"]))
    return idx


def assign_interval(arrays, region, start, end, code):
    ov_start = max(start, region["start"])
    ov_end = min(end, region["end"])
    if ov_end <= ov_start:
        return
    rel_start = ov_start - region["start"]
    rel_end = ov_end - region["start"]
    view = arrays[rel_start:rel_end]
    view[view == 0] = code
    arrays[rel_start:rel_end] = view


def add_tandem(arrays, regions, idx, genome, tandem_bed):
    if not tandem_bed or pd.isna(tandem_bed) or not Path(tandem_bed).exists():
        return 0
    assigned = 0
    with open(tandem_bed) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            if not re.match(r"^Chr[0-9]+$", chrom):
                continue
            start = int(float(fields[1]))
            end = int(float(fields[2]))
            for region_i, _, _ in idx.get((genome, chrom), []):
                before = int((arrays[region_i] == CAT_TO_CODE["Tandem repeat"]).sum())
                assign_interval(arrays[region_i], regions[region_i], start, end, CAT_TO_CODE["Tandem repeat"])
                after = int((arrays[region_i] == CAT_TO_CODE["Tandem repeat"]).sum())
                assigned += after - before
    return assigned


def add_te(arrays, regions, idx, genome, te_gff):
    assigned_by_cat = defaultdict(int)
    with open(te_gff) as handle:
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
            category = classify_te(feature, attrs.get("classification", "unknown"))
            if category is None:
                continue
            code = CAT_TO_CODE[category]
            start0 = int(start) - 1
            end0 = int(end)
            for region_i, _, _ in idx.get((genome, chrom), []):
                before = int((arrays[region_i] == code).sum())
                assign_interval(arrays[region_i], regions[region_i], start0, end0, code)
                after = int((arrays[region_i] == code).sum())
                assigned_by_cat[category] += after - before
    return assigned_by_cat


def summarize(regions, arrays):
    rows = []
    for region, arr in zip(regions, arrays):
        size = int(region["size_bp"])
        row = {
            "genome": region["genome"],
            "label": region["label"],
            "centromere": region["centromere"],
            "chr": region["chr"],
            "start": region["start"],
            "end": region["end"],
            "size_bp": size,
        }
        known_total = 0
        for category, code in CATEGORIES:
            length = int((arr == code).sum())
            known_total += length
            prefix = category.replace("/", "_").replace(" ", "_")
            row[f"{prefix}_length_bp"] = length
            row[f"{prefix}_ratio_pct"] = length / size * 100 if size else 0
        others = int((arr == 0).sum())
        row["Others_length_bp"] = others
        row["Others_ratio_pct"] = others / size * 100 if size else 0
        row["assigned_length_bp"] = known_total
        row["assigned_ratio_pct"] = known_total / size * 100 if size else 0
        rows.append(row)
    return pd.DataFrame(rows)


def write_long(summary, out_path):
    rows = []
    for _, row in summary.iterrows():
        for category, _ in CATEGORIES:
            prefix = category.replace("/", "_").replace(" ", "_")
            rows.append(
                {
                    "genome": row["genome"],
                    "label": row["label"],
                    "centromere": row["centromere"],
                    "chr": row["chr"],
                    "start": row["start"],
                    "end": row["end"],
                    "size_bp": row["size_bp"],
                    "category": category,
                    "length_bp": row[f"{prefix}_length_bp"],
                    "ratio_pct": row[f"{prefix}_ratio_pct"],
                }
            )
        rows.append(
            {
                "genome": row["genome"],
                "label": row["label"],
                "centromere": row["centromere"],
                "chr": row["chr"],
                "start": row["start"],
                "end": row["end"],
                "size_bp": row["size_bp"],
                "category": "Others",
                "length_bp": row["Others_length_bp"],
                "ratio_pct": row["Others_ratio_pct"],
            }
        )
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--regions", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    regions = load_regions(args.regions)
    config = load_config(args.config)
    idx = region_index(regions)
    arrays = [np.zeros(region["size_bp"], dtype=np.uint8) for region in regions]

    for genome in sorted({r["genome"] for r in regions}):
        genome_regions = [r for r in regions if r["genome"] == genome]
        if genome not in config:
            raise SystemExit(f"Missing config for {genome}")
        cfg = config[genome]
        tandem_bp = add_tandem(arrays, regions, idx, genome, cfg.get("tandem_bed", ""))
        te_bp = add_te(arrays, regions, idx, genome, cfg["teanno_gff"])
        print(f"{genome}: {len(genome_regions)} regions, tandem assigned {tandem_bp:,} bp, TE assigned {sum(te_bp.values()):,} bp")

    summary = summarize(regions, arrays)
    wide_path = outdir / "functional_centromere_sequence_composition.by_chr.wide.tsv"
    long_path = outdir / "functional_centromere_sequence_composition.by_chr.long.tsv"
    avg_path = outdir / "functional_centromere_sequence_composition.average_by_species.tsv"
    summary.to_csv(wide_path, sep="\t", index=False, float_format="%.4f")
    write_long(summary, long_path)

    metric_cols = [c for c in summary.columns if c.endswith("_length_bp") or c.endswith("_ratio_pct") or c == "size_bp"]
    avg = summary.groupby(["genome", "label"], as_index=False)[metric_cols].mean()
    avg.insert(2, "centromere_count", summary.groupby(["genome", "label"]).size().to_numpy())
    avg.to_csv(avg_path, sep="\t", index=False, float_format="%.4f")

    print(f"Wrote: {wide_path}")
    print(f"Wrote: {long_path}")
    print(f"Wrote: {avg_path}")


if __name__ == "__main__":
    main()
