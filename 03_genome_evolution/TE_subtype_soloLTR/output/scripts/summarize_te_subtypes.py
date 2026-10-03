#!/usr/bin/env python3
import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote


UNKNOWN_VALUES = {"", "?", "na", "nan", "none", "unknown", "unclassified"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize EDTA TE annotation lengths by TEsorter subtype."
    )
    parser.add_argument("--base", required=True, help="39.TE_type_soloLTR base directory")
    parser.add_argument("--db", default="rexdb-plant", help="TEsorter database name used in prefixes")
    parser.add_argument("--output-dir", required=True, help="Directory for summary TSV files")
    return parser.parse_args()


def first_file(directory, patterns):
    for pattern in patterns:
        hits = [hit for hit in sorted(Path(directory).glob(pattern)) if hit.is_file() and hit.stat().st_size > 0]
        if hits:
            return hits[0]
    return None


def fasta_lengths(path):
    lengths = {}
    name = None
    total = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0]
                total = 0
            else:
                total += len(line)
        if name is not None:
            lengths[name] = total
    return lengths


def parse_fasta_classes(path):
    classes = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            token = line[1:].strip().split()[0]
            if "#" in token:
                seq_id, cls = token.split("#", 1)
            else:
                seq_id, cls = token, ""
            classes[seq_id] = cls
    return classes


def parse_attrs(attr_text):
    attrs = {}
    for part in attr_text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            continue
        attrs[key] = unquote(value)
    return attrs


def attr(attrs, key, default=""):
    return attrs.get(key, attrs.get(key.lower(), default))


def normalize_name(name):
    return name.split("#", 1)[0].strip()


def companion_ids(name):
    name = normalize_name(name)
    candidates = [name]
    if name.endswith("_LTR"):
        candidates.append(name[:-4] + "_INT")
        candidates.append(name[:-4])
    elif name.endswith("_INT"):
        candidates.append(name[:-4] + "_LTR")
        candidates.append(name[:-4])
    else:
        candidates.extend([name + "_INT", name + "_LTR"])
    return candidates


def load_tesorter_cls(path):
    rows = {}
    if not path or not Path(path).exists():
        return rows
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            te_id = normalize_name(row.get("#TE") or row.get("TE") or "")
            if not te_id:
                continue
            rows[te_id] = {
                "order": (row.get("Order") or "").strip(),
                "superfamily": (row.get("Superfamily") or "").strip(),
                "clade": (row.get("Clade") or "").strip(),
                "complete": (row.get("Complete") or "").strip(),
                "strand": (row.get("Strand") or "").strip(),
                "domains": (row.get("Domains") or "").strip(),
            }
    return rows


def cls_for_name(name, cls_rows):
    for candidate in companion_ids(name):
        if candidate in cls_rows:
            return cls_rows[candidate]
    return None


def edta_group(edta_class):
    if "/" in edta_class:
        return edta_class.split("/", 1)[0]
    return edta_class or "Unknown"


def subtype_for(name, edta_class, cls_rows):
    row = cls_for_name(name, cls_rows)
    edta_class = edta_class or "Unknown"
    group = edta_group(edta_class)

    if row and row["order"].lower() == "ltr":
        clade = row["clade"].strip()
        superfamily = row["superfamily"].strip()
        if clade.lower() not in UNKNOWN_VALUES:
            return clade, "LTR-RT", row
        if superfamily.lower() in {"copia", "gypsy"}:
            return "Unclassified LTR-RTs", "LTR-RT", row
        return "Unknown", "LTR-RT", row

    if edta_class.startswith("LTR/"):
        if row:
            return "Unclassified LTR-RTs", "LTR-RT", row
        if edta_class.lower() == "ltr/unknown":
            return "Unknown", "LTR-RT", row
        return "Unclassified LTR-RTs", "LTR-RT", row

    if "/" in edta_class:
        return edta_class.split("/", 1)[1], group, row
    return edta_class, group, row


def add_interval(intervals, chrom, start, end, value):
    if end < start:
        return
    intervals[value][chrom].append((start, end))


def merged_bp(chrom_intervals):
    total = 0
    for intervals in chrom_intervals.values():
        if not intervals:
            continue
        intervals = sorted(intervals)
        cur_start, cur_end = intervals[0]
        for start, end in intervals[1:]:
            if start <= cur_end + 1:
                cur_end = max(cur_end, end)
            else:
                total += cur_end - cur_start + 1
                cur_start, cur_end = start, end
        total += cur_end - cur_start + 1
    return total


def percent(numerator, denominator):
    if not denominator:
        return 0.0
    return numerator * 100.0 / denominator


def write_table(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_wide(path, rows, value_field, row_field="species", col_field="subtype"):
    row_names = sorted({row[row_field] for row in rows})
    col_names = sorted({row[col_field] for row in rows})
    lookup = {(row[row_field], row[col_field]): row[value_field] for row in rows}
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([row_field] + col_names)
        for row_name in row_names:
            writer.writerow([row_name] + [lookup.get((row_name, col), "0.000000") for col in col_names])


def summarize_species(species_dir, species, db, out_dir):
    telib = first_file(species_dir, ["*.EDTA.TElib.fa"])
    gff = first_file(species_dir, ["*.EDTA.TEanno.gff3", "*.EDTA.anno/*.EDTA.TEanno.gff3"])
    genome = first_file(species_dir, ["*.genome.fa", "*.ctg.final.fa", "*.fasta.mod", "*.fa.mod", "*.fasta", "*.fa"])
    cls = Path(out_dir) / "tesorter" / species / f"{species}.{db}.cls.tsv"

    if not telib or not gff or not genome:
        missing = []
        if not telib:
            missing.append("TElib")
        if not gff:
            missing.append("TEanno.gff3")
        if not genome:
            missing.append("genome fasta")
        raise FileNotFoundError(f"{species}: missing {', '.join(missing)}")

    genome_bp = sum(fasta_lengths(genome).values())
    fasta_classes = parse_fasta_classes(telib)
    cls_rows = load_tesorter_cls(cls)

    subtype_intervals = defaultdict(lambda: defaultdict(list))
    ltr_intervals = defaultdict(lambda: defaultdict(list))
    all_te_intervals = defaultdict(list)
    bp_sum = Counter()
    counts = Counter()
    scope_by_subtype = {}
    class_counts = Counter()
    matched_tesorter = 0
    ltr_annotation_count = 0

    with open(gff, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, start_s, end_s, attrs_s = fields[0], fields[3], fields[4], fields[8]
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            attrs = parse_attrs(attrs_s)
            name = normalize_name(attr(attrs, "Name") or attr(attrs, "ID") or "")
            edta_class = attr(attrs, "Classification") or fasta_classes.get(name) or "Unknown"
            subtype, scope, row = subtype_for(name, edta_class, cls_rows)
            length = end - start + 1

            add_interval(subtype_intervals, chrom, start, end, subtype)
            all_te_intervals[chrom].append((start, end))
            bp_sum[subtype] += length
            counts[subtype] += 1
            scope_by_subtype[subtype] = scope
            class_counts[edta_class] += 1

            if scope == "LTR-RT":
                add_interval(ltr_intervals, chrom, start, end, subtype)
                ltr_annotation_count += 1
            if row:
                matched_tesorter += 1

    all_ltr_intervals = defaultdict(list)
    for subtype_chroms in ltr_intervals.values():
        for chrom, intervals in subtype_chroms.items():
            all_ltr_intervals[chrom].extend(intervals)

    total_te_bp = merged_bp(all_te_intervals)
    total_ltr_bp = merged_bp(all_ltr_intervals)

    subtype_rows = []
    for subtype in sorted(subtype_intervals):
        bp = merged_bp(subtype_intervals[subtype])
        subtype_rows.append({
            "species": species,
            "subtype": subtype,
            "scope": scope_by_subtype.get(subtype, "Unknown"),
            "bp_merged": str(bp),
            "bp_annotation_sum": str(bp_sum[subtype]),
            "annotation_count": str(counts[subtype]),
            "genome_bp": str(genome_bp),
            "total_te_merged_bp": str(total_te_bp),
            "percent_of_genome": f"{percent(bp, genome_bp):.6f}",
            "percent_of_all_TE": f"{percent(bp, total_te_bp):.6f}",
        })

    ltr_rows = []
    for subtype in sorted(ltr_intervals):
        bp = merged_bp(ltr_intervals[subtype])
        ltr_rows.append({
            "species": species,
            "subtype": subtype,
            "bp_merged": str(bp),
            "bp_annotation_sum": str(bp_sum[subtype]),
            "annotation_count": str(counts[subtype]),
            "genome_bp": str(genome_bp),
            "total_ltr_rt_merged_bp": str(total_ltr_bp),
            "percent_of_genome": f"{percent(bp, genome_bp):.6f}",
            "percent_of_ltr_rt": f"{percent(bp, total_ltr_bp):.6f}",
        })

    total_row = {
        "species": species,
        "genome_fasta": str(genome),
        "telib_fasta": str(telib),
        "teanno_gff3": str(gff),
        "tesorter_cls": str(cls),
        "genome_bp": str(genome_bp),
        "total_te_merged_bp": str(total_te_bp),
        "total_ltr_rt_merged_bp": str(total_ltr_bp),
        "te_percent_of_genome": f"{percent(total_te_bp, genome_bp):.6f}",
        "ltr_rt_percent_of_genome": f"{percent(total_ltr_bp, genome_bp):.6f}",
        "subtype_count": str(len(subtype_rows)),
        "ltr_rt_subtype_count": str(len(ltr_rows)),
        "annotations_matched_to_tesorter": str(matched_tesorter),
        "ltr_rt_annotation_count": str(ltr_annotation_count),
    }
    return subtype_rows, ltr_rows, total_row


def main():
    args = parse_args()
    base = Path(args.base)
    data_dir = base / "data"
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    all_ltr_rows = []
    totals = []
    errors = []

    for species_dir in sorted(data_dir.glob("*_EDTA")):
        species = species_dir.name[:-5]
        try:
            rows, ltr_rows, total_row = summarize_species(species_dir, species, args.db, out_dir)
        except Exception as exc:
            errors.append({"species": species, "error": str(exc)})
            continue
        all_rows.extend(rows)
        all_ltr_rows.extend(ltr_rows)
        totals.append(total_row)

    all_rows.sort(key=lambda r: (r["species"], -float(r["percent_of_all_TE"]), r["subtype"]))
    all_ltr_rows.sort(key=lambda r: (r["species"], -float(r["percent_of_ltr_rt"]), r["subtype"]))
    totals.sort(key=lambda r: r["species"])

    write_table(
        out_dir / "te_subtype_percent.tsv",
        [
            "species", "subtype", "scope", "bp_merged", "bp_annotation_sum",
            "annotation_count", "genome_bp", "total_te_merged_bp",
            "percent_of_genome", "percent_of_all_TE",
        ],
        all_rows,
    )
    write_table(
        out_dir / "ltr_rt_subtype_percent.tsv",
        [
            "species", "subtype", "bp_merged", "bp_annotation_sum",
            "annotation_count", "genome_bp", "total_ltr_rt_merged_bp",
            "percent_of_genome", "percent_of_ltr_rt",
        ],
        all_ltr_rows,
    )
    write_table(
        out_dir / "species_totals.tsv",
        [
            "species", "genome_fasta", "telib_fasta", "teanno_gff3", "tesorter_cls",
            "genome_bp", "total_te_merged_bp", "total_ltr_rt_merged_bp",
            "te_percent_of_genome", "ltr_rt_percent_of_genome",
            "subtype_count", "ltr_rt_subtype_count",
            "annotations_matched_to_tesorter", "ltr_rt_annotation_count",
        ],
        totals,
    )
    write_wide(out_dir / "te_subtype_percent_of_all_TE.wide.tsv", all_rows, "percent_of_all_TE")
    write_wide(out_dir / "te_subtype_percent_of_genome.wide.tsv", all_rows, "percent_of_genome")
    write_wide(out_dir / "ltr_rt_percent.wide.tsv", all_ltr_rows, "percent_of_ltr_rt")

    top_rows = []
    for species in sorted({row["species"] for row in all_rows}):
        rows = [row for row in all_rows if row["species"] == species]
        top_rows.extend(rows[:20])
    write_table(
        out_dir / "top20_te_subtypes.tsv",
        [
            "species", "subtype", "scope", "bp_merged", "bp_annotation_sum",
            "annotation_count", "genome_bp", "total_te_merged_bp",
            "percent_of_genome", "percent_of_all_TE",
        ],
        top_rows,
    )

    if errors:
        write_table(out_dir / "summary_errors.tsv", ["species", "error"], errors)
        raise SystemExit(f"Finished with {len(errors)} species errors; see summary_errors.tsv")


if __name__ == "__main__":
    main()
