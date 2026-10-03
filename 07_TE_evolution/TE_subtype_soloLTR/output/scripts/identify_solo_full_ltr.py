#!/usr/bin/env python3
import argparse
import bisect
import csv
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote


def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify solo-LTR candidates and summarize full-length LTR-RT counts from EDTA output."
    )
    parser.add_argument("--base", required=True, help="39.TE_type_soloLTR base directory")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--merge-gap", type=int, default=100, help="Merge adjacent same-family LTR hits within this gap")
    parser.add_argument("--int-window", type=int, default=20000, help="Exclude solo candidates with an INT hit within this bp window")
    parser.add_argument("--min-ltr-length", type=int, default=80, help="Minimum merged LTR length kept as solo candidate")
    return parser.parse_args()


def parse_attrs(text):
    attrs = {}
    for part in text.rstrip().split(";"):
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


def locate_first(directory, patterns):
    for pattern in patterns:
        hits = [hit for hit in sorted(Path(directory).glob(pattern)) if hit.is_file() and hit.stat().st_size > 0]
        if hits:
            return hits[0]
    return None


def fasta_total_bp(path):
    total = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total


def is_ltr_class(classification):
    return classification.startswith("LTR/")


def family_base(name):
    if name.endswith("_LTR") or name.endswith("_INT"):
        return name[:-4]
    return name


def ltr_group(classification):
    if classification == "LTR/Copia":
        return "LTR/Copia"
    if classification == "LTR/Gypsy":
        return "LTR/Gypsy"
    if classification.startswith("LTR/"):
        return "LTR/unknown"
    return "LTR/unknown"


def gff_records(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            attrs = parse_attrs(fields[8])
            yield {
                "chrom": fields[0],
                "source": fields[1],
                "type": fields[2],
                "start": start,
                "end": end,
                "score": fields[5],
                "strand": fields[6],
                "phase": fields[7],
                "attrs": attrs,
                "line": line.rstrip("\n"),
            }


def load_full_ltr(intact_gff):
    full = []
    for rec in gff_records(intact_gff):
        feature_type = rec["type"]
        attrs = rec["attrs"]
        classification = attr(attrs, "Classification")
        rec_id = attr(attrs, "ID")
        if not is_ltr_class(classification):
            continue
        if "LTR_retrotransposon" not in feature_type and not rec_id.startswith("LTRRT"):
            continue
        full.append({
            "chrom": rec["chrom"],
            "start": rec["start"],
            "end": rec["end"],
            "strand": rec["strand"],
            "name": attr(attrs, "Name", rec_id),
            "classification": classification,
            "group": ltr_group(classification),
            "id": rec_id,
            "length": rec["end"] - rec["start"] + 1,
        })
    return full


def load_teanno_ltr_hits(teanno_gff):
    ltr_hits = []
    int_hits = []
    for rec in gff_records(teanno_gff):
        attrs = rec["attrs"]
        name = attr(attrs, "Name")
        classification = attr(attrs, "Classification")
        if not is_ltr_class(classification):
            continue
        item = {
            "chrom": rec["chrom"],
            "start": rec["start"],
            "end": rec["end"],
            "strand": rec["strand"],
            "name": name,
            "family": family_base(name),
            "classification": classification,
            "group": ltr_group(classification),
            "length": rec["end"] - rec["start"] + 1,
        }
        if name.endswith("_LTR"):
            ltr_hits.append(item)
        elif name.endswith("_INT"):
            int_hits.append(item)
    return ltr_hits, int_hits


def merge_ltr_hits(hits, merge_gap):
    merged = []
    hits = sorted(hits, key=lambda x: (x["chrom"], x["name"], x["strand"], x["start"], x["end"]))
    current = None
    class_bp = Counter()
    raw_count = 0
    for hit in hits:
        if (
            current
            and hit["chrom"] == current["chrom"]
            and hit["name"] == current["name"]
            and hit["strand"] == current["strand"]
            and hit["start"] <= current["end"] + merge_gap + 1
        ):
            current["end"] = max(current["end"], hit["end"])
            raw_count += 1
            class_bp[hit["classification"]] += hit["length"]
        else:
            if current:
                current["raw_hit_count"] = raw_count
                current["classification"] = class_bp.most_common(1)[0][0]
                current["group"] = ltr_group(current["classification"])
                current["length"] = current["end"] - current["start"] + 1
                merged.append(current)
            current = {
                "chrom": hit["chrom"],
                "start": hit["start"],
                "end": hit["end"],
                "strand": hit["strand"],
                "name": hit["name"],
                "family": hit["family"],
            }
            class_bp = Counter({hit["classification"]: hit["length"]})
            raw_count = 1
    if current:
        current["raw_hit_count"] = raw_count
        current["classification"] = class_bp.most_common(1)[0][0]
        current["group"] = ltr_group(current["classification"])
        current["length"] = current["end"] - current["start"] + 1
        merged.append(current)
    return merged


def build_index(intervals):
    by_chrom = defaultdict(list)
    for item in intervals:
        by_chrom[item["chrom"]].append((item["start"], item["end"], item))
    index = {}
    for chrom, values in by_chrom.items():
        values.sort(key=lambda value: (value[0], value[1]))
        starts = [value[0] for value in values]
        index[chrom] = (starts, values)
    return index


def overlaps(index, chrom, start, end):
    if chrom not in index:
        return False
    starts, values = index[chrom]
    pos = bisect.bisect_right(starts, end)
    for i in range(max(0, pos - 1), -1, -1):
        s, e, _ = values[i]
        if e < start:
            break
        if s <= end and e >= start:
            return True
    return False


def write_bed(path, records, label_prefix):
    with open(path, "w", encoding="utf-8") as handle:
        for i, rec in enumerate(records, start=1):
            name = f"{label_prefix}_{i}|{rec['name']}|{rec['classification']}"
            fields = [
                rec["chrom"],
                str(rec["start"] - 1),
                str(rec["end"]),
                name,
                "0",
                rec.get("strand", "."),
            ]
            handle.write("\t".join(fields) + "\n")


def ratio(numerator, denominator):
    if denominator == 0:
        return "NA"
    return f"{numerator / denominator:.6f}"


def pct(numerator, denominator):
    if denominator == 0:
        return "0.000000"
    return f"{numerator * 100 / denominator:.6f}"


def summarize_species(species_dir, species, out_dir, args):
    teanno = locate_first(species_dir, ["*.EDTA.TEanno.gff3", "*.EDTA.anno/*.EDTA.TEanno.gff3"])
    intact = locate_first(species_dir, ["*.EDTA.final/*.EDTA.intact.gff3", "*.EDTA.intact.gff3"])
    genome = locate_first(species_dir, ["*.genome.fa", "*.ctg.final.fa", "*.fasta.mod", "*.fa.mod", "*.fasta", "*.fa"])
    if not teanno or not intact or not genome:
        raise FileNotFoundError(f"{species}: missing TEanno/intact/genome file")

    genome_bp = fasta_total_bp(genome)
    full = load_full_ltr(intact)
    ltr_hits, int_hits = load_teanno_ltr_hits(teanno)
    merged_ltr = merge_ltr_hits(ltr_hits, args.merge_gap)
    full_index = build_index(full)
    int_index = build_index(int_hits)

    solo = []
    excluded = Counter()
    for rec in merged_ltr:
        if rec["length"] < args.min_ltr_length:
            excluded["short"] += 1
            continue
        if overlaps(full_index, rec["chrom"], rec["start"], rec["end"]):
            excluded["overlap_full_length_ltr_rt"] += 1
            continue
        if overlaps(int_index, rec["chrom"], rec["start"] - args.int_window, rec["end"] + args.int_window):
            excluded["near_ltr_internal"] += 1
            continue
        solo.append(rec)

    species_out = out_dir / "beds"
    species_out.mkdir(parents=True, exist_ok=True)
    write_bed(species_out / f"{species}.solo_ltr.bed", solo, f"{species}_soloLTR")
    write_bed(species_out / f"{species}.full_length_ltr_rt.bed", full, f"{species}_fullLTRRT")

    full_count = len(full)
    full_bp = sum(item["length"] for item in full)
    solo_count = len(solo)
    solo_bp = sum(item["length"] for item in solo)
    total_count = solo_count + full_count
    total_bp = solo_bp + full_bp

    summary = {
        "species": species,
        "genome_bp": str(genome_bp),
        "teanno_gff3": str(teanno),
        "intact_gff3": str(intact),
        "raw_ltr_hits": str(len(ltr_hits)),
        "merged_ltr_hits": str(len(merged_ltr)),
        "full_length_ltr_rt_count": str(full_count),
        "full_length_ltr_rt_bp": str(full_bp),
        "solo_ltr_count": str(solo_count),
        "solo_ltr_bp": str(solo_bp),
        "solo_count_percent_of_solo_plus_full": pct(solo_count, total_count),
        "full_count_percent_of_solo_plus_full": pct(full_count, total_count),
        "solo_bp_percent_of_solo_plus_full_bp": pct(solo_bp, total_bp),
        "full_bp_percent_of_solo_plus_full_bp": pct(full_bp, total_bp),
        "solo_ltr_percent_of_genome_bp": pct(solo_bp, genome_bp),
        "full_length_ltr_rt_percent_of_genome_bp": pct(full_bp, genome_bp),
        "solo_to_full_count_ratio": ratio(solo_count, full_count),
        "solo_to_full_bp_ratio": ratio(solo_bp, full_bp),
        "excluded_short": str(excluded["short"]),
        "excluded_overlap_full_length_ltr_rt": str(excluded["overlap_full_length_ltr_rt"]),
        "excluded_near_ltr_internal": str(excluded["near_ltr_internal"]),
        "merge_gap_bp": str(args.merge_gap),
        "int_window_bp": str(args.int_window),
        "min_ltr_length_bp": str(args.min_ltr_length),
    }

    by_class = []
    for label, records in [("solo_ltr", solo), ("full_length_ltr_rt", full)]:
        count_by_group = Counter(item["group"] for item in records)
        bp_by_group = Counter()
        for item in records:
            bp_by_group[item["group"]] += item["length"]
        for group in sorted(set(count_by_group) | set(bp_by_group)):
            by_class.append({
                "species": species,
                "category": label,
                "ltr_class": group,
                "count": str(count_by_group[group]),
                "bp": str(bp_by_group[group]),
                "percent_of_category_count": pct(count_by_group[group], len(records)),
                "percent_of_category_bp": pct(bp_by_group[group], sum(item["length"] for item in records)),
                "percent_of_genome_bp": pct(bp_by_group[group], genome_bp),
            })

    return summary, by_class


def write_table(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    base = Path(args.base)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir = base / "data"

    summaries = []
    class_rows = []
    errors = []
    for species_dir in sorted(data_dir.glob("*_EDTA")):
        species = species_dir.name[:-5]
        try:
            summary, by_class = summarize_species(species_dir, species, out_dir, args)
        except Exception as exc:
            errors.append({"species": species, "error": str(exc)})
            continue
        summaries.append(summary)
        class_rows.extend(by_class)

    summary_fields = [
        "species", "genome_bp", "teanno_gff3", "intact_gff3",
        "raw_ltr_hits", "merged_ltr_hits",
        "full_length_ltr_rt_count", "full_length_ltr_rt_bp",
        "solo_ltr_count", "solo_ltr_bp",
        "solo_count_percent_of_solo_plus_full", "full_count_percent_of_solo_plus_full",
        "solo_bp_percent_of_solo_plus_full_bp", "full_bp_percent_of_solo_plus_full_bp",
        "solo_ltr_percent_of_genome_bp", "full_length_ltr_rt_percent_of_genome_bp",
        "solo_to_full_count_ratio", "solo_to_full_bp_ratio",
        "excluded_short", "excluded_overlap_full_length_ltr_rt", "excluded_near_ltr_internal",
        "merge_gap_bp", "int_window_bp", "min_ltr_length_bp",
    ]
    write_table(out_dir / "solo_full_ltr_summary.tsv", summary_fields, summaries)
    write_table(
        out_dir / "solo_full_ltr_by_class.tsv",
        [
            "species", "category", "ltr_class", "count", "bp",
            "percent_of_category_count", "percent_of_category_bp", "percent_of_genome_bp",
        ],
        class_rows,
    )
    if errors:
        write_table(out_dir / "solo_full_ltr_errors.tsv", ["species", "error"], errors)
        raise SystemExit(f"Finished with {len(errors)} errors; see solo_full_ltr_errors.tsv")


if __name__ == "__main__":
    main()
