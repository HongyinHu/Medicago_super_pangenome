#!/usr/bin/env python3
import csv
import os
import re

BASE = "path/to/project/N_4.pod_spiny"
OUT = os.path.join(BASE, "13_Chr23997_spiny_intron_window_comparison_20260708")
TE_DATA = "path/to/project/39.TE_type_soloLTR/data"

HOMOLOGY_TSV = os.path.join(OUT, "results", "Chr23997_spiny_window_homology.tsv")

EDTA_CANDIDATES = {
    "genome_Msa": [
        os.path.join(BASE, "00_data/3.two_ref/genome_Msa_EDTA"),
        os.path.join(TE_DATA, "genome_Msa_T2T_EDTA"),
        os.path.join(TE_DATA, "genome_Msa_EDTA"),
    ],
    "genome_R108": [
        os.path.join(BASE, "00_data/3.two_ref/genome_R108_EDTA"),
        os.path.join(TE_DATA, "genome_R108_EDTA"),
    ],
    "genome_410": [
        os.path.join(TE_DATA, "genome_410_EDTA"),
    ],
    "genome_474a": [
        os.path.join(TE_DATA, "genome_474_T2T_EDTA"),
        os.path.join(TE_DATA, "genome_474_EDTA"),
    ],
    "genome_474b": [
        os.path.join(TE_DATA, "genome_474_T2T_EDTA"),
        os.path.join(TE_DATA, "genome_474_EDTA"),
    ],
    "genome_Mpo": [
        os.path.join(TE_DATA, "genome_MPO_EDTA"),
    ],
}


def find_files(root, pattern):
    out = []
    if not root or not os.path.isdir(root):
        return out
    root = os.path.realpath(root)
    for cur, _, files in os.walk(root):
        depth = os.path.relpath(cur, root).count(os.sep)
        if depth > 4:
            continue
        for f in files:
            if re.fullmatch(pattern, f):
                out.append(os.path.join(cur, f))
    return sorted(out)


def parse_coord(coord):
    if not coord or ":" not in coord or "-" not in coord:
        return "", 0, 0
    seqid, rest = coord.split(":", 1)
    start, end = rest.split("-", 1)
    return seqid, int(start), int(end)


def parse_attrs(attrs):
    parsed = {}
    for part in re.split(r";\s*", attrs.strip()):
        if not part:
            continue
        if "=" in part:
            key, val = part.split("=", 1)
        elif " " in part:
            key, val = part.split(" ", 1)
        else:
            key, val = part, ""
        parsed[key.strip()] = val.strip().strip('"')
    return parsed


def short_label(attrs):
    parsed = parse_attrs(attrs)
    vals = []
    for key in ("Classification", "classification", "Class", "class", "Name", "ID", "Target", "Motif", "Note"):
        if parsed.get(key):
            vals.append(parsed[key])
    if vals:
        return ";".join(vals[:5])
    return attrs[:180]


def overlap_gff(path, query):
    rows = []
    qseq, qs, qe = parse_coord(query["target_coord"])
    if not qseq:
        return rows
    with open(path, errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[0] != qseq:
                continue
            try:
                s = int(parts[3])
                e = int(parts[4])
            except ValueError:
                continue
            if s <= qe and e >= qs:
                rows.append(
                    {
                        "genome": query["genome"],
                        "species": query["species"],
                        "role": query["role"],
                        "window": query["window"],
                        "target_coord": query["target_coord"],
                        "target_len": query["target_len"],
                        "annotation_file": path,
                        "seqid": parts[0],
                        "source": parts[1],
                        "feature_type": parts[2],
                        "te_start": s,
                        "te_end": e,
                        "overlap_bp": min(e, qe) - max(s, qs) + 1,
                        "strand": parts[6],
                        "annotation": short_label(parts[8]),
                        "attributes": parts[8],
                    }
                )
    return rows


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main():
    os.makedirs(os.path.join(OUT, "results"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "summary"), exist_ok=True)

    with open(HOMOLOGY_TSV) as handle:
        queries = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("status") in ("ok", "reference")]

    source_rows = []
    overlap_rows = []
    for genome, roots in EDTA_CANDIDATES.items():
        files = []
        for root in roots:
            files.extend((root, "TEanno_gff3", p) for p in find_files(root, r".*TEanno\.gff3"))
            files.extend((root, "RepeatMasker_out_gff", p) for p in find_files(root, r".*\.out\.gff"))
        seen = set()
        uniq = []
        for root, kind, path in files:
            real = os.path.realpath(path)
            if real in seen:
                continue
            seen.add(real)
            uniq.append((root, kind, path))
        for root, kind, path in uniq:
            source_rows.append({"genome": genome, "edta_root": root, "kind": kind, "annotation_file": path})
        for query in [q for q in queries if q["genome"] == genome]:
            for root, kind, path in uniq:
                for row in overlap_gff(path, query):
                    row["edta_root"] = root
                    row["annotation_kind"] = kind
                    overlap_rows.append(row)

    fields = [
        "genome",
        "species",
        "role",
        "window",
        "target_coord",
        "target_len",
        "edta_root",
        "annotation_kind",
        "annotation_file",
        "seqid",
        "source",
        "feature_type",
        "te_start",
        "te_end",
        "overlap_bp",
        "strand",
        "annotation",
        "attributes",
    ]
    overlap_tsv = os.path.join(OUT, "results", "Chr23997_spiny_window_EDTA_RepeatMasker_overlap.tsv")
    write_tsv(overlap_tsv, overlap_rows, fields)

    source_tsv = os.path.join(OUT, "results", "Chr23997_spiny_window_EDTA_sources.tsv")
    write_tsv(source_tsv, source_rows, ["genome", "edta_root", "kind", "annotation_file"])

    summary_md = os.path.join(OUT, "summary", "Chr23997_spiny_window_TE_overlap_summary.md")
    with open(summary_md, "w") as handle:
        handle.write("# Chr23997 spiny intron-window local TE annotation overlap\n\n")
        handle.write(f"Homology table: `{HOMOLOGY_TSV}`\n\n")
        handle.write(f"Overlap table: `{overlap_tsv}`\n\n")
        handle.write(f"EDTA source table: `{source_tsv}`\n\n")
        if not overlap_rows:
            handle.write("No EDTA/RepeatMasker feature overlaps were found in the extracted windows.\n")
        else:
            for row in overlap_rows:
                handle.write(
                    f"- {row['genome']} {row['window']} {row['target_coord']}: "
                    f"{row['annotation_kind']} {row['feature_type']} {row['seqid']}:{row['te_start']}-{row['te_end']}, "
                    f"overlap={row['overlap_bp']} bp, annotation={row['annotation']}\n"
                )

    print("SOURCES", source_tsv)
    print("OVERLAP_TSV", overlap_tsv)
    print("SUMMARY", summary_md)
    print("OVERLAPS", len(overlap_rows))
    for row in overlap_rows:
        print(
            "OV",
            row["genome"],
            row["window"],
            row["target_coord"],
            row["annotation_kind"],
            row["feature_type"],
            f"{row['seqid']}:{row['te_start']}-{row['te_end']}",
            row["overlap_bp"],
            row["annotation"][:120],
        )


if __name__ == "__main__":
    main()
