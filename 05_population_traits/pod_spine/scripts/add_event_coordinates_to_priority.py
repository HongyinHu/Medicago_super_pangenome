#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


EVENT_COLUMNS = [
    "RefA_CHROM",
    "RefA_POS",
    "RefA_END",
    "RefB_CHROM",
    "RefB_POS",
    "RefB_END",
]


def load_events(path):
    events = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            events[row["MetaSV_ID"]] = row
    return events


def locus(row):
    if row.get("RefA_CHROM") not in {"", ".", "NA", None}:
        return f"RefA:{row.get('RefA_CHROM')}:{row.get('RefA_POS')}-{row.get('RefA_END')}"
    if row.get("RefB_CHROM") not in {"", ".", "NA", None}:
        return f"RefB:{row.get('RefB_CHROM')}:{row.get('RefB_POS')}-{row.get('RefB_END')}"
    return "NA"


def add_coords(input_path, event_path, output_path):
    events = load_events(event_path)
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        base_fields = reader.fieldnames or []
        fields = []
        inserted = False
        for field in base_fields:
            fields.append(field)
            if field == "MetaSV_ID":
                fields.extend(EVENT_COLUMNS)
                fields.append("candidate_locus")
                inserted = True
        if not inserted:
            fields = ["candidate_locus"] + EVENT_COLUMNS + fields

        rows = []
        missing = 0
        for row in reader:
            event = events.get(row.get("MetaSV_ID", ""))
            if event is None:
                missing += 1
                for col in EVENT_COLUMNS:
                    row[col] = "NA"
            else:
                for col in EVENT_COLUMNS:
                    row[col] = event.get(col, "NA")
            row["candidate_locus"] = locus(row)
            rows.append(row)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows), missing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", required=True)
    parser.add_argument("--inputs", nargs="+", required=True)
    args = parser.parse_args()

    for item in args.inputs:
        src = Path(item)
        if src.suffix != ".tsv":
            raise SystemExit(f"Expected .tsv input: {src}")
        out = src.with_name(src.stem + ".with_coords.tsv")
        n, missing = add_coords(src, args.events, out)
        print(f"{out}\trows={n}\tmissing_events={missing}")


if __name__ == "__main__":
    main()
