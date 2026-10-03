#!/usr/bin/env python3
import argparse
import csv
import os
from collections import Counter

META_COLS = {
    "SV_ID", "source_methods", "reference", "CHROM", "POS", "END", "SVTYPE", "SVLEN",
    "read_SV_ID", "svgap_ids", "final_confidence"
}


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


def main():
    ap = argparse.ArgumentParser(description="Per-species shared and specific panSV PAV distribution.")
    ap.add_argument("--pav", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--exclude-sample", default="genome_Msa")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)

    shared = Counter()
    specific = Counter()
    total = Counter()
    freq = Counter()
    skipped_no_remaining_presence = 0

    with open(args.pav, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        sample_cols = [c for c in reader.fieldnames if c not in META_COLS and c != args.exclude_sample]
        for row in reader:
            present = [s for s in sample_cols if presence_value(row.get(s))]
            n = len(present)
            if n == 0:
                skipped_no_remaining_presence += 1
                continue
            freq[n] += 1
            if n == 1:
                specific[present[0]] += 1
                total[present[0]] += 1
            else:
                for s in present:
                    shared[s] += 1
                    total[s] += 1

    out_species = args.out_prefix + ".per_species.tsv"
    with open(out_species, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([
            "species", "shared_count", "specific_count", "total_present_count",
            "shared_K", "specific_K", "total_K", "shared_percent", "specific_percent"
        ])
        rows = []
        for s in sample_cols:
            t = total[s]
            rows.append([
                s, shared[s], specific[s], t, shared[s] / 1000.0, specific[s] / 1000.0, t / 1000.0,
                shared[s] / t * 100 if t else 0, specific[s] / t * 100 if t else 0
            ])
        rows.sort(key=lambda x: (-x[2], -x[1], x[0]))
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], f"{r[4]:.6f}", f"{r[5]:.6f}", f"{r[6]:.6f}", f"{r[7]:.6f}", f"{r[8]:.6f}"])

    out_freq = args.out_prefix + ".presence_frequency.tsv"
    with open(out_freq, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["present_species_count", "event_count"])
        for n in range(1, len(sample_cols) + 1):
            w.writerow([n, freq[n]])

    out_run = args.out_prefix + ".run_summary.tsv"
    with open(out_run, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerow(["input_pav", args.pav])
        w.writerow(["excluded_sample", args.exclude_sample])
        w.writerow(["sample_columns", ",".join(sample_cols)])
        w.writerow(["remaining_species_count", len(sample_cols)])
        w.writerow(["events_with_remaining_presence", sum(freq.values())])
        w.writerow(["events_skipped_no_remaining_presence", skipped_no_remaining_presence])
        w.writerow(["definition_specific", "present in exactly one remaining species"])
        w.writerow(["definition_shared", "present in this species and at least one other remaining species"])

    print(f"DONE species={len(sample_cols)} events={sum(freq.values())} skipped={skipped_no_remaining_presence}")
    print(out_species)
    print(out_freq)


if __name__ == "__main__":
    main()
