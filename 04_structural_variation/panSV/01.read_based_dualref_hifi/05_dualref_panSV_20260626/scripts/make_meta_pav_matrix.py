#!/usr/bin/env python3
import argparse
import os
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a meta-panSV PAV matrix from cross-reference events and RefA/RefB PAV matrices."
    )
    parser.add_argument("--events", required=True)
    parser.add_argument("--refa-matrix", required=True)
    parser.add_argument("--refb-matrix", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--stats", required=True)
    return parser.parse_args()


def load_pav(path):
    data = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if len(header) <= 6 or header[0] != "SV_ID":
            raise SystemExit(f"Unexpected PAV matrix header in {path}")
        samples = header[6:]
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            values = row[6:]
            if len(values) < len(samples):
                values += ["NA"] * (len(samples) - len(values))
            data[row[0]] = values[: len(samples)]
    return samples, data


def align_values(values, source_samples, output_samples):
    index = {sample: i for i, sample in enumerate(source_samples)}
    return [
        values[index[sample]] if sample in index and index[sample] < len(values) else "NA"
        for sample in output_samples
    ]


def get_values(sv_id, data, source_samples, output_samples):
    if not sv_id or sv_id == ".":
        return None
    values = data.get(sv_id)
    if values is None:
        return ["NA"] * len(output_samples)
    if source_samples == output_samples:
        return values
    return align_values(values, source_samples, output_samples)


def combine_values(refa_values, refb_values, sample_count):
    if refa_values is None:
        refa_values = [None] * sample_count
    if refb_values is None:
        refb_values = [None] * sample_count

    out = []
    for refa_value, refb_value in zip(refa_values, refb_values):
        values = []
        if refa_value is not None:
            values.append(refa_value)
        if refb_value is not None:
            values.append(refb_value)

        if any(value == "1" for value in values):
            out.append("1")
        elif any(value == "0" for value in values):
            out.append("0")
        else:
            out.append("NA")
    return out


def write_stats(path, samples, row_count, counters, type_counts):
    all_types = sorted({svtype for sample_counts in type_counts.values() for svtype in sample_counts})
    with open(path, "w", encoding="utf-8") as out:
        out.write("sample\tpresent\tabsent\tmissing\tmissing_rate\tsingleton_present")
        for svtype in all_types:
            out.write("\tSVTYPE_" + svtype)
        out.write("\n")
        denominator = float(row_count) if row_count else 1.0
        for sample in samples:
            counts = counters[sample]
            out.write(
                f"{sample}\t{counts['present']}\t{counts['absent']}\t"
                f"{counts['missing']}\t{counts['missing'] / denominator:.6f}\t"
                f"{counts['singleton_present']}"
            )
            for svtype in all_types:
                out.write("\t" + str(type_counts[sample][svtype]))
            out.write("\n")


def main():
    args = parse_args()
    refa_samples, refa_pav = load_pav(args.refa_matrix)
    refb_samples, refb_pav = load_pav(args.refb_matrix)
    if set(refa_samples) != set(refb_samples):
        raise SystemExit("RefA and RefB PAV matrices have different sample sets")

    samples = refa_samples
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    tmp_out = args.out + ".tmp"
    tmp_stats = args.stats + ".tmp"

    counters = {sample: Counter() for sample in samples}
    type_counts = {sample: Counter() for sample in samples}
    row_count = 0

    with open(args.events, encoding="utf-8", errors="replace") as events, open(
        tmp_out, "w", encoding="utf-8"
    ) as out:
        header = events.readline().rstrip("\n").split("\t")
        column = {name: i for i, name in enumerate(header)}
        required = [
            "MetaSV_ID",
            "RefA_SV_ID",
            "RefB_SV_ID",
            "SVTYPE",
            "SVLEN",
            "match_status",
            "match_method",
            "confidence",
        ]
        missing = [name for name in required if name not in column]
        if missing:
            raise SystemExit("Missing columns in events file: " + ",".join(missing))

        out.write(
            "MetaSV_ID\tRefA_SV_ID\tRefB_SV_ID\tSVTYPE\tSVLEN\tmatch_status\t"
            "match_method\tconfidence\t"
            + "\t".join(samples)
            + "\n"
        )

        for line in events:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")

            def field(name):
                i = column[name]
                return row[i] if i < len(row) else "."

            meta_id = field("MetaSV_ID")
            refa_id = field("RefA_SV_ID")
            refb_id = field("RefB_SV_ID")
            svtype = field("SVTYPE")
            svlen = field("SVLEN")

            values = combine_values(
                get_values(refa_id, refa_pav, refa_samples, samples),
                get_values(refb_id, refb_pav, refb_samples, samples),
                len(samples),
            )
            out.write(
                "\t".join(
                    [
                        meta_id,
                        refa_id,
                        refb_id,
                        svtype,
                        svlen,
                        field("match_status"),
                        field("match_method"),
                        field("confidence"),
                    ]
                    + values
                )
                + "\n"
            )
            row_count += 1

            present_samples = []
            for sample, value in zip(samples, values):
                if value == "1":
                    counters[sample]["present"] += 1
                    type_counts[sample][svtype] += 1
                    present_samples.append(sample)
                elif value == "0":
                    counters[sample]["absent"] += 1
                else:
                    counters[sample]["missing"] += 1
            if len(present_samples) == 1:
                counters[present_samples[0]]["singleton_present"] += 1

    write_stats(tmp_stats, samples, row_count, counters, type_counts)
    os.replace(tmp_out, args.out)
    os.replace(tmp_stats, args.stats)


if __name__ == "__main__":
    main()
