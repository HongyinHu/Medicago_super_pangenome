#!/usr/bin/env python3
import collections
import glob
import gzip
import os


SV_BASE = (
    "path/to/project/N_4.pod_spiny/"
    "09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/"
    "results/consensus_per_sample"
)
TARGET_CHROM = "Chr4"
TARGET_START = 88_969_818
TARGET_END = 88_992_331


def parse_info(value):
    info = {}
    for item in value.split(";"):
        if "=" in item:
            key, item_value = item.split("=", 1)
            info[key] = item_value
        elif item:
            info[item] = True
    return info


files = sorted(glob.glob(os.path.join(SV_BASE, "*", "*.vcf.gz")))
type_counts = collections.Counter()
caller_counts = collections.Counter()
gt_counts = collections.Counter()
sample_counts = {}
target_records = []

for path in files:
    sample_id = os.path.basename(os.path.dirname(path))
    record_count = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            record_count += 1
            fields = line.rstrip("\n").split("\t")
            info = parse_info(fields[7])
            svtype = info.get("SVTYPE", "NA")
            type_counts[svtype] += 1
            caller_counts[info.get("CALLERS", "NA")] += 1
            if len(fields) > 9:
                sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
                gt_counts[sample.get("GT", "NA")] += 1

            pos = int(fields[1])
            end_text = str(info.get("END", pos))
            end = int(end_text) if end_text.lstrip("-").isdigit() else pos
            if (
                fields[0] == TARGET_CHROM
                and max(pos, end) >= TARGET_START
                and min(pos, end) <= TARGET_END
            ):
                target_records.append(
                    (
                        sample_id,
                        fields[0],
                        pos,
                        end,
                        svtype,
                        info.get("SVLEN", "."),
                        info.get("CALLERS", "."),
                        fields[2],
                    )
                )
    sample_counts[sample_id] = record_count

counts = sorted(sample_counts.values())
print(f"files\t{len(files)}")
print(f"total_records\t{sum(counts)}")
if counts:
    print(f"per_sample_min\t{counts[0]}")
    print(f"per_sample_median\t{counts[len(counts) // 2]}")
    print(f"per_sample_max\t{counts[-1]}")
print("svtypes\t" + ",".join(f"{key}:{value}" for key, value in sorted(type_counts.items())))
print("gt\t" + ",".join(f"{key}:{value}" for key, value in sorted(gt_counts.items())))
print("caller_combos")
for key, value in caller_counts.most_common():
    print(f"{key}\t{value}")
print(f"chr23997_overlaps\t{len(target_records)}")
print(f"chr23997_samples\t{len({row[0] for row in target_records})}")
print("sample_id\tchrom\tpos\tend\tsvtype\tsvlen\tcallers\trecord_id")
for row in target_records:
    print("\t".join(map(str, row)))
