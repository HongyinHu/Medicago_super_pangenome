#!/usr/bin/env python3
import collections
import csv
import glob
import gzip
import os


GWAS = "path/to/project/N_4.pod_spiny/28_GWAS_spiny"
SV_BASE = (
    "path/to/project/N_4.pod_spiny/"
    "09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/"
    "results/consensus_per_sample"
)
PHENO = os.path.join(GWAS, "01_input_audit/summary/phenotype_144.tsv")
OUTDIR = os.path.join(GWAS, "01_input_audit/summary")
TARGET_CHROM = "Chr4"
GENE_START = 88_979_818
GENE_END = 88_982_331
WINDOW_START = GENE_START - 10_000
WINDOW_END = GENE_END + 10_000


def parse_info(text):
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        elif item:
            result[item] = True
    return result


phenotypes = {}
with open(PHENO, encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        phenotypes[row["Sample_ID"]] = row

size_bins = [50, 500, 5_000, 50_000, 100_000, 1_000_000]
counts = collections.Counter()
support_counts = collections.Counter()
per_sample_filtered = collections.Counter()
local_records = []

for path in sorted(glob.glob(os.path.join(SV_BASE, "*", "*.vcf.gz"))):
    sample_id = os.path.basename(os.path.dirname(path))
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = parse_info(fields[7])
            svtype = info.get("SVTYPE", "NA")
            support = int(info.get("SUPP_EXT", info.get("SUPP", 0)) or 0)
            support_counts[(svtype, support)] += 1
            pos = int(fields[1])
            end_text = str(info.get("END", pos))
            end = int(end_text) if end_text.lstrip("-").isdigit() else pos
            length_text = str(info.get("SVLEN", ""))
            if length_text and length_text.split(",")[0].lstrip("-").isdigit():
                length = abs(int(length_text.split(",")[0]))
            elif svtype not in {"BND", "TRA"}:
                length = abs(end - pos)
            else:
                length = None

            counts[(svtype, "raw")] += 1
            if length is not None:
                for upper in size_bins:
                    if 50 <= length <= upper:
                        counts[(svtype, f"50_to_{upper}")] += 1
            clean_interval = (
                svtype in {"DEL", "INS", "DUP", "INV"}
                and length is not None
                and 50 <= length <= 1_000_000
                and support >= 2
                and fields[6] in {"PASS", "."}
            )
            if clean_interval:
                counts[(svtype, "clean_interval")] += 1
                per_sample_filtered[sample_id] += 1

            local_breakpoint = (
                fields[0] == TARGET_CHROM
                and length is not None
                and length <= 50_000
                and (
                    WINDOW_START <= pos <= WINDOW_END
                    or WINDOW_START <= end <= WINDOW_END
                )
            )
            if local_breakpoint:
                phenotype = phenotypes[sample_id]
                local_records.append(
                    {
                        "sample_id": sample_id,
                        "pod_spine_binary": phenotype[
                            "Pod_spine_binary_spiny1_spineless0"
                        ],
                        "pod_spine": phenotype["Pod_spine"],
                        "latin_name": phenotype["Latin_name"],
                        "chrom": fields[0],
                        "pos": pos,
                        "end": end,
                        "svtype": svtype,
                        "svlen_abs": length,
                        "supp": support,
                        "supp_vec": info.get("SUPP_VEC_EXT", info.get("SUPP_VEC", ".")),
                        "idlist": info.get("IDLIST_EXT", info.get("IDLIST", ".")),
                        "record_id": fields[2],
                    }
                )

with open(os.path.join(OUTDIR, "sv_qc_counts.tsv"), "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["svtype", "filter", "record_count"])
    for (svtype, filter_name), value in sorted(counts.items()):
        writer.writerow([svtype, filter_name, value])

with open(os.path.join(OUTDIR, "sv_support_counts.tsv"), "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["svtype", "caller_support", "record_count"])
    for (svtype, support), value in sorted(support_counts.items()):
        writer.writerow([svtype, support, value])

with open(
    os.path.join(OUTDIR, "sv_clean_counts_per_sample.tsv"), "w", newline=""
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample_id", "clean_interval_sv_count"])
    for sample_id in sorted(phenotypes):
        writer.writerow([sample_id, per_sample_filtered[sample_id]])

local_headers = [
    "sample_id",
    "pod_spine_binary",
    "pod_spine",
    "latin_name",
    "chrom",
    "pos",
    "end",
    "svtype",
    "svlen_abs",
    "supp",
    "supp_vec",
    "idlist",
    "record_id",
]
with open(os.path.join(OUTDIR, "Chr23997_local_calls.tsv"), "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=local_headers, delimiter="\t")
    writer.writeheader()
    writer.writerows(sorted(local_records, key=lambda row: (row["pos"], row["sample_id"])))

print(f"clean_interval_records={sum(per_sample_filtered.values())}")
print(f"local_breakpoint_records={len(local_records)}")
print(f"local_breakpoint_samples={len({row['sample_id'] for row in local_records})}")
