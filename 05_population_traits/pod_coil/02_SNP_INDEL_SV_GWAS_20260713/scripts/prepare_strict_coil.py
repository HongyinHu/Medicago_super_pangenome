#!/usr/bin/env python3
import csv
import sys
from pathlib import Path


if len(sys.argv) != 4:
    raise SystemExit("Usage: prepare_strict_coil.py phenotype.tsv source.fam outdir")

phenotype_path, fam_path, outdir = sys.argv[1:]
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)

rows = {}
with open(phenotype_path, encoding="utf-8-sig") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        rows[row["Sample_ID"]] = row

fam = []
with open(fam_path, encoding="utf-8") as handle:
    for line in handle:
        fields = line.split()
        if len(fields) < 2:
            raise RuntimeError(f"Malformed FAM line: {line.rstrip()}")
        fam.append((fields[0], fields[1]))

included = []
excluded = []
for fid, iid in fam:
    if iid not in rows:
        raise RuntimeError(f"Phenotype missing for genotype sample {iid}")
    row = rows[iid]
    state = row["Pod_spiral"].strip().lower()
    if state == "spiral":
        value = 1
    elif state == "non spiral":
        value = 0
    else:
        excluded.append((fid, iid, state, row.get("Latin_name", "")))
        continue
    included.append((fid, iid, value, row))

with open(out / "strict_samples.keep", "w", encoding="utf-8") as handle:
    for fid, iid, _, _ in included:
        handle.write(f"{fid}\t{iid}\n")

with open(out / "strict_sample_ids.txt", "w", encoding="utf-8") as handle:
    for _, iid, _, _ in included:
        handle.write(iid + "\n")

fields = [
    "FID",
    "IID",
    "Pod_spiral_binary_spiral1_nonspiral0",
    "Pod_spiral",
    "Latin_name",
    "Section",
    "Trait_priority",
    "Pod_trait_source",
]
with open(out / "coil_phenotype_strict.tsv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    for fid, iid, value, row in included:
        writer.writerow(
            {
                "FID": fid,
                "IID": iid,
                "Pod_spiral_binary_spiral1_nonspiral0": value,
                "Pod_spiral": row["Pod_spiral"],
                "Latin_name": row.get("Latin_name", ""),
                "Section": row.get("Section", ""),
                "Trait_priority": row.get("Trait_priority", ""),
                "Pod_trait_source": row.get("Pod_trait_source", ""),
            }
        )

with open(out / "excluded_weak_spiral.tsv", "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["FID", "IID", "Pod_spiral", "Latin_name"])
    writer.writerows(excluded)

counts = {0: 0, 1: 0}
for _, _, value, _ in included:
    counts[value] += 1
with open(out / "strict_sample_counts.tsv", "w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"source_genotype_samples\t{len(fam)}\n")
    handle.write(f"strict_samples\t{len(included)}\n")
    handle.write(f"non_spiral_0\t{counts[0]}\n")
    handle.write(f"spiral_1\t{counts[1]}\n")
    handle.write(f"excluded_weak_spiral\t{len(excluded)}\n")

if len(included) != 125 or counts != {0: 40, 1: 85}:
    raise RuntimeError(
        f"Unexpected strict phenotype composition: n={len(included)}, counts={counts}"
    )

print(
    f"Strict coil phenotype: n={len(included)}, non-spiral={counts[0]}, "
    f"spiral={counts[1]}, excluded weak={len(excluded)}"
)
