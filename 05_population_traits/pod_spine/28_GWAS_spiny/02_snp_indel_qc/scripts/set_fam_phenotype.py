#!/usr/bin/env python3
import csv
import sys


if len(sys.argv) != 3:
    raise SystemExit("Usage: set_fam_phenotype.py phenotype_144.tsv input.fam")

phenotype_path, fam_path = sys.argv[1:]
phenotypes = {}
with open(phenotype_path, encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        phenotypes[row["Sample_ID"]] = row[
            "Pod_spine_binary_spiny1_spineless0"
        ]

output_lines = []
with open(fam_path, encoding="utf-8") as handle:
    for line in handle:
        fields = line.rstrip("\n").split()
        if len(fields) < 6:
            raise RuntimeError(f"Malformed FAM line: {line.rstrip()}")
        sample_id = fields[1]
        if sample_id not in phenotypes:
            raise RuntimeError(f"No phenotype for FAM sample {sample_id}")
        # PLINK reserves 0 for a missing phenotype. Store the binary trait in
        # standard case/control coding here: 1 = spineless, 2 = spiny.
        fields[5] = str(int(phenotypes[sample_id]) + 1)
        output_lines.append(" ".join(fields))

with open(fam_path, "w", encoding="utf-8") as handle:
    handle.write("\n".join(output_lines) + "\n")

print(f"Updated {len(output_lines)} FAM phenotypes in {fam_path}")
