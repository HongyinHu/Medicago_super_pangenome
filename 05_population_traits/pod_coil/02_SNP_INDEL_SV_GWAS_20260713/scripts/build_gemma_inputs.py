#!/usr/bin/env python3
import csv
import sys
from pathlib import Path


if len(sys.argv) != 6:
    raise SystemExit(
        "Usage: build_gemma_inputs.py fam strict_pheno.tsv eigenvec outdir n_pcs"
    )

fam_path, phenotype_path, eigenvec_path, outdir, n_pcs_text = sys.argv[1:]
n_pcs = int(n_pcs_text)
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)

phenotypes = {}
with open(phenotype_path, encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        phenotypes[row["IID"]] = int(row["Pod_spiral_binary_spiral1_nonspiral0"])

fam = []
with open(fam_path, encoding="utf-8") as handle:
    for line in handle:
        fields = line.split()
        if len(fields) < 2:
            raise RuntimeError(f"Malformed FAM line: {line.rstrip()}")
        fam.append((fields[0], fields[1]))

pcs = {}
with open(eigenvec_path, encoding="utf-8") as handle:
    for line in handle:
        fields = line.split()
        if len(fields) < 2 + n_pcs:
            raise RuntimeError(f"Malformed eigenvec line: {line.rstrip()}")
        pcs[fields[1]] = [float(x) for x in fields[2 : 2 + n_pcs]]

missing_pheno = [iid for _, iid in fam if iid not in phenotypes]
missing_pcs = [iid for _, iid in fam if iid not in pcs]
if missing_pheno or missing_pcs:
    raise RuntimeError(f"Missing phenotypes={missing_pheno}; missing PCs={missing_pcs}")

with open(out / "phenotype.0_1.txt", "w", encoding="utf-8") as pheno, open(
    out / f"covariates.pc{n_pcs}.txt", "w", encoding="utf-8"
) as covar, open(
    out / "phenotype_with_ids.tsv", "w", encoding="utf-8"
) as with_ids:
    with_ids.write("FID\tIID\tPod_spiral_binary\n")
    for fid, iid in fam:
        pheno.write(f"{phenotypes[iid]}\n")
        covar.write("1 " + " ".join(str(x) for x in pcs[iid]) + "\n")
        with_ids.write(f"{fid}\t{iid}\t{phenotypes[iid]}\n")

counts = {0: 0, 1: 0}
for _, iid in fam:
    counts[phenotypes[iid]] += 1
if len(fam) != 125 or counts != {0: 40, 1: 85}:
    raise RuntimeError(f"Unexpected GEMMA sample composition: n={len(fam)}, {counts}")

print(f"GEMMA inputs: n={len(fam)}, non-spiral={counts[0]}, spiral={counts[1]}")
