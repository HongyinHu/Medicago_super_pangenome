#!/usr/bin/env python3
import csv
import sys
from collections import Counter
from pathlib import Path


if len(sys.argv) != 7:
    raise SystemExit(
        "Usage: make_model_metadata.py fam phenotype.tsv eigenvec outdir label n_pcs"
    )

fam_path, phenotype_path, eigenvec_path, outdir, label, n_pcs_text = sys.argv[1:]
n_pcs = int(n_pcs_text)
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)

phenotypes = {}
with open(phenotype_path, encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        phenotypes[row["IID"]] = int(row["coil_ge1_binary_spiral1_nonspiral0"])

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
        pcs[fields[1]] = [float(value) for value in fields[2 : 2 + n_pcs]]

missing_pheno = [iid for _, iid in fam if iid not in phenotypes]
missing_pcs = [iid for _, iid in fam if iid not in pcs]
if missing_pheno or missing_pcs:
    raise RuntimeError(f"Missing phenotypes={missing_pheno}; missing PCs={missing_pcs}")

model_path = out / f"{label}.model.tsv"
covar_path = out / f"{label}.plink.covar.tsv"
pheno_path = out / f"{label}.plink.pheno.tsv"

with open(model_path, "w", encoding="utf-8") as model, open(
    covar_path, "w", encoding="utf-8"
) as covar, open(pheno_path, "w", encoding="utf-8") as pheno:
    pc_names = [f"PC{i}" for i in range(1, n_pcs + 1)]
    model.write("id\ty\t" + "\t".join(pc_names) + "\n")
    covar.write("FID\tIID\t" + "\t".join(pc_names) + "\n")
    for fid, iid in fam:
        pc_values = pcs[iid]
        model.write(f"{iid}\t{phenotypes[iid]}\t" + "\t".join(map(str, pc_values)) + "\n")
        covar.write(f"{fid}\t{iid}\t" + "\t".join(map(str, pc_values)) + "\n")
        pheno.write(f"{fid}\t{iid}\t{phenotypes[iid] + 1}\n")

counts = Counter(phenotypes[iid] for _, iid in fam)
with open(out / f"{label}.sample_audit.tsv", "w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"samples_after_marker_qc\t{len(fam)}\n")
    handle.write(f"spiral_case_1\t{counts[1]}\n")
    handle.write(f"nonspiral_control_0\t{counts[0]}\n")
    handle.write(f"PC_covariates\t{n_pcs}\n")

if len(fam) < 100 or counts[0] < 20 or counts[1] < 20:
    raise RuntimeError(f"Insufficient post-QC sample composition: n={len(fam)}, {dict(counts)}")

print(f"{label} metadata: n={len(fam)}, spiral={counts[1]}, non-spiral={counts[0]}")

