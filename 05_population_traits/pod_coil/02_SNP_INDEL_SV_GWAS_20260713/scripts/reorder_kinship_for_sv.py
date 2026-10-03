#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

import numpy as np


if len(sys.argv) != 8:
    raise SystemExit(
        "Usage: reorder_kinship_for_sv.py source.fam source.kinship target.fam "
        "strict_pheno.tsv eigenvec outdir n_pcs"
    )

source_fam, kinship_path, target_fam, phenotype_path, eigenvec_path, outdir, n_pcs_text = sys.argv[1:]
n_pcs = int(n_pcs_text)
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)


def read_fam(path):
    samples = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                raise RuntimeError(f"Malformed FAM line in {path}: {line.rstrip()}")
            samples.append((fields[0], fields[1]))
    return samples


source = read_fam(source_fam)
target = read_fam(target_fam)
source_index = {iid: i for i, (_, iid) in enumerate(source)}
missing = [iid for _, iid in target if iid not in source_index]
if missing:
    raise RuntimeError(f"SV samples absent from SNP kinship: {missing}")

matrix = np.loadtxt(kinship_path)
if matrix.shape != (len(source), len(source)):
    raise RuntimeError(f"Kinship shape {matrix.shape} != {(len(source), len(source))}")
indices = [source_index[iid] for _, iid in target]
reordered = matrix[np.ix_(indices, indices)]
np.savetxt(out / "SV_kinship.reordered.cXX.txt", reordered, fmt="%.10g")

phenotypes = {}
with open(phenotype_path, encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        phenotypes[row["IID"]] = int(row["Pod_spiral_binary_spiral1_nonspiral0"])

pcs = {}
with open(eigenvec_path, encoding="utf-8") as handle:
    for line in handle:
        fields = line.split()
        pcs[fields[1]] = [float(x) for x in fields[2 : 2 + n_pcs]]

with open(out / "SV.phenotype.0_1.txt", "w", encoding="utf-8") as pheno, open(
    out / f"SV.covariates.pc{n_pcs}.txt", "w", encoding="utf-8"
) as covar:
    for _, iid in target:
        pheno.write(f"{phenotypes[iid]}\n")
        covar.write("1 " + " ".join(str(x) for x in pcs[iid]) + "\n")

print(f"Reordered SNP kinship and covariates for {len(target)} SV samples")
