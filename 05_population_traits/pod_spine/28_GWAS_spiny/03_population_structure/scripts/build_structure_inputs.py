#!/usr/bin/env python3
import csv
import sys
from pathlib import Path


if len(sys.argv) != 6:
    raise SystemExit(
        "Usage: build_structure_inputs.py fam phenotype.tsv eigenvec outdir n_pcs"
    )

fam_path, phenotype_path, eigenvec_path, outdir, n_pcs_text = sys.argv[1:]
n_pcs = int(n_pcs_text)
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)

phenotypes = {}
with open(phenotype_path, encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        phenotypes[row["Sample_ID"]] = int(
            row["Pod_spine_binary_spiny1_spineless0"]
        )

fam_samples = []
with open(fam_path, encoding="utf-8") as handle:
    for line in handle:
        fields = line.split()
        if len(fields) < 6:
            raise RuntimeError(f"Malformed FAM line: {line.rstrip()}")
        fam_samples.append((fields[0], fields[1]))

pcs = {}
with open(eigenvec_path, encoding="utf-8") as handle:
    for line in handle:
        fields = line.split()
        if len(fields) < 2 + n_pcs:
            raise RuntimeError(f"Malformed eigenvec line: {line.rstrip()}")
        pcs[fields[1]] = [float(value) for value in fields[2 : 2 + n_pcs]]

missing_pheno = [iid for _, iid in fam_samples if iid not in phenotypes]
missing_pcs = [iid for _, iid in fam_samples if iid not in pcs]
if missing_pheno or missing_pcs:
    raise RuntimeError(
        f"Missing phenotype={missing_pheno}; missing PCs={missing_pcs}"
    )

with open(out / "analysis_samples.keep", "w", encoding="utf-8") as keep, open(
    out / "phenotype.0_1.txt", "w", encoding="utf-8"
) as pheno, open(
    out / "phenotype_with_ids.tsv", "w", encoding="utf-8"
) as pheno_ids, open(
    out / f"covariates.pc{n_pcs}.txt", "w", encoding="utf-8"
) as covar:
    pheno_ids.write("FID\tIID\tPod_spine_binary\n")
    for fid, iid in fam_samples:
        value = phenotypes[iid]
        keep.write(f"{fid}\t{iid}\n")
        pheno.write(f"{value}\n")
        pheno_ids.write(f"{fid}\t{iid}\t{value}\n")
        covar.write("1 " + " ".join(str(x) for x in pcs[iid]) + "\n")

counts = {0: 0, 1: 0}
for _, iid in fam_samples:
    counts[phenotypes[iid]] += 1
with open(out / "analysis_sample_counts.tsv", "w", encoding="utf-8") as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"samples\t{len(fam_samples)}\n")
    handle.write(f"spineless_0\t{counts[0]}\n")
    handle.write(f"spiny_1\t{counts[1]}\n")
    handle.write(f"pcs\t{n_pcs}\n")

print(
    f"Prepared {len(fam_samples)} samples: "
    f"spineless={counts[0]}, spiny={counts[1]}, PCs={n_pcs}"
)
