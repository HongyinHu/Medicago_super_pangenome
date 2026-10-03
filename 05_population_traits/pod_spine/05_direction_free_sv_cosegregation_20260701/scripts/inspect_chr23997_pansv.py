#!/usr/bin/env python3
from pathlib import Path
import subprocess

TABIX = "path/to/home/anaconda3/bin/tabix"
REGION = "Chr4:88978818-88982780"
FILES = [
    "path/to/project/N_3.call_SV/01.read_based_dualref_hifi/05_dualref_panSV_20260626/results/RefA/panSV/panSV.RefA.discovery.vcf.gz",
    "path/to/project/N_3.call_SV/01.read_based_dualref_hifi/05_dualref_panSV_20260626/results/RefA/panSV/panSV.RefA.genotyped.18species.vcf.gz",
    "path/to/project/N_3.call_SV/01.read_based_dualref_hifi/05_dualref_panSV_20260626/results/RefA/panSV/panSV.RefA.PAV.matrix.tsv",
]

for f in FILES:
    path = Path(f)
    print("##", f)
    if not path.exists():
        print("MISSING")
        continue
    if f.endswith(".vcf.gz"):
        proc = subprocess.run([TABIX, f, REGION], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        lines = [x for x in proc.stdout.splitlines() if x and not x.startswith("#")]
        print("records", len(lines))
        for line in lines[:50]:
            parts = line.split("\t")
            print("\t".join(parts[:8]))
    else:
        # PAV matrix may be large and not indexed; search exact region token only lightly.
        print("exists", path.stat().st_size)
