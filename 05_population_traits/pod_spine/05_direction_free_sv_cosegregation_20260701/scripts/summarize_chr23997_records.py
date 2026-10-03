#!/usr/bin/env python3
from pathlib import Path
import csv
from collections import Counter

path = Path("path/to/project/N_4.pod_spiny/05_direction_free_sv_cosegregation_20260701/summary/Chr23997_caller_records_20260701.tsv")
rows = []
with path.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("near_manual_del") == "True":
            rows.append(row)

print("near_manual_DEL_records", len(rows))
print("sample\tphenotype\tcaller\tpos\tend\tsvlen\tfilter")
for r in sorted(rows, key=lambda x: (x["sample"], x["caller"], int(x["pos"]))):
    print("\t".join([r["sample"], r["phenotype"], r["caller"], r["pos"], r["end"], r["svlen"], r["filter"]]))

print("\ncluster_by_rounded_breakpoints")
c = Counter((round(int(r["pos"]) / 50) * 50, round(int(r["end"]) / 50) * 50, r["svlen"], r["caller"]) for r in rows)
for key, value in c.most_common(80):
    print("\t".join(map(str, key)), value, sep="\t")
