#!/usr/bin/env python3
import csv
import pathlib
import collections

base = pathlib.Path("path/to/project/39.TE_type_soloLTR/output")

print("Top summary:")
with open(base / "te_subtype_percent.tsv", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["subtype"] in {"Unclassified LTR-RTs", "Unknown"}:
            print(
                row["species"],
                row["subtype"],
                "TE%",
                row["percent_of_all_TE"],
                "Genome%",
                row["percent_of_genome"],
                "count",
                row["annotation_count"],
            )

print("\nTEsorter clade unknown / total LTR rows:")
for cls in sorted((base / "tesorter").glob("*/*.cls.tsv")):
    species = cls.parent.name
    total = ltr = unk = super_only = 0
    clades = collections.Counter()
    superfamilies = collections.Counter()
    with open(cls, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            if row["Order"] == "LTR":
                ltr += 1
                clade = row["Clade"] or ""
                superfamily = row["Superfamily"] or ""
                clades[clade] += 1
                superfamilies[superfamily] += 1
                if clade.lower() in {"", "?", "unknown", "none", "unclassified"}:
                    unk += 1
                    if superfamily in {"Copia", "Gypsy"}:
                        super_only += 1
    print(
        species,
        "LTR_rows",
        ltr,
        "unknown_clade",
        unk,
        "superfamily_only",
        super_only,
        "total_rows",
        total,
        "top_clades",
        clades.most_common(5),
        "superfamilies",
        superfamilies.most_common(),
    )
