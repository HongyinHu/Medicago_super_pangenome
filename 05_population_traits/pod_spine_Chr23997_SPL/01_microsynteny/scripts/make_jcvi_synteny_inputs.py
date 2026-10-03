#!/usr/bin/env python3
import csv
from pathlib import Path


ROOT = Path("path/to/project/41.Chr23997_realName")
SRC = ROOT / "01_microsynteny"
OUT = ROOT / "01_microsynteny" / "02_jcvi_plot"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = read_tsv(SRC / "neighboring_gene_table.tsv")
    blast = read_tsv(SRC / "reciprocal_blast_results.tsv")

    # Prefix seqids so the same chromosome name in Msa/R108 does not collide in JCVI.
    seqid_prefix = {"Msa": "Msa_Chr4", "Ath": "Ath_Chr1", "R108": "R108_Chr4"}
    species_order = {"Msa": 0, "Ath": 1, "R108": 2}

    with open(OUT / "all.bed", "w") as out:
        for g in sorted(genes, key=lambda x: (species_order[x["Species"]], int(x["Start"]))):
            out.write(
                "\t".join(
                    [
                        seqid_prefix[g["Species"]],
                        g["Start"],
                        g["End"],
                        g["Gene_ID"],
                        "0",
                        g["Strand"],
                    ]
                )
                + "\n"
            )

    with open(OUT / "layout.csv", "w") as out:
        out.write("#x, y, rotation, ha, va, color, ratio, label\n")
        out.write("0.54, 0.78, 0, left, top, #4c78a8, 1, Msa Chr4\n")
        out.write("0.54, 0.52, 0, left, bottom, #222222, 1, AtSPL10 Chr1\n")
        out.write("0.54, 0.25, 0, left, bottom, #4c78a8, 1, R108 Chr4\n")
        out.write("e, 0, 1, gainsboro\n")
        out.write("e, 1, 2, gainsboro\n")

    rows = []
    link_rows = []

    def add_link(prefix, row, category, note):
        rows.append((prefix, row))
        link_rows.append(
            {
                "Msa_gene": row[0],
                "Ath_gene": row[1],
                "R108_gene": row[2],
                "Category": category,
                "Note": note,
            }
        )

    # Red target comparison line: these are the focal loci, not asserted true RBH.
    add_link("r*", ["Chr23997", "AT1G27370", "R10821044"], "focal_candidate", "candidate SPL loci")

    # Draw every local-window BLASTP top hit for Msa-Ath and Ath-R108.
    # Green = reciprocal local top hit; grey = non-reciprocal but still within both local windows.
    for r in blast:
        qs, ss = r["Query_species"], r["Subject_species"]
        if r["Subject_in_local_window"] != "1":
            continue
        if set([qs, ss]) not in (set(["Msa", "Ath"]), set(["Ath", "R108"])):
            continue
        if not r["Subject_gene"]:
            continue
        row = [".", ".", "."]
        row[species_order[qs]] = r["Query_gene"]
        row[species_order[ss]] = r["Subject_gene"]
        if r["Reciprocal_top_hit"] == "1":
            prefix = "g*"
            category = "local_RBH"
        else:
            prefix = "#c7c7c7*"
            category = "local_top_hit_non_RBH"
        # Highlight the SPL11 ambiguity separately.
        if row[0] == "Chr23997" and row[1] == "AT1G27360":
            prefix = "#f28e2b*"
            category = "SPL11_local_top_hit"
        add_link(
            prefix,
            row,
            category,
            "pident={}; evalue={}; bitscore={}".format(r["Pident"], r["Evalue"], r["Bitscore"]),
        )

    # Singleton rows keep the full local windows visible without drawing direct Msa-R108 links.
    for sp in ("Msa", "Ath", "R108"):
        col = species_order[sp]
        for g in genes:
            if g["Species"] != sp:
                continue
            row = [".", ".", "."]
            row[col] = g["Gene_ID"]
            rows.append(("", row))

    seen = set()
    with open(OUT / "mcscan.txt", "w") as out:
        for prefix, row in rows:
            key = tuple(row)
            if key in seen:
                continue
            seen.add(key)
            out.write(prefix + "\t".join(row) + "\n")

    with open(OUT / "local_window_links.tsv", "w", newline="") as out:
        fields = ["Msa_gene", "Ath_gene", "R108_gene", "Category", "Note"]
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        emitted = set()
        for r in link_rows:
            key = (r["Msa_gene"], r["Ath_gene"], r["R108_gene"], r["Category"])
            if key in emitted:
                continue
            emitted.add(key)
            w.writerow(r)

    with open(OUT / "README.txt", "w") as out:
        out.write("JCVI microsynteny plot inputs for Chr23997 vs Arabidopsis SPL10 and R108.\n")
        out.write("Red: focal SPL candidate comparison. Green: local reciprocal top BLASTP hit. Grey: local top BLASTP hit that is not reciprocal. Orange: Chr23997-to-SPL11 local top-hit ambiguity. No direct Msa-R108 synteny edge is drawn.\n")
        out.write("Generated from neighboring_gene_table.tsv and reciprocal_blast_results.tsv.\n")


if __name__ == "__main__":
    main()
