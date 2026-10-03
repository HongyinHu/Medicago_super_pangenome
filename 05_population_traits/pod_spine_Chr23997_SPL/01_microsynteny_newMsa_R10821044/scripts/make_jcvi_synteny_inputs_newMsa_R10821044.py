#!/usr/bin/env python3
import csv
from pathlib import Path


ROOT = Path("path/to/project/41.Chr23997_realName")
SRC = ROOT / "01_microsynteny_newMsa_R10821044"
OUT = SRC / "02_jcvi_plot"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = read_tsv(SRC / "neighboring_gene_table.tsv")
    blast = read_tsv(SRC / "reciprocal_blast_results.tsv")

    targets = {g["Species"]: g["Gene_ID"] for g in genes if g["Is_target"] == "1"}
    species_order = {"Msa": 0, "Ath": 1, "R108": 2}

    # Prefix sequence IDs so shared chromosome names do not collide in JCVI.
    seqid_prefix = {
        "Msa": "newMsa_" + next(g["Chromosome"] for g in genes if g["Species"] == "Msa"),
        "Ath": "Ath_" + next(g["Chromosome"] for g in genes if g["Species"] == "Ath"),
        "R108": "R108_" + next(g["Chromosome"] for g in genes if g["Species"] == "R108"),
    }

    with open(OUT / "all.bed", "w") as out:
        for g in sorted(genes, key=lambda x: (species_order[x["Species"]], int(x["Start"]))):
            out.write("\t".join([
                seqid_prefix[g["Species"]],
                g["Start"],
                g["End"],
                g["Gene_ID"],
                "0",
                g["Strand"],
            ]) + "\n")

    with open(OUT / "layout.csv", "w") as out:
        out.write("#x, y, rotation, ha, va, color, ratio, label\n")
        out.write(f"0.54, 0.78, 0, left, top, #4c78a8, 1, new Msa {seqid_prefix['Msa'].replace('newMsa_', '')}\n")
        out.write(f"0.54, 0.52, 0, left, bottom, #222222, 1, Arabidopsis SPL10 {seqid_prefix['Ath'].replace('Ath_', '')}\n")
        out.write(f"0.54, 0.25, 0, left, bottom, #4c78a8, 1, R108 {seqid_prefix['R108'].replace('R108_', '')}\n")
        out.write("e, 0, 1, gainsboro\n")
        out.write("e, 1, 2, gainsboro\n")

    rows = []
    link_rows = []

    def add_link(prefix, row, category, note):
        rows.append((prefix, row))
        link_rows.append({
            "newMsa_gene": row[0],
            "Ath_gene": row[1],
            "R108_gene": row[2],
            "Category": category,
            "Note": note,
        })

    add_link(
        "r*",
        [targets["Msa"], targets["Ath"], targets["R108"]],
        "focal_candidate",
        "R10821044-anchored focal SPL candidate comparison",
    )

    # Draw every local-window BLASTP top hit for newMsa-Ath and Ath-R108.
    # Green: reciprocal local top hit. Grey: local top hit that is not reciprocal.
    for r in blast:
        qs, ss = r["Query_species"], r["Subject_species"]
        if r["Subject_in_local_window"] != "1":
            continue
        if frozenset([qs, ss]) not in (frozenset(["Msa", "Ath"]), frozenset(["Ath", "R108"])):
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
        add_link(
            prefix,
            row,
            category,
            "pident={}; evalue={}; bitscore={}".format(r["Pident"], r["Evalue"], r["Bitscore"]),
        )

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
        fields = ["newMsa_gene", "Ath_gene", "R108_gene", "Category", "Note"]
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        emitted = set()
        for r in link_rows:
            key = (r["newMsa_gene"], r["Ath_gene"], r["R108_gene"], r["Category"])
            if key in emitted:
                continue
            emitted.add(key)
            w.writerow(r)

    with open(OUT / "README.txt", "w") as out:
        out.write("JCVI microsynteny plot inputs for new Msa target anchored by R10821044, Arabidopsis SPL10, and R10821044.\n")
        out.write("Red: focal SPL candidate comparison. Green: local reciprocal top BLASTP hit. Grey: local top BLASTP hit that is not reciprocal. No direct newMsa-R108 synteny edge is drawn.\n")
        out.write("Generated from neighboring_gene_table.tsv and reciprocal_blast_results.tsv.\n")


if __name__ == "__main__":
    main()
