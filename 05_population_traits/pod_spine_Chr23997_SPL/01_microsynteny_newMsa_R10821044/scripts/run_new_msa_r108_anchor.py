#!/usr/bin/env python3
import argparse
import gzip
import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path


def xopen(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_attrs(s):
    out = {}
    for part in s.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
        else:
            continue
        out[k] = v
    return out


def norm_gene_id(raw):
    if raw is None:
        return None
    raw = raw.strip()
    for prefix in ("gene:", "transcript:", "CDS:"):
        if raw.startswith(prefix):
            raw = raw[len(prefix):]
    return raw


def parse_gff(path, species):
    genes = {}
    with xopen(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != "gene":
                continue
            a = parse_attrs(attrs)
            gid = norm_gene_id(a.get("gene_id") or a.get("ID") or a.get("Name"))
            if not gid:
                continue
            biotype = a.get("biotype") or a.get("gene_biotype") or ""
            if species == "Ath" and biotype and biotype != "protein_coding":
                continue
            genes[gid] = {
                "species": species,
                "gene_id": gid,
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "attrs": a,
            }
    return genes


def fasta_records(path):
    name = None
    desc = ""
    chunks = []
    with xopen(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks).replace(".", "")
                desc = line[1:]
                name = desc.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, desc, "".join(chunks).replace(".", "")


def gene_from_pep_header(seqid, desc):
    m = re.search(r"(?:^|\s)gene=([^\s;]+)", desc)
    if m:
        return norm_gene_id(m.group(1))
    m = re.search(r"(?:^|\s)gene:([^\s;]+)", desc)
    if m:
        return norm_gene_id(m.group(1))
    if "." in seqid:
        return seqid.split(".", 1)[0]
    return seqid


def parse_pep(path):
    proteins = {}
    gene_to_proteins = defaultdict(list)
    for seqid, desc, seq in fasta_records(path):
        gid = gene_from_pep_header(seqid, desc)
        proteins[seqid] = {"protein_id": seqid, "gene_id": gid, "desc": desc, "seq": seq}
        gene_to_proteins[gid].append(seqid)
    return proteins, gene_to_proteins


def longest_protein(gid, gene_to_proteins, proteins):
    ids = gene_to_proteins.get(gid, [])
    if not ids:
        return None
    return max(ids, key=lambda x: len(proteins[x]["seq"]))


def write_longest_proteome(genes, gene_to_proteins, proteins, out_fa):
    with open(out_fa, "w") as out:
        for gid in sorted(genes):
            pid = longest_protein(gid, gene_to_proteins, proteins)
            if not pid:
                continue
            seq = proteins[pid]["seq"]
            out.write(f">{pid} gene={gid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def write_gene_fasta(gene_ids, gene_to_proteins, proteins, out_fa):
    with open(out_fa, "w") as out:
        for gid in gene_ids:
            pid = longest_protein(gid, gene_to_proteins, proteins)
            if not pid:
                continue
            seq = proteins[pid]["seq"]
            out.write(f">{pid} gene={gid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def run(cmd, log_path):
    with open(log_path, "a") as log:
        log.write(" ".join(map(str, cmd)) + "\n")
    subprocess.run(cmd, check=True)


def make_db(makeblastdb, fasta, db_prefix, log_path):
    if Path(str(db_prefix) + ".pin").exists():
        return
    run([makeblastdb, "-in", str(fasta), "-dbtype", "prot", "-out", str(db_prefix)], log_path)


def blast_top(blastp, query, db_prefix, out_tsv, log_path, evalue="1e-5"):
    run([
        blastp, "-query", str(query), "-db", str(db_prefix), "-evalue", evalue,
        "-max_target_seqs", "1", "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-out", str(out_tsv), "-num_threads", "4",
    ], log_path)


def read_blast(path, protein_to_gene):
    hits = {}
    if not Path(path).exists():
        return hits
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            q, s = p[0], p[1]
            if q in hits:
                continue
            hits[q] = {
                "query_protein": q,
                "subject_protein": s,
                "subject_gene": protein_to_gene.get(s, s.split(".", 1)[0]),
                "pident": p[2],
                "length": p[3],
                "evalue": p[10],
                "bitscore": p[11],
            }
    return hits


def window(genes, gene_to_proteins, target_gene, flank=5):
    if target_gene not in genes:
        raise SystemExit(f"Target gene not found in GFF: {target_gene}")
    chrom = genes[target_gene]["chrom"]
    rows = [g for g in genes.values() if g["chrom"] == chrom and g["gene_id"] in gene_to_proteins]
    rows.sort(key=lambda x: (x["start"], x["end"], x["gene_id"]))
    idx = next(i for i, g in enumerate(rows) if g["gene_id"] == target_gene)
    sub = rows[max(0, idx - flank): idx + flank + 1]
    out = []
    for g in sub:
        rel = rows.index(g) - idx
        label = target_gene if rel == 0 else str(rel)
        out.append((rel, label, g))
    return out


def position(g):
    return f"{g['chrom']}:{g['start']}-{g['end']}"


def find_new_msa_target_from_r108(blastp, makeblastdb, out, maps, log_path):
    r108_db = out / "blastdb" / "R108.longest"
    msa_db = out / "blastdb" / "Msa.longest"
    make_db(makeblastdb, out / "R108.longest.pep.fa", r108_db, log_path)
    make_db(makeblastdb, out / "Msa.longest.pep.fa", msa_db, log_path)

    r108_query = out / "blast" / "R10821044.query.fa"
    write_gene_fasta(["R10821044"], maps["R108"]["gene_to_proteins"], maps["R108"]["proteins"], r108_query)
    top = out / "blast" / "R10821044_vs_newMsa.top.tsv"
    blast_top(blastp, r108_query, msa_db, top, log_path)
    hits = read_blast(top, maps["Msa"]["protein_to_gene"])
    if not hits:
        raise SystemExit("No new Msa BLASTP hit found for R10821044")
    hit = next(iter(hits.values()))
    msa_gene = hit["subject_gene"]

    msa_query = out / "blast" / "newMsa_target.query.fa"
    write_gene_fasta([msa_gene], maps["Msa"]["gene_to_proteins"], maps["Msa"]["proteins"], msa_query)
    back = out / "blast" / "newMsa_target_vs_R108.top.tsv"
    blast_top(blastp, msa_query, r108_db, back, log_path)
    back_hits = read_blast(back, maps["R108"]["protein_to_gene"])
    back_hit = next(iter(back_hits.values())) if back_hits else None
    reciprocal = bool(back_hit and back_hit["subject_gene"] == "R10821044")
    return msa_gene, reciprocal, hit, back_hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = Path(args.data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for d in ("blastdb", "blast", "logs"):
        (out / d).mkdir(exist_ok=True)
    log_path = out / "commands_log.txt"
    blastp = os.environ["BLASTP"]
    makeblastdb = os.environ["MAKEBLASTDB"]

    specs = {
        "Msa": {
            "label": "Medicago new Msa",
            "gff": data / "new_genome_Msa" / "genome_Msa.anno.gff",
            "pep": data / "new_genome_Msa" / "genome_Msa.anno.pep",
            "target": None,
            "target_label": None,
        },
        "Ath": {
            "label": "Arabidopsis thaliana",
            "gff": data / "Arabidopsis_thaliana.TAIR10.63.gff3.gz",
            "pep": data / "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz",
            "target": "AT1G27370",
            "target_label": "AtSPL10",
        },
        "R108": {
            "label": "Medicago R108",
            "gff": data / "genome_R108.gff",
            "pep": data / "genome_R108.pep",
            "target": "R10821044",
            "target_label": "R10821044",
        },
    }

    maps = {}
    for sp, s in specs.items():
        genes = parse_gff(s["gff"], sp)
        proteins, gene_to_proteins = parse_pep(s["pep"])
        protein_to_gene = {pid: rec["gene_id"] for pid, rec in proteins.items()}
        maps[sp] = {
            "genes": genes,
            "proteins": proteins,
            "gene_to_proteins": gene_to_proteins,
            "protein_to_gene": protein_to_gene,
        }
        write_longest_proteome(genes, gene_to_proteins, proteins, out / f"{sp}.longest.pep.fa")

    msa_gene, msa_is_rbh, msa_hit, back_hit = find_new_msa_target_from_r108(
        blastp, makeblastdb, out, maps, log_path
    )
    specs["Msa"]["target"] = msa_gene
    specs["Msa"]["target_label"] = msa_gene

    windows = {}
    neighbor_rows = []
    for sp in ("Msa", "Ath", "R108"):
        w = window(maps[sp]["genes"], maps[sp]["gene_to_proteins"], specs[sp]["target"], flank=5)
        windows[sp] = w
        write_gene_fasta(
            [g["gene_id"] for _, _, g in w],
            maps[sp]["gene_to_proteins"],
            maps[sp]["proteins"],
            out / f"{sp}.window.pep.fa",
        )
        for rel, label, g in w:
            gid = g["gene_id"]
            pid = longest_protein(gid, maps[sp]["gene_to_proteins"], maps[sp]["proteins"])
            display = specs[sp]["target_label"] if rel == 0 else gid
            neighbor_rows.append({
                "Species": sp,
                "Species_label": specs[sp]["label"],
                "Gene_order": label,
                "Order_index": str(rel + 6),
                "Gene_ID": gid,
                "Display_label": display,
                "Chromosome": g["chrom"],
                "Start": str(g["start"]),
                "End": str(g["end"]),
                "Position": position(g),
                "Strand": g["strand"],
                "Protein_sequence_ID": pid or "",
                "Is_target": "1" if rel == 0 else "0",
            })

    for sp in ("Msa", "Ath", "R108"):
        make_db(makeblastdb, out / f"{sp}.longest.pep.fa", out / "blastdb" / f"{sp}.longest", log_path)

    pair_jobs = [
        ("Msa", "Ath"), ("Ath", "Msa"),
        ("R108", "Ath"), ("Ath", "R108"),
        ("Msa", "R108"), ("R108", "Msa"),
    ]
    blast_maps = {}
    for qsp, ssp in pair_jobs:
        out_tsv = out / "blast" / f"{qsp}_window_vs_{ssp}.top.tsv"
        blast_top(blastp, out / f"{qsp}.window.pep.fa", out / "blastdb" / f"{ssp}.longest", out_tsv, log_path)
        blast_maps[(qsp, ssp)] = read_blast(out_tsv, maps[ssp]["protein_to_gene"])

    window_gene_sets = {sp: {g["gene_id"] for _, _, g in windows[sp]} for sp in windows}
    rows = []
    links = {}
    group_by_ath = {}
    palette_idx = 1
    for qsp, ssp in [
        ("Msa", "Ath"), ("R108", "Ath"),
        ("Ath", "Msa"), ("Ath", "R108"),
        ("Msa", "R108"), ("R108", "Msa"),
    ]:
        rev = blast_maps.get((ssp, qsp), {})
        for rel, label, g in windows[qsp]:
            qgene = g["gene_id"]
            qpid = longest_protein(qgene, maps[qsp]["gene_to_proteins"], maps[qsp]["proteins"])
            h = blast_maps[(qsp, ssp)].get(qpid)
            if not h:
                rows.append([qsp, qgene, qpid, ssp, "", "", "", "", "", "", "0", "0"])
                continue
            sgene = h["subject_gene"]
            spid = h["subject_protein"]
            rh = rev.get(spid)
            reciprocal = bool(rh and rh["subject_gene"] == qgene)
            in_window = sgene in window_gene_sets[ssp]
            rows.append([
                qsp, qgene, qpid, ssp, sgene, spid, h["pident"], h["length"],
                h["evalue"], h["bitscore"], "1" if reciprocal else "0", "1" if in_window else "0",
            ])
            if reciprocal and in_window and ((qsp == "Msa" and ssp == "Ath") or (qsp == "R108" and ssp == "Ath")):
                ath_gene = sgene
                if ath_gene not in group_by_ath:
                    group_by_ath[ath_gene] = f"OG{palette_idx:02d}"
                    palette_idx += 1
                og = group_by_ath[ath_gene]
                links[(qsp, qgene)] = (og, ath_gene)
                links[(ssp, sgene)] = (og, ath_gene)

    for sp in ("Msa", "Ath", "R108"):
        links[(sp, specs[sp]["target"])] = ("SPL_candidate", "AT1G27370")

    with open(out / "neighboring_gene_table.tsv", "w") as fh:
        fields = [
            "Species", "Species_label", "Gene_order", "Order_index", "Gene_ID", "Display_label",
            "Chromosome", "Start", "End", "Position", "Strand", "Protein_sequence_ID", "Is_target",
        ]
        fh.write("\t".join(fields) + "\n")
        for r in neighbor_rows:
            fh.write("\t".join(r[f] for f in fields) + "\n")

    with open(out / "reciprocal_blast_results.tsv", "w") as fh:
        fh.write("\t".join([
            "Query_species", "Query_gene", "Query_protein", "Subject_species", "Subject_gene",
            "Subject_protein", "Pident", "Alignment_length", "Evalue", "Bitscore",
            "Reciprocal_top_hit", "Subject_in_local_window",
        ]) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")

    with open(out / "synteny_links.tsv", "w") as fh:
        fh.write("Species\tGene_ID\tOrthogroup\tArabidopsis_anchor\n")
        for (sp, gid), (og, anchor) in sorted(links.items()):
            fh.write(f"{sp}\t{gid}\t{og}\t{anchor}\n")

    def count_conserved(sp):
        return len({gid for (s, gid), (og, anchor) in links.items() if s == sp and og != "SPL_candidate"})

    msa_conserved = count_conserved("Msa")
    r108_conserved = count_conserved("R108")
    msa_target_vs_ath = next((r for r in rows if r[0] == "Msa" and r[1] == msa_gene and r[3] == "Ath"), None)
    atspl10_vs_msa = next((r for r in rows if r[0] == "Ath" and r[1] == "AT1G27370" and r[3] == "Msa"), None)
    atspl10_vs_r108 = next((r for r in rows if r[0] == "Ath" and r[1] == "AT1G27370" and r[3] == "R108"), None)
    case = "B"
    conclusion = "ambiguous microsynteny; use SPL10-like transcription factor"
    if msa_conserved >= 2:
        case = "A"
        conclusion = f"{msa_gene} is supported as the new Msa local orthologous SPL candidate anchored by R10821044"

    with open(out / "interpretation_report.md", "w") as fh:
        fh.write("# New Msa R10821044-anchored microsynteny interpretation\n\n")
        fh.write("## Locus summary\n\n")
        for sp in ("Msa", "Ath", "R108"):
            g = maps[sp]["genes"][specs[sp]["target"]]
            pid = longest_protein(specs[sp]["target"], maps[sp]["gene_to_proteins"], maps[sp]["proteins"])
            fh.write(f"- {specs[sp]['label']}: {specs[sp]['target']} ({specs[sp]['target_label']}), "
                     f"{position(g)}, strand {g['strand']}, protein {pid}\n")
        fh.write(
            f"- New Msa target selection: R10821044 top BLASTP hit was {msa_gene}; "
            f"reciprocal back to R10821044: {msa_is_rbh}.\n"
        )
        if back_hit:
            fh.write(f"- Reciprocal check hit: {back_hit['subject_gene']} ({back_hit['subject_protein']}).\n")
        fh.write("\n## Microsynteny evidence\n\n")
        fh.write(f"- Conserved Arabidopsis-window RBH genes around new Msa target, excluding focal SPL gene: {msa_conserved}\n")
        fh.write(f"- Conserved Arabidopsis-window RBH genes around R10821044, excluding focal SPL gene: {r108_conserved}\n")
        if msa_target_vs_ath:
            fh.write(
                f"- New Msa target's top Arabidopsis hit was {msa_target_vs_ath[4]} "
                f"({msa_target_vs_ath[5]}); reciprocal top hit: {msa_target_vs_ath[10]}; "
                f"in local window: {msa_target_vs_ath[11]}.\n"
            )
        if atspl10_vs_msa:
            fh.write(
                f"- AtSPL10's top new Msa hit was {atspl10_vs_msa[4]} "
                f"({atspl10_vs_msa[5]}); reciprocal top hit: {atspl10_vs_msa[10]}; "
                f"in local window: {atspl10_vs_msa[11]}.\n"
            )
        if atspl10_vs_r108:
            fh.write(
                f"- AtSPL10's top R108 hit was {atspl10_vs_r108[4]} "
                f"({atspl10_vs_r108[5]}); reciprocal top hit: {atspl10_vs_r108[10]}; "
                f"in local window: {atspl10_vs_r108[11]}.\n"
            )
        fh.write("- Gene order was assessed in genomic coordinate order using five protein-coding genes on each side of the target locus.\n")
        fh.write("- Orientation conservation can be inspected in `neighboring_gene_table.tsv`; the JCVI plot encodes strand by arrow direction.\n\n")
        fh.write("## Conclusion\n\n")
        fh.write(f"Case {case}: {conclusion}.\n\n")
        fh.write("## Results paragraph\n\n")
        fh.write(
            f"Using R10821044 as the anchor, we identified {msa_gene} as the corresponding locus in the updated Msa genome. "
            "Microsynteny was then compared between the updated Msa locus, Arabidopsis AT1G27370/SPL10, and the R10821044 locus, "
            "using five protein-coding genes on each side of the focal gene. Reciprocal BLASTP and local-window filtering identified "
            f"{msa_conserved} conserved neighboring gene(s) between the updated Msa and Arabidopsis windows and {r108_conserved} "
            "conserved neighboring gene(s) between the Arabidopsis and R108 windows, excluding the focal SPL candidate itself. "
            "These results provide the basis for assigning the updated Msa locus relative to the Arabidopsis SPL10 neighborhood.\n\n"
        )
        fh.write("## Methods paragraph\n\n")
        fh.write(
            "Microsynteny was assessed using the updated Msa genome annotation, the Medicago truncatula R108 annotation, and the Arabidopsis TAIR10 annotation. "
            "The updated Msa target was first identified by BLASTP using R10821044 as the query, followed by reciprocal BLASTP back to the R108 proteome. "
            "For each target locus, protein-coding genes were sorted by genomic coordinate on the target chromosome, and the five nearest genes upstream and downstream of the target were extracted. "
            "The longest protein isoform per gene was used for homology searches. Reciprocal BLASTP searches were performed between local windows and whole proteomes using an e-value threshold of 1e-5. "
            "Conserved syntenic genes were defined as reciprocal top BLASTP hits that also fell within the compared local windows.\n\n"
        )
        fh.write("## Supplementary figure legend\n\n")
        fh.write(
            "Supplementary Fig. X. R10821044-anchored microsynteny analysis of the SPL10 region in Arabidopsis and Medicago. "
            "The Arabidopsis SPL10 locus (AT1G27370) is shown between the updated Msa locus identified from R10821044 and the R10821044 locus. "
            "Arrows indicate protein-coding genes in genomic order and strand orientation. Gene IDs are displayed above or below each arrow. "
            "Curves connect local-window BLASTP relationships between adjacent tracks only; no direct updated-Msa-to-R108 links are drawn. "
            "The focal SPL candidate comparison is highlighted separately from neighboring-gene homology links.\n"
        )

    print(f"New Msa target from R10821044: {msa_gene}, reciprocal_to_R10821044={msa_is_rbh}")
    print(f"Msa conserved neighboring RBH genes excluding focal SPL: {msa_conserved}")
    print(f"R108 conserved neighboring RBH genes excluding focal SPL: {r108_conserved}")
    print(f"Conclusion case {case}: {conclusion}")


if __name__ == "__main__":
    main()
