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
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def fasta_records(path):
    name, desc, chunks = None, "", []
    with xopen(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks).replace("*", "")
                desc = line[1:]
                name = desc.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, desc, "".join(chunks).replace("*", "")


def write_fasta(records, path):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def sanitize(name):
    return re.sub(r"[^A-Za-z0-9_.|:-]+", "_", name)


def gene_from_header(seqid, desc):
    m = re.search(r"(?:^|\s)gene=([^\s;]+)", desc)
    if m:
        return m.group(1)
    m = re.search(r"(?:^|\s)gene:([^\s;]+)", desc)
    if m:
        return m.group(1).replace("gene:", "")
    if "|" in seqid:
        first = seqid.split("|", 1)[0]
        return first.rsplit(".", 1)[0] if "." in first else first
    return seqid.rsplit(".", 1)[0] if "." in seqid else seqid


def parse_pep(path):
    proteins = {}
    gene_to_proteins = defaultdict(list)
    for seqid, desc, seq in fasta_records(path):
        gid = gene_from_header(seqid, desc)
        proteins[seqid] = {"gene": gid, "desc": desc, "seq": seq}
        gene_to_proteins[gid].append(seqid)
    return proteins, gene_to_proteins


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


def parse_gff(path):
    genes = {}
    with xopen(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attrs = parse_attrs(p[8])
            gid = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
            if not gid:
                continue
            gid = gid.replace("gene:", "")
            genes[gid] = {
                "chrom": p[0],
                "start": p[3],
                "end": p[4],
                "strand": p[6],
            }
    return genes


def run(cmd, log):
    with open(log, "a") as out:
        out.write(" ".join(map(str, cmd)) + "\n")
    subprocess.run(cmd, check=True)


def make_db(makeblastdb, fasta, prefix, log):
    if Path(str(prefix) + ".pin").exists():
        return
    run([makeblastdb, "-in", str(fasta), "-dbtype", "prot", "-out", str(prefix)], log)


def blast(blastp, query, db, out, log, threads=8):
    run([
        blastp, "-query", str(query), "-db", str(db), "-evalue", "1e-5",
        "-max_target_seqs", "200", "-num_threads", str(threads),
        "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs",
        "-out", str(out),
    ], log)


def load_blast_hits(path, subject_proteins, min_qcov):
    best_by_gene = {}
    all_hits = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            rec = {
                "query": p[0],
                "subject_protein": p[1],
                "subject_gene": subject_proteins[p[1]]["gene"],
                "pident": float(p[2]),
                "aln_len": int(p[3]),
                "qlen": int(p[4]),
                "slen": int(p[5]),
                "evalue": float(p[10]),
                "bitscore": float(p[11]),
                "qcovs": float(p[12]),
            }
            if rec["qcovs"] < min_qcov:
                continue
            all_hits.append(rec)
            gid = rec["subject_gene"]
            old = best_by_gene.get(gid)
            if old is None or (rec["bitscore"], -rec["evalue"], rec["qcovs"]) > (old["bitscore"], -old["evalue"], old["qcovs"]):
                best_by_gene[gid] = rec
    return best_by_gene, all_hits


def classify_target(qid):
    if qid.startswith("out_"):
        return "Outgroup"
    if qid.startswith("AtSPL"):
        return "At"
    if qid.startswith("Medtr"):
        return "Medtr"
    return "Target"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-qcov", type=float, default=25.0)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    data = Path(args.data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "blastdb").mkdir(exist_ok=True)
    (out / "blast").mkdir(exist_ok=True)
    log = out / "commands_log.txt"
    blastp = os.environ.get("BLASTP", "blastp")
    makeblastdb = os.environ.get("MAKEBLASTDB", "makeblastdb")

    target_fa = data / "Target.SPL.pep.fa"
    specs = {
        "Msa": {"pep": data / "genome_Msa.pep", "gff": data / "genome_Msa.gff"},
        "R108": {"pep": data / "genome_R108.pep", "gff": data / "genome_R108.gff"},
    }

    target_records = []
    target_query_records = []
    for seqid, desc, seq in fasta_records(target_fa):
        sid = sanitize(seqid)
        target_records.append((f"Ref|{classify_target(sid)}|{sid}", seq))
        if not sid.startswith("out_"):
            target_query_records.append((sid, seq))
    write_fasta(target_query_records, out / "Target.SPL.no_outgroup.query.fa")

    summary_rows = []
    homolog_records = []
    homolog_counts = {}
    for sp, spec in specs.items():
        proteins, gene_to_proteins = parse_pep(spec["pep"])
        genes = parse_gff(spec["gff"])
        db = out / "blastdb" / f"{sp}.pep"
        make_db(makeblastdb, spec["pep"], db, log)
        blast_out = out / "blast" / f"TargetSPL_vs_{sp}.blastp.tsv"
        blast(blastp, out / "Target.SPL.no_outgroup.query.fa", db, blast_out, log, args.threads)
        best_by_gene, all_hits = load_blast_hits(blast_out, proteins, args.min_qcov)
        homolog_counts[sp] = len(best_by_gene)

        with open(out / f"{sp}.all_filtered_blast_hits.tsv", "w") as fh:
            fh.write("Query\tSubject_gene\tSubject_protein\tPident\tAlignment_length\tQuery_length\tSubject_length\tEvalue\tBitscore\tQcovs\n")
            for h in sorted(all_hits, key=lambda x: (x["query"], -x["bitscore"])):
                fh.write("\t".join(map(str, [
                    h["query"], h["subject_gene"], h["subject_protein"], h["pident"], h["aln_len"],
                    h["qlen"], h["slen"], h["evalue"], h["bitscore"], h["qcovs"],
                ])) + "\n")

        sp_records = []
        for gid, h in sorted(best_by_gene.items(), key=lambda x: x[0]):
            prot = h["subject_protein"]
            coord = genes.get(gid, {"chrom": "", "start": "", "end": "", "strand": ""})
            summary_rows.append({
                "Species": sp,
                "Gene_ID": gid,
                "Protein_ID": prot,
                "Chromosome": coord["chrom"],
                "Start": coord["start"],
                "End": coord["end"],
                "Strand": coord["strand"],
                "Best_reference_hit": h["query"],
                "Pident": f"{h['pident']:.3f}",
                "Alignment_length": str(h["aln_len"]),
                "Query_length": str(h["qlen"]),
                "Subject_length": str(h["slen"]),
                "Evalue": f"{h['evalue']:.3g}",
                "Bitscore": f"{h['bitscore']:.1f}",
                "Qcovs": f"{h['qcovs']:.1f}",
            })
            label = f"{sp}|{gid}|{prot}|best={h['query']}"
            sp_records.append((label, proteins[prot]["seq"]))
            homolog_records.append((label, proteins[prot]["seq"]))
        write_fasta(sp_records, out / f"{sp}.TargetSPL_homologs.fa")

    fields = [
        "Species", "Gene_ID", "Protein_ID", "Chromosome", "Start", "End", "Strand",
        "Best_reference_hit", "Pident", "Alignment_length", "Query_length",
        "Subject_length", "Evalue", "Bitscore", "Qcovs",
    ]
    with open(out / "target_spl_homolog_summary.tsv", "w") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in sorted(summary_rows, key=lambda x: (x["Species"], x["Chromosome"], int(x["Start"] or 0), x["Gene_ID"])):
            fh.write("\t".join(r[f] for f in fields) + "\n")

    write_fasta(target_records + homolog_records, out / "Target_Msa_R108_SPL_for_tree.fa")
    with open(out / "README.txt", "w") as fh:
        fh.write("Target SPL homolog search in chr23997-version genome_Msa and genome_R108.\n")
        fh.write(f"Queries: Target.SPL.pep.fa, excluding out_FLC_AT5G10140 for BLAST search and retaining it in the tree as outgroup.\n")
        fh.write(f"Filtering: BLASTP e-value <= 1e-5 by command plus qcovs >= {args.min_qcov}.\n")
        fh.write(f"Msa homolog genes: {homolog_counts['Msa']}; R108 homolog genes: {homolog_counts['R108']}.\n")
        fh.write("Tree FASTA: Target_Msa_R108_SPL_for_tree.fa.\n")

    print(f"Msa homolog genes: {homolog_counts['Msa']}")
    print(f"R108 homolog genes: {homolog_counts['R108']}")
    print(out / "target_spl_homolog_summary.tsv")
    print(out / "Target_Msa_R108_SPL_for_tree.fa")


if __name__ == "__main__":
    main()
