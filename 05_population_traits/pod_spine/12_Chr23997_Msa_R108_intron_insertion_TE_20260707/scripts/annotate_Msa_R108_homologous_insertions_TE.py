#!/usr/bin/env python3
import csv
import os
import re
from collections import defaultdict

BASE = "path/to/project/N_4.pod_spiny"
OUT = os.path.join(BASE, "12_Chr23997_Msa_R108_intron_insertion_TE_20260707")

MSA_EDTA = os.path.join(BASE, "00_data/3.two_ref/genome_Msa_EDTA")
R108_EDTA = os.path.join(BASE, "00_data/3.two_ref/genome_R108_EDTA")
M22_EDTA = "path/to/project/39.TE_type_soloLTR/data/genome_M22_EDTA"

QUERIES = [
    {"genome": "genome_Msa", "name": "Msa_intron_insertion", "seqid": "Chr4", "start": 88980831, "end": 88981052},
    {"genome": "genome_R108", "name": "R108_homologous_intron_insertion", "seqid": "Chr4", "start": 45421964, "end": 45422252},
]


def fasta_iter(path):
    name = None
    chunks = []
    with open(path, errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(re.sub(r"[^ACGTNacgtn]", "", line))
        if name is not None:
            yield name, "".join(chunks).upper()


def load_query_fasta():
    path = os.path.join(OUT, "query", "Chr23997_Msa_R108_homologous_intron_insertions.fa")
    seqs = {}
    for name, seq in fasta_iter(path):
        if name.startswith("genome_Msa"):
            seqs["genome_Msa"] = seq
        elif name.startswith("genome_R108"):
            seqs["genome_R108"] = seq
    return path, seqs


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def find_files(root, regex):
    root = os.path.realpath(root)
    out = []
    if not os.path.isdir(root):
        return out
    for cur, _, files in os.walk(root):
        depth = os.path.relpath(cur, root).count(os.sep)
        if depth > 3:
            continue
        for f in files:
            if re.fullmatch(regex, f):
                out.append(os.path.join(cur, f))
    return sorted(out)


def pick(root, regex):
    files = find_files(root, regex)
    return files[0] if files else ""


def parse_attrs(attrs):
    parsed = {}
    for part in re.split(r";\s*", attrs.strip()):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
        else:
            k, v = part, ""
        parsed[k.strip()] = v.strip().strip('"')
    return parsed


def label(attrs):
    parsed = parse_attrs(attrs)
    vals = []
    for key in ("classification", "Classification", "class", "Class", "Name", "ID", "Target", "Motif"):
        if parsed.get(key):
            vals.append(parsed[key])
    return ";".join(vals[:5]) if vals else attrs[:160]


def overlap(gff, query, dataset):
    rows = []
    if not gff or not os.path.exists(gff):
        return rows
    with open(gff, errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != query["seqid"]:
                continue
            try:
                s, e = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            if s <= query["end"] and e >= query["start"]:
                rows.append(
                    {
                        "genome": query["genome"],
                        "query_name": query["name"],
                        "query": f"{query['seqid']}:{query['start']}-{query['end']}",
                        "dataset": dataset,
                        "annotation_file": gff,
                        "seqid": parts[0],
                        "source": parts[1],
                        "feature_type": parts[2],
                        "start": s,
                        "end": e,
                        "overlap_bp": min(e, query["end"]) - max(s, query["start"]) + 1,
                        "strand": parts[6],
                        "annotation": label(parts[8]),
                        "attributes": parts[8],
                    }
                )
    return rows


def class_from_header(header):
    if "#" in header:
        return header.split("#", 1)[1].split()[0]
    for token in re.split(r"\s+", header):
        if "/" in token:
            return token
    return "NA"


def scan_library(lib, lib_label, genome, seq, k):
    if not lib or not os.path.exists(lib):
        return []
    q = seq.upper().replace("N", "")
    qr = revcomp(q)
    qsets = {
        "+": {q[i : i + k] for i in range(max(0, len(q) - k + 1))},
        "-": {qr[i : i + k] for i in range(max(0, len(qr) - k + 1))},
    }
    denom = max(1, len(q) - k + 1)
    rows = []
    for header, target in fasta_iter(lib):
        if len(target) < k:
            continue
        plus, minus = set(), set()
        for i in range(len(target) - k + 1):
            token = target[i : i + k]
            if token in qsets["+"]:
                plus.add(token)
            if token in qsets["-"]:
                minus.add(token)
        if plus or minus:
            if len(plus) >= len(minus):
                strand, n = "+", len(plus)
            else:
                strand, n = "-", len(minus)
            rows.append(
                {
                    "genome": genome,
                    "library": lib_label,
                    "class": class_from_header(header),
                    "strand": strand,
                    "query_k": k,
                    "query_kmers": denom,
                    "matched_unique_kmers": n,
                    "kmer_coverage": n / denom,
                    "subject_length": len(target),
                    "header": header[:220],
                    "file": lib,
                }
            )
    rows.sort(key=lambda r: (r["kmer_coverage"], r["matched_unique_kmers"]), reverse=True)
    return rows[:25]


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main():
    os.makedirs(os.path.join(OUT, "results"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "summary"), exist_ok=True)
    query_fa, seqs = load_query_fasta()

    gffs = {
        "genome_Msa_EDTA_TEanno": pick(MSA_EDTA, r".*TEanno\.gff3"),
        "genome_Msa_RepeatMasker_gff": pick(MSA_EDTA, r".*\.out\.gff"),
        "genome_R108_EDTA_TEanno": pick(R108_EDTA, r".*TEanno\.gff3"),
        "genome_R108_RepeatMasker_gff": pick(R108_EDTA, r".*\.out\.gff"),
    }
    libs = {
        "Msa_EDTA_TElib": pick(MSA_EDTA, r".*TElib\.fa"),
        "R108_EDTA_TElib": pick(R108_EDTA, r".*TElib\.fa"),
        "M22_EDTA_TElib": pick(M22_EDTA, r".*TElib\.fa"),
    }

    overlap_rows = []
    for q in QUERIES:
        if q["genome"] == "genome_Msa":
            overlap_rows += overlap(gffs["genome_Msa_EDTA_TEanno"], q, "EDTA_TEanno")
            overlap_rows += overlap(gffs["genome_Msa_RepeatMasker_gff"], q, "RepeatMasker_gff")
        else:
            overlap_rows += overlap(gffs["genome_R108_EDTA_TEanno"], q, "EDTA_TEanno")
            overlap_rows += overlap(gffs["genome_R108_RepeatMasker_gff"], q, "RepeatMasker_gff")

    overlap_tsv = os.path.join(OUT, "results", "Chr23997_Msa_R108_homologous_insertion_EDTA_RepeatMasker_overlap.tsv")
    write_tsv(
        overlap_tsv,
        overlap_rows,
        [
            "genome",
            "query_name",
            "query",
            "dataset",
            "annotation_file",
            "seqid",
            "source",
            "feature_type",
            "start",
            "end",
            "overlap_bp",
            "strand",
            "annotation",
            "attributes",
        ],
    )

    kmer_rows = []
    for genome, seq in seqs.items():
        for k in (17, 21):
            for lib_label, lib_path in libs.items():
                kmer_rows += scan_library(lib_path, lib_label, genome, seq, k)
    kmer_tsv = os.path.join(OUT, "results", "Chr23997_Msa_R108_homologous_insertion_TElib_kmer_hits.tsv")
    write_tsv(
        kmer_tsv,
        kmer_rows,
        [
            "genome",
            "library",
            "class",
            "strand",
            "query_k",
            "query_kmers",
            "matched_unique_kmers",
            "kmer_coverage",
            "subject_length",
            "header",
            "file",
        ],
    )

    best = defaultdict(list)
    for row in sorted(kmer_rows, key=lambda r: (int(r["query_k"]), float(r["kmer_coverage"]), int(r["matched_unique_kmers"])), reverse=True):
        if len(best[row["genome"]]) < 8:
            best[row["genome"]].append(row)

    summary = os.path.join(OUT, "summary", "Chr23997_Msa_R108_homologous_insertion_TE_annotation.md")
    with open(summary, "w") as handle:
        handle.write("# Chr23997 Msa/R108 homologous intron insertion TE annotation\n\n")
        handle.write(f"Query FASTA: `{query_fa}`\n\n")
        handle.write("## Query intervals\n\n")
        for q in QUERIES:
            seq = seqs.get(q["genome"], "")
            handle.write(f"- {q['genome']}: `{q['seqid']}:{q['start']}-{q['end']}`, length={len(seq)} bp\n")
        handle.write("\n## EDTA / RepeatMasker overlaps\n\n")
        if overlap_rows:
            for row in overlap_rows:
                handle.write(
                    f"- {row['genome']} {row['dataset']}: {row['feature_type']} {row['seqid']}:{row['start']}-{row['end']}, overlap={row['overlap_bp']} bp, annotation={row['annotation']}\n"
                )
        else:
            handle.write("No exact interval overlap with EDTA/RepeatMasker TE features.\n")
        handle.write("\n## Local EDTA TE-library k-mer hits\n\n")
        if best:
            for genome, rows in best.items():
                handle.write(f"### {genome}\n\n")
                for row in rows:
                    handle.write(
                        f"- {row['library']} k={row['query_k']}, coverage={float(row['kmer_coverage']):.3f}, matched={row['matched_unique_kmers']}/{row['query_kmers']}, class={row['class']}, header={row['header']}\n"
                    )
        else:
            handle.write("No k-mer support against local EDTA TE libraries.\n")
        handle.write("\n## Output files\n\n")
        handle.write(f"- `{overlap_tsv}`\n")
        handle.write(f"- `{kmer_tsv}`\n")

    print("QUERY_FASTA", query_fa)
    print("OVERLAP_TSV", overlap_tsv)
    print("KMER_TSV", kmer_tsv)
    print("SUMMARY", summary)
    print("OVERLAPS", len(overlap_rows))
    for row in overlap_rows:
        print("OV", row["genome"], row["dataset"], row["feature_type"], row["query"], row["start"], row["end"], row["overlap_bp"], row["annotation"])
    print("BEST")
    for genome, rows in best.items():
        for row in rows:
            print("HIT", genome, row["library"], "k", row["query_k"], "cov", f"{float(row['kmer_coverage']):.3f}", "class", row["class"], "header", row["header"][:120])


if __name__ == "__main__":
    main()
