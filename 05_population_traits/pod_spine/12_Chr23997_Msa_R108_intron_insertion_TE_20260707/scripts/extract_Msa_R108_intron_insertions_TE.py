#!/usr/bin/env python3
import csv
import os
import re
from collections import defaultdict

BASE = "path/to/project/N_4.pod_spiny"
OUT = os.path.join(BASE, "12_Chr23997_Msa_R108_intron_insertion_TE_20260707")

MSA_FASTA = os.path.join(BASE, "00_data/1.reference_genome/genome_Msa.fa")
R108_FASTA = os.path.join(BASE, "00_data/1.reference_genome/genome_R108.fa")
M22_FASTA = os.path.join(BASE, "00_data/1.reference_genome/genome_M22.fa")

MSA_EDTA = os.path.join(BASE, "00_data/3.two_ref/genome_Msa_EDTA")
R108_EDTA = os.path.join(BASE, "00_data/3.two_ref/genome_R108_EDTA")
M22_EDTA = "path/to/project/39.TE_type_soloLTR/data/genome_M22_EDTA"

CHROM = "Chr4"
MSA_START = 88980831
MSA_END = 88981052


def ensure_dirs():
    for sub in ("query", "results", "summary", "scripts"):
        os.makedirs(os.path.join(OUT, sub), exist_ok=True)


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


def extract_interval(fasta, seqid, start, end):
    for name, seq in fasta_iter(fasta):
        if name == seqid:
            return seq[start - 1 : end]
    raise RuntimeError(f"sequence {seqid} not found in {fasta}")


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def exact_hits(fasta, query):
    hits = []
    rc = revcomp(query)
    for name, seq in fasta_iter(fasta):
        start = 0
        while True:
            pos = seq.find(query, start)
            if pos < 0:
                break
            hits.append({"seqid": name, "start": pos + 1, "end": pos + len(query), "strand": "+", "sequence": query})
            start = pos + 1
        start = 0
        while True:
            pos = seq.find(rc, start)
            if pos < 0:
                break
            hits.append({"seqid": name, "start": pos + 1, "end": pos + len(query), "strand": "-", "sequence": rc})
            start = pos + 1
    return hits


def find_files(root, patterns):
    out = []
    root = os.path.realpath(root)
    if not os.path.isdir(root):
        return out
    for cur, _, files in os.walk(root):
        depth = os.path.relpath(cur, root).count(os.sep)
        if depth > 3:
            continue
        for file_name in files:
            for pat in patterns:
                if re.fullmatch(pat, file_name):
                    out.append(os.path.join(cur, file_name))
                    break
    return sorted(out)


def pick_file(root, patterns, prefer=None):
    files = find_files(root, patterns)
    if prefer:
        for file_path in files:
            if prefer in os.path.basename(file_path):
                return file_path
    return files[0] if files else ""


def parse_attrs(attrs):
    parsed = {}
    for part in re.split(r";\s*", attrs.strip()):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            key, value = part, ""
        parsed[key.strip()] = value.strip().strip('"')
    return parsed


def annotation_label(attrs):
    parsed = parse_attrs(attrs)
    vals = []
    for key in ("classification", "Classification", "class", "Class", "Name", "ID", "Target", "Motif"):
        if parsed.get(key):
            vals.append(parsed[key])
    return ";".join(vals[:5]) if vals else attrs[:180]


def overlap_gff(gff, genome, seqid, start, end, query_name):
    rows = []
    if not gff or not os.path.exists(gff):
        return rows
    with open(gff, errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != seqid:
                continue
            try:
                s = int(parts[3])
                e = int(parts[4])
            except ValueError:
                continue
            if s <= end and e >= start:
                rows.append(
                    {
                        "genome": genome,
                        "query_name": query_name,
                        "query": f"{seqid}:{start}-{end}",
                        "annotation_file": gff,
                        "seqid": parts[0],
                        "source": parts[1],
                        "feature_type": parts[2],
                        "start": s,
                        "end": e,
                        "overlap_bp": min(e, end) - max(s, start) + 1,
                        "strand": parts[6],
                        "annotation": annotation_label(parts[8]),
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


def scan_te_library(lib_path, lib_label, genome, seq, k):
    if not lib_path or not os.path.exists(lib_path):
        return []
    q = seq.upper().replace("N", "")
    qr = revcomp(q)
    qsets = {
        "+": {q[i : i + k] for i in range(max(0, len(q) - k + 1))},
        "-": {qr[i : i + k] for i in range(max(0, len(qr) - k + 1))},
    }
    denom = max(1, len(q) - k + 1)
    rows = []
    for header, s in fasta_iter(lib_path):
        if len(s) < k:
            continue
        plus = set()
        minus = set()
        for i in range(len(s) - k + 1):
            token = s[i : i + k]
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
                    "subject_length": len(s),
                    "header": header[:220],
                    "file": lib_path,
                }
            )
    rows.sort(key=lambda row: (row["kmer_coverage"], row["matched_unique_kmers"]), reverse=True)
    return rows[:25]


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def main():
    ensure_dirs()

    msa_seq = extract_interval(MSA_FASTA, CHROM, MSA_START, MSA_END)
    r108_hits = exact_hits(R108_FASTA, msa_seq)
    m22_hits = exact_hits(M22_FASTA, msa_seq)

    if r108_hits:
        r108_hit = r108_hits[0]
        r108_seq = r108_hit["sequence"]
    else:
        r108_hit = {"seqid": "NA", "start": "", "end": "", "strand": "", "sequence": ""}
        r108_seq = ""

    records = [
        (f"genome_Msa_Chr23997_intron_insertion|{CHROM}:{MSA_START}-{MSA_END}|len={len(msa_seq)}", msa_seq)
    ]
    if r108_seq:
        records.append(
            (
                f"genome_R108_Chr23997_intron_insertion_exact_Msa_match|{r108_hit['seqid']}:{r108_hit['start']}-{r108_hit['end']}({r108_hit['strand']})|len={len(r108_seq)}",
                r108_seq,
            )
        )

    query_fa = os.path.join(OUT, "query", "Chr23997_Msa_R108_intron_insertions.fa")
    write_fasta(query_fa, records)

    seq_summary = [
        {
            "genome": "genome_Msa",
            "seqid": CHROM,
            "start": MSA_START,
            "end": MSA_END,
            "strand": "+",
            "length": len(msa_seq),
            "sequence": msa_seq,
            "note": "extracted from genome_Msa reference coordinate",
        }
    ]
    if r108_seq:
        seq_summary.append(
            {
                "genome": "genome_R108",
                "seqid": r108_hit["seqid"],
                "start": r108_hit["start"],
                "end": r108_hit["end"],
                "strand": r108_hit["strand"],
                "length": len(r108_seq),
                "sequence": r108_seq,
                "note": f"exact match to genome_Msa insertion; total_R108_exact_hits={len(r108_hits)}",
            }
        )
    else:
        seq_summary.append(
            {
                "genome": "genome_R108",
                "seqid": "NA",
                "start": "",
                "end": "",
                "strand": "",
                "length": 0,
                "sequence": "",
                "note": "no exact match to genome_Msa insertion found",
            }
        )
    seq_summary.append(
        {
            "genome": "genome_M22",
            "seqid": "NA",
            "start": "",
            "end": "",
            "strand": "",
            "length": 0,
            "sequence": "",
            "note": f"exact matches to Msa insertion in M22={len(m22_hits)}",
        }
    )
    seq_tsv = os.path.join(OUT, "results", "Chr23997_Msa_R108_intron_insertion_sequences.tsv")
    write_tsv(seq_tsv, seq_summary, ["genome", "seqid", "start", "end", "strand", "length", "sequence", "note"])

    msa_teanno = pick_file(MSA_EDTA, [r".*TEanno\.gff3"], prefer="TEanno")
    msa_rm = pick_file(MSA_EDTA, [r".*\.out\.gff"], prefer=".out")
    r108_teanno = pick_file(R108_EDTA, [r".*TEanno\.gff3"], prefer="TEanno")
    r108_rm = pick_file(R108_EDTA, [r".*\.out\.gff"], prefer=".out")

    overlap_rows = []
    overlap_rows.extend(overlap_gff(msa_teanno, "genome_Msa", CHROM, MSA_START, MSA_END, "Msa_insertion_exact"))
    overlap_rows.extend(overlap_gff(msa_rm, "genome_Msa", CHROM, MSA_START, MSA_END, "Msa_insertion_exact"))
    if r108_seq:
        overlap_rows.extend(
            overlap_gff(
                r108_teanno,
                "genome_R108",
                r108_hit["seqid"],
                int(r108_hit["start"]),
                int(r108_hit["end"]),
                "R108_insertion_exact_match",
            )
        )
        overlap_rows.extend(
            overlap_gff(
                r108_rm,
                "genome_R108",
                r108_hit["seqid"],
                int(r108_hit["start"]),
                int(r108_hit["end"]),
                "R108_insertion_exact_match",
            )
        )
    overlap_tsv = os.path.join(OUT, "results", "Chr23997_Msa_R108_intron_insertion_EDTA_RepeatMasker_overlap.tsv")
    write_tsv(
        overlap_tsv,
        overlap_rows,
        [
            "genome",
            "query_name",
            "query",
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

    msa_lib = pick_file(MSA_EDTA, [r".*TElib\.fa"], prefer="TElib")
    r108_lib = pick_file(R108_EDTA, [r".*TElib\.fa"], prefer="TElib")
    m22_lib = pick_file(M22_EDTA, [r".*TElib\.fa"], prefer="TElib")
    lib_rows = []
    for genome, seq in (("genome_Msa", msa_seq), ("genome_R108", r108_seq)):
        if not seq:
            continue
        for k in (17, 21):
            lib_rows.extend(scan_te_library(msa_lib, "Msa_EDTA_TElib", genome, seq, k))
            lib_rows.extend(scan_te_library(r108_lib, "R108_EDTA_TElib", genome, seq, k))
            lib_rows.extend(scan_te_library(m22_lib, "M22_EDTA_TElib", genome, seq, k))
    lib_tsv = os.path.join(OUT, "results", "Chr23997_Msa_R108_intron_insertion_TElib_kmer_hits.tsv")
    write_tsv(
        lib_tsv,
        lib_rows,
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

    best_by_genome = defaultdict(list)
    for row in sorted(lib_rows, key=lambda r: (int(r["query_k"]), float(r["kmer_coverage"])), reverse=True):
        if len(best_by_genome[row["genome"]]) < 5:
            best_by_genome[row["genome"]].append(row)

    summary_md = os.path.join(OUT, "summary", "Chr23997_Msa_R108_intron_insertion_TE_annotation.md")
    with open(summary_md, "w") as handle:
        handle.write("# Chr23997 Msa/R108 intron insertion TE annotation\n\n")
        handle.write(f"Msa insertion coordinate: `{CHROM}:{MSA_START}-{MSA_END}`; length `{len(msa_seq)}` bp.\n\n")
        if r108_seq:
            handle.write(
                f"R108 insertion coordinate: `{r108_hit['seqid']}:{r108_hit['start']}-{r108_hit['end']}({r108_hit['strand']})`; length `{len(r108_seq)}` bp; exact matches in R108: `{len(r108_hits)}`.\n\n"
            )
            handle.write(f"Msa and R108 insertion sequences identical: `{msa_seq == r108_seq}`.\n\n")
        else:
            handle.write("No exact R108 match to the Msa insertion sequence was found.\n\n")
        handle.write(f"Exact matches in M22: `{len(m22_hits)}`.\n\n")
        handle.write("## EDTA / RepeatMasker overlap\n\n")
        if overlap_rows:
            for row in overlap_rows:
                handle.write(
                    f"- {row['genome']} {row['query_name']}: {row['source']} {row['feature_type']} {row['seqid']}:{row['start']}-{row['end']}, overlap={row['overlap_bp']} bp, annotation={row['annotation']}\n"
                )
        else:
            handle.write("No EDTA/RepeatMasker overlap was detected for the exact insertion sequences.\n")
        handle.write("\n## Local EDTA TE-library k-mer support\n\n")
        if best_by_genome:
            for genome, rows in best_by_genome.items():
                handle.write(f"### {genome}\n\n")
                for row in rows:
                    handle.write(
                        f"- {row['library']} k={row['query_k']}, coverage={float(row['kmer_coverage']):.3f}, matched={row['matched_unique_kmers']}/{row['query_kmers']}, class={row['class']}, header={row['header']}\n"
                    )
        else:
            handle.write("No local EDTA TE-library k-mer hits were detected.\n")
        handle.write("\n## Output files\n\n")
        for path in (query_fa, seq_tsv, overlap_tsv, lib_tsv):
            handle.write(f"- `{path}`\n")

    print("OUT", OUT)
    print("QUERY_FASTA", query_fa)
    print("SEQUENCE_TSV", seq_tsv)
    print("OVERLAP_TSV", overlap_tsv)
    print("KMER_TSV", lib_tsv)
    print("SUMMARY", summary_md)
    print("MSA_LEN", len(msa_seq))
    print("R108_EXACT_HITS", len(r108_hits))
    if r108_seq:
        print("R108_COORD", f"{r108_hit['seqid']}:{r108_hit['start']}-{r108_hit['end']}({r108_hit['strand']})")
        print("MSA_R108_IDENTICAL", msa_seq == r108_seq)
    print("M22_EXACT_HITS", len(m22_hits))
    print("OVERLAPS", len(overlap_rows))
    for row in overlap_rows[:12]:
        print("OV", row["genome"], row["source"], row["feature_type"], row["query"], row["start"], row["end"], row["overlap_bp"], row["annotation"])
    print("BEST_KMER")
    for genome, rows in best_by_genome.items():
        for row in rows[:5]:
            print("HIT", genome, row["library"], "k", row["query_k"], "cov", f"{float(row['kmer_coverage']):.3f}", "class", row["class"], "header", row["header"][:100])


if __name__ == "__main__":
    main()
