#!/usr/bin/env python3
import csv
import os
import re

BASE = "path/to/project/N_4.pod_spiny"
TEBASE = "path/to/project/39.TE_type_soloLTR/data"
OUT = os.path.join(BASE, "11_Chr23997_intron2_TE_check_20260707")

CHROM = "Chr4"
DEL_START = 88980831
DEL_END = 88981052
INTRON_START = 88980584
INTRON_END = 88981373
GENE_START = 88979818
GENE_END = 88982331

FALLBACK_SEQ = (
    "ACTGTTTTCTTGTGCCTTTAAGTAAATCCACTATACAGCCCACAACTCTAATAGCAATCCCTTCTAGCAACTCTAATATGACCATTACTCTGCATGGTTATTATGAAATTGATACTGCTAGAACTCTTGAACATGCTTTATCTTATTCGTGTCACAAAGGATTAAGCATCTTTATCTTTAGGATGTATGACAATTGAAAATCCAAAAAGATCCATGAAATAA"
)


def ensure_dirs():
    for sub in ("query", "results", "summary", "scripts"):
        os.makedirs(os.path.join(OUT, sub), exist_ok=True)


def load_query_sequence():
    records_tsv = os.path.join(
        BASE,
        "05_direction_free_sv_cosegregation_20260701/summary/Chr23997_caller_records_20260701.tsv",
    )
    if os.path.exists(records_tsv):
        with open(records_tsv, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                event_id = row.get("id", "")
                same_interval = (
                    row.get("chrom") == CHROM
                    and row.get("pos") == str(DEL_START)
                    and row.get("end") == str(DEL_END)
                    and row.get("svtype") == "DEL"
                )
                if "0_0_pbsv.DEL.130205" in event_id or same_interval:
                    seq = (row.get("ref") or "").strip().upper()
                    if seq and set(seq) <= set("ACGTN") and len(seq) >= 50:
                        return seq, records_tsv
    return FALLBACK_SEQ, "fallback_from_prior_pbsv_record"


def write_fasta(seq):
    old_path = os.path.join(OUT, "query", "Chr23997_M22_deleted_intron2_221bp.Msa_ref.fa")
    if os.path.exists(old_path):
        os.remove(old_path)
    path = os.path.join(OUT, "query", "Chr23997_M22_deleted_intron2_222bp.Msa_ref.fa")
    with open(path, "w") as handle:
        handle.write(
            f">Chr23997_M22_deleted_intron2_221bp|{CHROM}:{DEL_START}-{DEL_END}|present_in_R108_Msa_ref_absent_in_M22\n"
        )
        for i in range(0, len(seq), 80):
            handle.write(seq[i : i + 80] + "\n")
    return path


def real_path(*parts):
    return os.path.realpath(os.path.join(*parts))


PATHS = {
    "Msa_EDTA_TEanno": real_path(
        TEBASE,
        "genome_Msa_T2T_EDTA/genome_Msa.T2T.ctg.final.fa.mod.EDTA.TEanno.gff3",
    ),
    "Msa_RepeatMasker_gff": real_path(
        TEBASE, "genome_Msa_T2T_EDTA/RepeatMasker_softmask/genome_Msa.T2T.ctg.final.fa.out.gff"
    ),
    "R108_EDTA_TEanno": real_path(
        TEBASE,
        "genome_R108_EDTA/genome_R108.T2T.ctg.final.fa.mod.EDTA.TEanno.gff3",
    ),
    "R108_RepeatMasker_gff": real_path(
        TEBASE, "genome_R108_EDTA/RepeatMasker_softmask/genome_R108.T2T.ctg.final.fa.out.gff"
    ),
    "M22_EDTA_TEanno": real_path(
        TEBASE, "genome_M22_EDTA/X_genome_M22.genome.fa.mod.EDTA.TEanno.gff3"
    ),
    "M22_TElib": real_path(TEBASE, "genome_M22_EDTA/X_genome_M22.genome.fa.mod.EDTA.TElib.fa"),
    "Msa_TElib": real_path(
        TEBASE, "genome_Msa_T2T_EDTA/genome_Msa.T2T.ctg.final.fa.mod.EDTA.TElib.fa"
    ),
    "R108_TElib": real_path(
        TEBASE, "genome_R108_EDTA/genome_R108.T2T.ctg.final.fa.mod.EDTA.TElib.fa"
    ),
}


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


def feature_label(attrs):
    parsed = parse_attrs(attrs)
    labels = []
    for key in ("Classification", "class", "Class", "Name", "ID", "Target", "Method", "Motif"):
        if parsed.get(key):
            labels.append(parsed[key])
    return ";".join(labels[:4]) if labels else attrs[:180]


def overlap_rows(gff_path, qname, qchrom, qstart, qend, dataset):
    rows = []
    if not os.path.exists(gff_path):
        return rows
    with open(gff_path, errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != qchrom:
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            if start <= qend and end >= qstart:
                overlap_bp = min(end, qend) - max(start, qstart) + 1
                rows.append(
                    {
                        "dataset": dataset,
                        "query_name": qname,
                        "query": f"{qchrom}:{qstart}-{qend}",
                        "seqid": parts[0],
                        "source": parts[1],
                        "type": parts[2],
                        "start": start,
                        "end": end,
                        "overlap_bp": max(0, overlap_bp),
                        "strand": parts[6],
                        "annotation": feature_label(parts[8]),
                        "attributes": parts[8],
                        "file": gff_path,
                    }
                )
    return rows


def collect_overlaps():
    queries = [
        ("candidate_DEL_221bp", CHROM, DEL_START, DEL_END),
        ("candidate_DEL_plus_500bp", CHROM, DEL_START - 500, DEL_END + 500),
        ("intron2_interval", CHROM, INTRON_START, INTRON_END),
        ("gene_body", CHROM, GENE_START, GENE_END),
    ]
    datasets = [
        ("Msa_EDTA_TEanno", PATHS["Msa_EDTA_TEanno"]),
        ("Msa_RepeatMasker_gff", PATHS["Msa_RepeatMasker_gff"]),
        ("R108_EDTA_TEanno_same_coord_check", PATHS["R108_EDTA_TEanno"]),
        ("R108_RepeatMasker_same_coord_check", PATHS["R108_RepeatMasker_gff"]),
    ]
    rows = []
    for qname, qchrom, qstart, qend in queries:
        for dataset, path in datasets:
            rows.extend(overlap_rows(path, qname, qchrom, qstart, qend, dataset))
    return rows


def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


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
                name = line[1:]
                chunks = []
            else:
                chunks.append(re.sub(r"[^ACGTNacgtn]", "", line))
        if name is not None:
            yield name, "".join(chunks).upper()


def class_from_header(header):
    if "#" in header:
        return header.split("#", 1)[1].split()[0]
    for token in re.split(r"\s+", header):
        if "/" in token and any(word in token for word in ("LTR", "DNA", "LINE", "SINE", "Helitron", "TIR")):
            return token
    return "NA"


def kmer_scan(lib_path, lib_label, query_seq, k):
    if not os.path.exists(lib_path):
        return []
    q_forward = query_seq.upper().replace("N", "")
    q_reverse = revcomp(query_seq).replace("N", "")
    qsets = {
        "+": {q_forward[i : i + k] for i in range(0, max(0, len(q_forward) - k + 1))},
        "-": {q_reverse[i : i + k] for i in range(0, max(0, len(q_reverse) - k + 1))},
    }
    n_query_kmers = max(1, len(q_forward) - k + 1)
    hits = []
    for header, seq in fasta_iter(lib_path):
        if len(seq) < k:
            continue
        plus = set()
        minus = set()
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i : i + k]
            if kmer in qsets["+"]:
                plus.add(kmer)
            if kmer in qsets["-"]:
                minus.add(kmer)
        if plus or minus:
            if len(plus) >= len(minus):
                strand, matched = "+", len(plus)
            else:
                strand, matched = "-", len(minus)
            hits.append(
                {
                    "library": lib_label,
                    "class": class_from_header(header),
                    "strand": strand,
                    "query_k": k,
                    "query_kmers": n_query_kmers,
                    "matched_unique_kmers": matched,
                    "kmer_coverage": matched / n_query_kmers,
                    "subject_length": len(seq),
                    "header": header[:220],
                    "file": lib_path,
                }
            )
    hits.sort(key=lambda row: (row["kmer_coverage"], row["matched_unique_kmers"]), reverse=True)
    return hits[:25]


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main():
    ensure_dirs()
    query_seq, seq_source = load_query_sequence()
    query_fasta = write_fasta(query_seq)

    overlaps = collect_overlaps()
    overlap_tsv = os.path.join(OUT, "results", "Chr23997_M22_deleted_segment_TE_overlap.tsv")
    write_tsv(
        overlap_tsv,
        overlaps,
        [
            "dataset",
            "query_name",
            "query",
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "overlap_bp",
            "strand",
            "annotation",
            "attributes",
            "file",
        ],
    )

    kmer_hits = []
    for k in (17, 21):
        for label in ("Msa_TElib", "R108_TElib", "M22_TElib"):
            kmer_hits.extend(kmer_scan(PATHS[label], label, query_seq, k))
    kmer_tsv = os.path.join(OUT, "results", "Chr23997_M22_deleted_segment_TElib_kmer_hits.tsv")
    write_tsv(
        kmer_tsv,
        kmer_hits,
        [
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

    exact_msa = [
        row
        for row in overlaps
        if row["query_name"] == "candidate_DEL_221bp" and row["dataset"].startswith("Msa_")
    ]
    nearby_msa = [row for row in overlaps if row["dataset"].startswith("Msa_")]
    best_hits = sorted(
        kmer_hits,
        key=lambda row: (int(row["query_k"]), float(row["kmer_coverage"]), int(row["matched_unique_kmers"])),
        reverse=True,
    )[:20]

    summary_path = os.path.join(OUT, "summary", "Chr23997_M22_intron2_deleted_segment_TE_check.md")
    with open(summary_path, "w") as handle:
        handle.write("# Chr23997 M22 intron deletion TE check\n\n")
        handle.write(
            f"Query segment: `{CHROM}:{DEL_START}-{DEL_END}` ({len(query_seq)} bp), present in the R108/Msa-intact state and absent in M22 by read-depth/genotype evidence.\n\n"
        )
        handle.write(f"Query sequence source: `{seq_source}`\n\n")
        handle.write("## Direct TE annotation overlap\n\n")
        if exact_msa:
            for row in exact_msa[:30]:
                handle.write(
                    f"- {row['dataset']}: {row['seqid']}:{row['start']}-{row['end']}; overlap={row['overlap_bp']} bp; type={row['type']}; annotation={row['annotation']}\n"
                )
        else:
            handle.write("No Msa EDTA/RepeatMasker feature directly overlaps the exact 221 bp interval.\n")
        handle.write("\n## Nearby Msa TE annotation features\n\n")
        if nearby_msa:
            for row in nearby_msa[:40]:
                handle.write(
                    f"- {row['dataset']} {row['query_name']} {row['query']}: {row['seqid']}:{row['start']}-{row['end']}; overlap={row['overlap_bp']} bp; type={row['type']}; annotation={row['annotation']}\n"
                )
        else:
            handle.write("No Msa TE annotation found in the queried gene/nearby windows.\n")
        handle.write("\n## TE-library k-mer similarity\n\n")
        if best_hits:
            for row in best_hits:
                handle.write(
                    f"- {row['library']} k={row['query_k']} coverage={float(row['kmer_coverage']):.3f} matched={row['matched_unique_kmers']}/{row['query_kmers']} class={row['class']} header={row['header']}\n"
                )
        else:
            handle.write("No 17/21-mer support against the Msa/R108/M22 EDTA TE libraries.\n")
        handle.write("\n## Output files\n\n")
        handle.write(f"- `{query_fasta}`\n")
        handle.write(f"- `{overlap_tsv}`\n")
        handle.write(f"- `{kmer_tsv}`\n")

    print("OUT", OUT)
    print("QUERY_FASTA", query_fasta)
    print("QUERY_LEN", len(query_seq))
    print("QUERY_SOURCE", seq_source)
    print("OVERLAP_TSV", overlap_tsv)
    print("KMER_TSV", kmer_tsv)
    print("SUMMARY", summary_path)
    print("EXACT_MSA_OVERLAPS", len(exact_msa))
    for row in exact_msa[:10]:
        print(
            "EXACT",
            row["dataset"],
            row["seqid"],
            row["start"],
            row["end"],
            row["overlap_bp"],
            row["type"],
            row["annotation"],
        )
    print("BEST_KMER_HITS")
    for row in best_hits[:10]:
        print(
            "HIT",
            row["library"],
            "k",
            row["query_k"],
            "cov",
            f"{float(row['kmer_coverage']):.3f}",
            "matched",
            f"{row['matched_unique_kmers']}/{row['query_kmers']}",
            "class",
            row["class"],
            "header",
            row["header"][:120],
        )


if __name__ == "__main__":
    main()
