#!/usr/bin/env python3
import csv
import os
import re
from collections import Counter, defaultdict

BASE = "path/to/project/N_4.pod_spiny"
REFDIR = os.path.join(BASE, "00_data/1.reference_genome")
OUT = os.path.join(BASE, "15_Chr23997_Msa_reference_similarity_RBH_20260708")
ORTHO_SUMMARY = os.path.join(
    BASE,
    "14_Chr23997_ortholog_sequences_20260708_with_A17",
    "summary",
    "Chr23997_all_best_hit_candidates.summary.tsv",
)

CHROM = "Chr4"
REF_GENOME = "genome_Msa"
REF_FASTA = os.path.join(REFDIR, "genome_Msa.fa")
REF_START = 88979818
REF_END = 88982331
REF_LABEL = "Chr23997"
REF_STRAND = "-"

SEED_K = 31
SEED_STEP = 10
BIN_SIZE = 100
PAD = 1200
FLANK = 500
SMOOTH = 60
STEP = 10

EXONS = [
    ("exon4", 88979818, 88980149),
    ("exon3", 88980368, 88980583),
    ("exon2", 88981374, 88981510),
    ("exon1", 88981709, 88982331),
]

EVENTS = [
    ("second_intron", 88980584, 88981373, "second intron"),
    ("core_INS_DEL", 88980831, 88981052, "222-bp INS/DEL"),
]

SAMPLES = [
    ("genome_R108", "genome_R108", "spiny", "strict_RBH", "main spiny"),
    ("genome_Mpo", "genome_Mpo", "spiny", "strict_RBH", "main spiny"),
    ("genome_474", "genome_474a", "spiny", "strict_RBH", "main spiny; genome_474a representative"),
    ("genome_410", "genome_410", "spiny", "strict_RBH", "main spiny"),
    ("genome_457", "genome_457", "spiny", "strict_RBH", "spiny but low assembly"),
    ("genome_436", "genome_436", "spiny", "strict_RBH", "weak/ambiguous spiny"),
    ("genome_395", "genome_395", "spineless", "strict_RBH", "main spineless"),
    ("genome_454", "genome_454", "spineless", "strict_RBH", "pubescent not true spiny"),
    ("genome_461", "genome_461", "spineless", "strict_RBH", "main spineless"),
    ("genome_468", "genome_468", "spineless", "strict_RBH", "main spineless"),
    ("genome_472", "genome_472", "spineless", "strict_RBH", "main spineless"),
    ("genome_M22", "genome_M22", "spineless", "strict_RBH", "main spineless"),
    ("genome_Mar", "genome_Mar", "spineless", "strict_RBH", "main spineless"),
    ("genome_Mru", "genome_Mru", "spineless", "strict_RBH", "main spineless; opposite strand"),
    ("genome_ZM4", "genome_ZM4", "spineless", "strict_RBH", "main spineless; contig naming may differ"),
    ("genome_482", "genome_482", "candidate", "best_hit_non_RBH", "best-hit candidate only; reciprocal best hit is not Chr23997"),
    ("genome_M46", "genome_M46", "candidate", "best_hit_non_RBH", "best-hit candidate only; reciprocal best hit is not Chr23997"),
]


def clean(line):
    return re.sub(r"[^ACGTNacgtn]", "", line).upper()


def read_seqid(fasta, seqid):
    found = False
    chunks = []
    with open(fasta, errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if found and name != seqid:
                    break
                found = name == seqid
                continue
            if found:
                chunks.append(clean(line.strip()))
    if not chunks:
        raise RuntimeError(f"{seqid} not found in {fasta}")
    return "".join(chunks)


def subseq(seq, start, end):
    return seq[start - 1 : end]


def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def read_ortholog_summary(path):
    rows = {}
    with open(path, newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["genome_id"]] = row
    return rows


def add_seeds(seed_map, ref_seq):
    for i in range(0, len(ref_seq) - SEED_K + 1, SEED_STEP):
        token = ref_seq[i : i + SEED_K]
        if "N" not in token:
            seed_map[token].append(i)


def scan_chr(seq, seed_map):
    seed_set = set(seed_map)
    bins = Counter()
    details = defaultdict(list)
    for pos in range(0, len(seq) - SEED_K + 1):
        token = seq[pos : pos + SEED_K]
        if token not in seed_set:
            continue
        for qpos in seed_map[token]:
            predicted = pos - qpos + 1
            if predicted < 1:
                continue
            b = round(predicted / BIN_SIZE) * BIN_SIZE
            bins[b] += 1
            if len(details[b]) < 300:
                details[b].append(predicted)
    if not bins:
        return None, 0
    best_bin, hits = bins.most_common(1)[0]
    starts = sorted(details[best_bin])
    return starts[len(starts) // 2], hits


def smith_waterman(a, b):
    match, mismatch, gap = 2, -2, -3
    n, m = len(a), len(b)
    prev = [0] * (m + 1)
    ptr = [bytearray(m + 1) for _ in range(n + 1)]
    best = (0, 0, 0)
    for i in range(1, n + 1):
        cur = [0] * (m + 1)
        ai = a[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match if ai == b[j - 1] else mismatch)
            up = prev[j] + gap
            left = cur[j - 1] + gap
            val = max(0, diag, up, left)
            cur[j] = val
            if val == 0:
                ptr[i][j] = 0
            elif val == diag:
                ptr[i][j] = 1
            elif val == up:
                ptr[i][j] = 2
            else:
                ptr[i][j] = 3
            if val > best[0]:
                best = (val, i, j)
        prev = cur
    score, i, j = best
    end_a, end_b = i, j
    matches = aligned = 0
    while i > 0 and j > 0 and ptr[i][j] != 0:
        d = ptr[i][j]
        if d == 1:
            if a[i - 1] == b[j - 1]:
                matches += 1
            aligned += 1
            i -= 1
            j -= 1
        elif d == 2:
            aligned += 1
            i -= 1
        else:
            aligned += 1
            j -= 1
    return {
        "score": score,
        "query_start": i + 1,
        "query_end": end_a,
        "target_start": j + 1,
        "target_end": end_b,
        "identity": matches / aligned if aligned else 0,
        "aligned": aligned,
    }


def global_identity_marks(ref, target):
    match, mismatch, gap = 2, -1, -2
    n, m = len(ref), len(target)
    prev = [j * gap for j in range(m + 1)]
    ptr = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        ptr[0][j] = 3
    for i in range(1, n + 1):
        cur = [0] * (m + 1)
        cur[0] = i * gap
        ptr[i][0] = 2
        ai = ref[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match if ai == target[j - 1] else mismatch)
            up = prev[j] + gap
            left = cur[j - 1] + gap
            if diag >= up and diag >= left:
                cur[j] = diag
                ptr[i][j] = 1
            elif up >= left:
                cur[j] = up
                ptr[i][j] = 2
            else:
                cur[j] = left
                ptr[i][j] = 3
        prev = cur
    i, j = n, m
    marks = []
    aligned_ref = aligned = matches = 0
    while i > 0 or j > 0:
        d = ptr[i][j] if i >= 0 and j >= 0 else 0
        if i > 0 and j > 0 and d == 1:
            val = 1 if ref[i - 1] == target[j - 1] else 0
            marks.append(val)
            aligned_ref += 1
            aligned += 1
            matches += val
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or d == 2):
            marks.append(0)
            aligned_ref += 1
            aligned += 1
            i -= 1
        else:
            aligned += 1
            j -= 1
    marks.reverse()
    return marks, matches, aligned_ref, aligned


def smooth_marks(marks, window=SMOOTH, step=STEP):
    out = []
    n = len(marks)
    for start in range(0, n, step):
        end = min(n, start + window)
        if end <= start:
            continue
        vals = marks[start:end]
        out.append((start + 1 + (end - start) / 2, 100 * sum(vals) / len(vals)))
    return out


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main():
    os.makedirs(os.path.join(OUT, "scripts"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "summary"), exist_ok=True)

    msa_chr = read_seqid(REF_FASTA, CHROM)
    ref_seq = subseq(msa_chr, REF_START, REF_END)
    orthologs = read_ortholog_summary(ORTHO_SUMMARY)

    profile_rows = []
    region_rows = []
    sample_rows = []

    for order, (sample, ortholog_id, group, confidence, note) in enumerate(SAMPLES, 1):
        row = orthologs.get(ortholog_id)
        if row is None:
            sample_rows.append(
                {
                    "sample": sample,
                    "ortholog_id": ortholog_id,
                    "fasta": "",
                    "group": group,
                    "confidence": confidence,
                    "plot_order": order,
                    "note": f"{note}; ortholog row missing",
                }
            )
            region_rows.append({"sample": sample, "group": group, "confidence": confidence, "status": "missing_ortholog_row", "note": note})
            continue
        fasta = row["genome_fasta"]
        fasta_name = os.path.basename(fasta)
        sample_rows.append(
            {
                "sample": sample,
                "ortholog_id": ortholog_id,
                "fasta": fasta_name,
                "group": group,
                "confidence": confidence,
                "plot_order": order,
                "note": note,
            }
        )
        if not os.path.exists(fasta):
            region_rows.append({"sample": sample, "group": group, "confidence": confidence, "status": "missing_fasta", "fasta": fasta})
            continue
        chrom = row["chrom"]
        start = int(row["gene_start"])
        end = int(row["gene_end"])
        strand = row["strand"]
        print("LOAD", sample, ortholog_id, chrom, f"{start}-{end}", strand, confidence, flush=True)
        try:
            chrom_seq = read_seqid(fasta, chrom)
        except RuntimeError as exc:
            region_rows.append(
                {
                    "sample": sample,
                    "group": group,
                    "confidence": confidence,
                    "status": "missing_chrom",
                    "fasta": fasta,
                    "target_coord": f"{chrom}:{start}-{end}",
                    "strand": strand,
                    "strict_ortholog": row.get("strict_ortholog", ""),
                    "note": f"{note}; {exc}",
                }
            )
            print("SKIP", sample, "missing", chrom, flush=True)
            continue
        target_seq = subseq(chrom_seq, start, end)
        orientation = "direct"
        if strand != REF_STRAND:
            target_seq = revcomp(target_seq)
            orientation = "reverse_complemented_to_Msa"
        marks, matches, ref_aligned, aligned = global_identity_marks(ref_seq, target_seq)
        ident = matches / ref_aligned if ref_aligned else 0
        for x, sim in smooth_marks(marks):
            profile_rows.append(
                {
                    "sample": sample,
                    "group": group,
                    "confidence": confidence,
                    "x": f"{x:.1f}",
                    "similarity": f"{sim:.3f}",
                }
            )
        region_rows.append(
            {
                "sample": sample,
                "group": group,
                "confidence": confidence,
                "status": "ok_ortholog_coord",
                "ortholog_id": ortholog_id,
                "transcript": row.get("homolog_transcript_id", ""),
                "gene": row.get("homolog_gene_id", ""),
                "target_coord": f"{chrom}:{start}-{end}",
                "strand": strand,
                "orientation": orientation,
                "target_len": len(target_seq),
                "ref_global_identity": f"{ident:.4f}",
                "aligned_columns": aligned,
                "blast_pident": row.get("blast_pident", ""),
                "blast_qcovs": row.get("blast_qcovs", ""),
                "reciprocal_best_msa_hit": row.get("reciprocal_best_msa_hit", ""),
                "strict_ortholog": row.get("strict_ortholog", ""),
                "fasta": fasta,
                "note": note,
            }
        )

    write_tsv(os.path.join(OUT, "data", "samples.tsv"), sample_rows, ["sample", "ortholog_id", "fasta", "group", "confidence", "plot_order", "note"])
    write_tsv(os.path.join(OUT, "data", "similarity_profile.tsv"), profile_rows, ["sample", "group", "confidence", "x", "similarity"])
    write_tsv(
        os.path.join(OUT, "data", "target_regions.tsv"),
        region_rows,
        [
            "sample",
            "group",
            "confidence",
            "status",
            "ortholog_id",
            "transcript",
            "gene",
            "target_coord",
            "strand",
            "orientation",
            "target_len",
            "ref_global_identity",
            "aligned_columns",
            "blast_pident",
            "blast_qcovs",
            "reciprocal_best_msa_hit",
            "strict_ortholog",
            "fasta",
            "note",
        ],
    )

    gene_rows = []
    for feature, s, e in EXONS:
        gene_rows.append({"feature": feature, "start": s - REF_START + 1, "end": e - REF_START + 1, "type": "exon"})
    write_tsv(os.path.join(OUT, "data", "gene_model.tsv"), gene_rows, ["feature", "start", "end", "type"])

    event_rows = []
    for name, s, e, label in EVENTS:
        event_rows.append({"event": name, "start": s - REF_START + 1, "end": e - REF_START + 1, "label": label})
    write_tsv(os.path.join(OUT, "data", "event_regions.tsv"), event_rows, ["event", "start", "end", "label"])

    summary = os.path.join(OUT, "summary", "data_preparation_summary.txt")
    with open(summary, "w") as handle:
        handle.write(f"Reference: {REF_GENOME} {CHROM}:{REF_START}-{REF_END}, length={len(ref_seq)} bp\n")
        handle.write(f"Ortholog summary: {ORTHO_SUMMARY}\n")
        handle.write(f"Samples requested: {len(SAMPLES)}\n")
        handle.write(f"Samples plotted: {len(set(row['sample'] for row in profile_rows))}\n")
        handle.write(f"Profile rows: {len(profile_rows)}\n")
    print("OUT", OUT)
    print("PROFILE", os.path.join(OUT, "data", "similarity_profile.tsv"))
    print("REGIONS", os.path.join(OUT, "data", "target_regions.tsv"))


if __name__ == "__main__":
    main()
