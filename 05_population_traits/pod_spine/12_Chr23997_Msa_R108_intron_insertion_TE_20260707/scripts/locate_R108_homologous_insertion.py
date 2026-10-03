#!/usr/bin/env python3
import os
import re
from collections import Counter, defaultdict

BASE = "path/to/project/N_4.pod_spiny"
OUT = os.path.join(BASE, "12_Chr23997_Msa_R108_intron_insertion_TE_20260707")

MSA_FASTA = os.path.join(BASE, "00_data/1.reference_genome/genome_Msa.fa")
R108_FASTA = os.path.join(BASE, "00_data/1.reference_genome/genome_R108.fa")
CHROM = "Chr4"
INS_START = 88980831
INS_END = 88981052
FLANK = 1000
SEED_K = 25
SEED_STEP = 5


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


def extract(fasta, seqid, start, end):
    for name, seq in fasta_iter(fasta):
        if name == seqid:
            return seq[start - 1 : end]
    raise RuntimeError(f"{seqid} not found in {fasta}")


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def add_seeds(label, seq, offset, seed_map):
    for i in range(0, len(seq) - SEED_K + 1, SEED_STEP):
        seed = seq[i : i + SEED_K]
        if "N" in seed:
            continue
        seed_map[seed].append((label, offset + i, "+"))
        seed_map[revcomp(seed)].append((label, offset + i, "-"))


def smith_waterman(a, b):
    # Small local alignment for flank-to-window placement.
    match = 2
    mismatch = -2
    gap = -3
    n = len(a)
    m = len(b)
    prev = [0] * (m + 1)
    best = (0, 0, 0)
    pointer_rows = []
    # Store traceback directions compactly for these small windows.
    ptr = [[0] * (m + 1) for _ in range(n + 1)]
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
    matches = 0
    aligned = 0
    while i > 0 and j > 0 and ptr[i][j] != 0:
        d = ptr[i][j]
        if d == 1:
            matches += 1 if a[i - 1] == b[j - 1] else 0
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
        "matches": matches,
        "aligned_length": aligned,
        "identity": matches / aligned if aligned else 0,
    }


def main():
    os.makedirs(os.path.join(OUT, "results"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "query"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "summary"), exist_ok=True)

    context_start = INS_START - FLANK
    context_end = INS_END + FLANK
    msa_context = extract(MSA_FASTA, CHROM, context_start, context_end)
    left_flank = msa_context[:FLANK]
    insertion = msa_context[FLANK : FLANK + (INS_END - INS_START + 1)]
    right_flank = msa_context[FLANK + len(insertion) :]

    seed_map = defaultdict(list)
    add_seeds("left_flank", left_flank, 0, seed_map)
    add_seeds("insertion", insertion, FLANK, seed_map)
    add_seeds("right_flank", right_flank, FLANK + len(insertion), seed_map)
    seed_set = set(seed_map)
    bins = Counter()
    detail = defaultdict(list)
    context_len = len(msa_context)

    for seqid, seq in fasta_iter(R108_FASTA):
        for pos in range(0, len(seq) - SEED_K + 1):
            token = seq[pos : pos + SEED_K]
            if token not in seed_set:
                continue
            for label, qpos, strand in seed_map[token]:
                if strand != "+":
                    continue
                predicted_start = pos - qpos + 1
                if predicted_start < 1:
                    continue
                key = (seqid, round(predicted_start / 100) * 100)
                bins[key] += 1
                if len(detail[key]) < 200:
                    detail[key].append((label, qpos + 1, pos + 1, predicted_start))

    best = bins.most_common(10)
    seed_tsv = os.path.join(OUT, "results", "R108_homologous_window_seed_hits.tsv")
    with open(seed_tsv, "w") as handle:
        handle.write("rank\tseqid\tpredicted_start_bin\tseed_hits\tlabels\n")
        for rank, ((seqid, start_bin), count) in enumerate(best, 1):
            labels = Counter(x[0] for x in detail[(seqid, start_bin)])
            handle.write(f"{rank}\t{seqid}\t{start_bin}\t{count}\t{dict(labels)}\n")

    if not best:
        raise RuntimeError("No R108 homologous window found by seed scan")

    best_seqid, best_bin = best[0][0]
    # Use median predicted start from the best bin details.
    starts = sorted(x[3] for x in detail[(best_seqid, best_bin)])
    predicted_start = starts[len(starts) // 2]
    win_start = max(1, predicted_start - 500)
    win_end = predicted_start + context_len + 500
    r108_window = extract(R108_FASTA, best_seqid, win_start, win_end)

    left_aln = smith_waterman(left_flank, r108_window)
    right_aln = smith_waterman(right_flank, r108_window)

    r108_ins_start = win_start + left_aln["target_end"]
    r108_ins_end = win_start + right_aln["target_start"] - 2
    if r108_ins_start <= r108_ins_end:
        r108_insertion = extract(R108_FASTA, best_seqid, r108_ins_start, r108_ins_end)
    else:
        r108_insertion = ""

    out_fa = os.path.join(OUT, "query", "Chr23997_Msa_R108_homologous_intron_insertions.fa")
    with open(out_fa, "w") as handle:
        handle.write(f">genome_Msa_Chr23997_intron_insertion|{CHROM}:{INS_START}-{INS_END}|len={len(insertion)}\n")
        handle.write(insertion + "\n")
        handle.write(
            f">genome_R108_Chr23997_homologous_intron_insertion|{best_seqid}:{r108_ins_start}-{r108_ins_end}|len={len(r108_insertion)}\n"
        )
        handle.write(r108_insertion + "\n")

    summary = os.path.join(OUT, "results", "R108_homologous_insertion_location.tsv")
    with open(summary, "w") as handle:
        handle.write("item\tvalue\n")
        handle.write(f"msa_context\t{CHROM}:{context_start}-{context_end}\n")
        handle.write(f"msa_insertion\t{CHROM}:{INS_START}-{INS_END}\n")
        handle.write(f"msa_insertion_len\t{len(insertion)}\n")
        handle.write(f"r108_seed_window\t{best_seqid}:{win_start}-{win_end}\n")
        handle.write(f"r108_predicted_context_start\t{predicted_start}\n")
        handle.write(f"r108_insertion\t{best_seqid}:{r108_ins_start}-{r108_ins_end}\n")
        handle.write(f"r108_insertion_len\t{len(r108_insertion)}\n")
        handle.write(f"left_flank_alignment\t{left_aln}\n")
        handle.write(f"right_flank_alignment\t{right_aln}\n")
        handle.write(f"query_fasta\t{out_fa}\n")
        handle.write(f"seed_hits\t{seed_tsv}\n")

    print("QUERY_FASTA", out_fa)
    print("SUMMARY", summary)
    print("SEED_HITS", seed_tsv)
    print("BEST_WINDOW", f"{best_seqid}:{win_start}-{win_end}", "seed_hits", best[0][1])
    print("R108_INSERTION", f"{best_seqid}:{r108_ins_start}-{r108_ins_end}", "len", len(r108_insertion))
    print("LEFT_ALN", left_aln)
    print("RIGHT_ALN", right_aln)
    print("MSA_INSERTION_LEN", len(insertion))


if __name__ == "__main__":
    main()
