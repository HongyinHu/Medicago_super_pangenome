#!/usr/bin/env python3
import csv
import os
import re
from collections import Counter, defaultdict

BASE = "path/to/project/N_4.pod_spiny"
REFDIR = os.path.join(BASE, "00_data/1.reference_genome")
OUT = os.path.join(BASE, "13_Chr23997_spiny_intron_window_comparison_20260708")

CHROM = "Chr4"
FLANK = 1000
PAD = 500
SEED_K = 31
SEED_STEP = 10
BIN_SIZE = 100

WINDOWS = [
    {
        "window": "core_222",
        "label": "Msa_core_insertion_222bp",
        "start": 88980831,
        "end": 88981052,
        "description": "Msa core intron insertion absent/deleted in M22",
    },
    {
        "window": "core_plus_100",
        "label": "Msa_core_insertion_plus_100bp",
        "start": 88980731,
        "end": 88981152,
        "description": "core insertion with 100 bp flanks on both sides",
    },
    {
        "window": "intron2_full",
        "label": "Msa_second_intron_full",
        "start": 88980584,
        "end": 88981373,
        "description": "whole intron between exon2 and exon3 of Chr23997.1",
    },
]

TARGETS = [
    {
        "genome": "genome_Msa",
        "species": "Medicago_sativa",
        "fasta": os.path.join(REFDIR, "genome_Msa.fa"),
        "role": "reference_window",
    },
    {
        "genome": "genome_R108",
        "species": "Medicago_truncatula_R108",
        "fasta": os.path.join(REFDIR, "genome_R108.fa"),
        "role": "spiny_primary",
    },
    {
        "genome": "genome_410",
        "species": "Medicago_praecox",
        "fasta": os.path.join(REFDIR, "genome_410.fa"),
        "role": "spiny_primary",
    },
    {
        "genome": "genome_474a",
        "species": "Medicago_carstiensis",
        "fasta": os.path.join(REFDIR, "genome_474a.fa"),
        "role": "spiny_primary_474_assembly_a",
    },
    {
        "genome": "genome_474b",
        "species": "Medicago_carstiensis",
        "fasta": os.path.join(REFDIR, "genome_474b.fa"),
        "role": "spiny_primary_474_assembly_b",
    },
    {
        "genome": "genome_Mpo",
        "species": "Medicago_polymorpha",
        "fasta": os.path.join(REFDIR, "genome_Mpo.fa"),
        "role": "spiny_primary",
    },
]


def clean_seq(line):
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
                chunks.append(clean_seq(line.strip()))
    if not chunks:
        raise RuntimeError(f"{seqid} not found in {fasta}")
    return "".join(chunks)


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def subseq(seq, start, end):
    start = max(1, start)
    end = min(len(seq), end)
    if start > end:
        return ""
    return seq[start - 1 : end]


def wrap(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def add_flank_seeds(seed_map, win_name, label, seq, offset):
    for i in range(0, len(seq) - SEED_K + 1, SEED_STEP):
        token = seq[i : i + SEED_K]
        if "N" in token:
            continue
        seed_map[token].append((win_name, label, offset + i))


def smith_waterman(a, b):
    match = 2
    mismatch = -2
    gap = -3
    n = len(a)
    m = len(b)
    prev = [0] * (m + 1)
    ptr = [[0] * (m + 1) for _ in range(n + 1)]
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
    matches = 0
    aligned = 0
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
        "matches": matches,
        "aligned_length": aligned,
        "identity": matches / aligned if aligned else 0,
        "query_coverage": aligned / max(1, len(a)),
    }


def scan_target_chr4(seq, seed_map):
    bins = Counter()
    details = defaultdict(list)
    seed_set = set(seed_map)
    for pos in range(0, len(seq) - SEED_K + 1):
        token = seq[pos : pos + SEED_K]
        if token not in seed_set:
            continue
        for win_name, label, qpos in seed_map[token]:
            predicted_start = pos - qpos + 1
            if predicted_start < 1:
                continue
            start_bin = round(predicted_start / BIN_SIZE) * BIN_SIZE
            key = (win_name, start_bin)
            bins[key] += 1
            if len(details[key]) < 300:
                details[key].append((label, qpos + 1, pos + 1, predicted_start))
    return bins, details


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main():
    os.makedirs(os.path.join(OUT, "scripts"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "query"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "results"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "summary"), exist_ok=True)

    msa_seq = read_seqid(os.path.join(REFDIR, "genome_Msa.fa"), CHROM)
    win_info = {}
    seed_map = defaultdict(list)

    ref_fa = os.path.join(OUT, "query", "Chr23997_Msa_reference_windows.fa")
    with open(ref_fa, "w") as handle:
        for win in WINDOWS:
            wseq = subseq(msa_seq, win["start"], win["end"])
            context_start = win["start"] - FLANK
            context_end = win["end"] + FLANK
            context = subseq(msa_seq, context_start, context_end)
            left = context[:FLANK]
            middle = subseq(msa_seq, win["start"], win["end"])
            right = context[FLANK + len(middle) :]
            win_info[win["window"]] = {
                **win,
                "seq": wseq,
                "context_start": context_start,
                "context_end": context_end,
                "context_len": len(context),
                "left": left,
                "right": right,
            }
            add_flank_seeds(seed_map, win["window"], "left_flank", left, 0)
            add_flank_seeds(seed_map, win["window"], "right_flank", right, FLANK + len(middle))
            handle.write(f">{win['window']}|{win['label']}|{CHROM}:{win['start']}-{win['end']}|len={len(wseq)}\n")
            handle.write(wrap(wseq) + "\n")

    rows = []
    seed_rows = []
    all_fa = os.path.join(OUT, "query", "Chr23997_spiny_homologous_windows.fa")
    with open(all_fa, "w") as fa:
        for target in TARGETS:
            if not os.path.exists(target["fasta"]):
                for win in WINDOWS:
                    rows.append(
                        {
                            "genome": target["genome"],
                            "species": target["species"],
                            "role": target["role"],
                            "window": win["window"],
                            "status": "missing_fasta",
                            "fasta": target["fasta"],
                        }
                    )
                continue

            tseq = read_seqid(target["fasta"], CHROM)
            print("LOADED", target["genome"], CHROM, len(tseq), flush=True)

            if target["genome"] == "genome_Msa":
                for win in WINDOWS:
                    info = win_info[win["window"]]
                    wseq = info["seq"]
                    fa.write(
                        f">{target['genome']}|{win['window']}|{CHROM}:{win['start']}-{win['end']}|len={len(wseq)}|status=reference\n"
                    )
                    fa.write(wrap(wseq) + "\n")
                    rows.append(
                        {
                            "genome": target["genome"],
                            "species": target["species"],
                            "role": target["role"],
                            "window": win["window"],
                            "status": "reference",
                            "query_coord": f"{CHROM}:{win['start']}-{win['end']}",
                            "target_coord": f"{CHROM}:{win['start']}-{win['end']}",
                            "query_len": len(wseq),
                            "target_len": len(wseq),
                            "length_delta_vs_Msa": 0,
                            "seed_hits": "",
                            "left_identity": 1,
                            "right_identity": 1,
                            "window_local_identity_vs_Msa": 1,
                            "window_local_coverage_vs_Msa": 1,
                            "notes": "Msa reference interval",
                            "fasta": target["fasta"],
                        }
                    )
                continue

            bins, details = scan_target_chr4(tseq, seed_map)
            for key, count in bins.items():
                labels = Counter(x[0] for x in details[key])
                seed_rows.append(
                    {
                        "genome": target["genome"],
                        "window": key[0],
                        "target_seqid": CHROM,
                        "predicted_context_start_bin": key[1],
                        "seed_hits": count,
                        "labels": dict(labels),
                    }
                )

            for win in WINDOWS:
                info = win_info[win["window"]]
                candidates = [(key, count) for key, count in bins.items() if key[0] == win["window"]]
                if not candidates:
                    rows.append(
                        {
                            "genome": target["genome"],
                            "species": target["species"],
                            "role": target["role"],
                            "window": win["window"],
                            "status": "no_seed_hits",
                            "query_coord": f"{CHROM}:{win['start']}-{win['end']}",
                            "query_len": len(info["seq"]),
                            "fasta": target["fasta"],
                        }
                    )
                    continue

                best_key, seed_hits = sorted(candidates, key=lambda x: x[1], reverse=True)[0]
                starts = sorted(x[3] for x in details[best_key])
                predicted_start = starts[len(starts) // 2]
                cand_start = max(1, predicted_start - PAD)
                cand_end = min(len(tseq), predicted_start + info["context_len"] + PAD)
                candidate = subseq(tseq, cand_start, cand_end)

                left_aln = smith_waterman(info["left"], candidate)
                right_aln = smith_waterman(info["right"], candidate)
                target_start = cand_start + left_aln["target_end"]
                target_end = cand_start + right_aln["target_start"] - 2

                status = "ok"
                notes = []
                if target_start > target_end:
                    status = "invalid_flank_order"
                    notes.append("right flank aligns before left flank")
                    target_seq = ""
                else:
                    target_seq = subseq(tseq, target_start, target_end)

                if left_aln["identity"] < 0.75 or right_aln["identity"] < 0.75:
                    notes.append("low flank identity")
                if left_aln["query_coverage"] < 0.5 or right_aln["query_coverage"] < 0.5:
                    notes.append("low flank coverage")

                if target_seq:
                    aln = smith_waterman(info["seq"], target_seq)
                    fa.write(
                        f">{target['genome']}|{win['window']}|{CHROM}:{target_start}-{target_end}|len={len(target_seq)}|status={status}\n"
                    )
                    fa.write(wrap(target_seq) + "\n")
                else:
                    aln = {"identity": 0, "query_coverage": 0, "score": 0, "aligned_length": 0}

                rows.append(
                    {
                        "genome": target["genome"],
                        "species": target["species"],
                        "role": target["role"],
                        "window": win["window"],
                        "status": status,
                        "query_coord": f"{CHROM}:{win['start']}-{win['end']}",
                        "target_coord": f"{CHROM}:{target_start}-{target_end}" if target_seq else "",
                        "query_len": len(info["seq"]),
                        "target_len": len(target_seq),
                        "length_delta_vs_Msa": len(target_seq) - len(info["seq"]),
                        "seed_hits": seed_hits,
                        "predicted_context_start": predicted_start,
                        "candidate_window": f"{CHROM}:{cand_start}-{cand_end}",
                        "left_identity": f"{left_aln['identity']:.4f}",
                        "left_coverage": f"{left_aln['query_coverage']:.4f}",
                        "right_identity": f"{right_aln['identity']:.4f}",
                        "right_coverage": f"{right_aln['query_coverage']:.4f}",
                        "window_local_identity_vs_Msa": f"{aln['identity']:.4f}",
                        "window_local_coverage_vs_Msa": f"{aln['query_coverage']:.4f}",
                        "window_local_alignment_len": aln["aligned_length"],
                        "notes": "; ".join(notes),
                        "fasta": target["fasta"],
                    }
                )

    homology_tsv = os.path.join(OUT, "results", "Chr23997_spiny_window_homology.tsv")
    write_tsv(
        homology_tsv,
        rows,
        [
            "genome",
            "species",
            "role",
            "window",
            "status",
            "query_coord",
            "target_coord",
            "query_len",
            "target_len",
            "length_delta_vs_Msa",
            "seed_hits",
            "predicted_context_start",
            "candidate_window",
            "left_identity",
            "left_coverage",
            "right_identity",
            "right_coverage",
            "window_local_identity_vs_Msa",
            "window_local_coverage_vs_Msa",
            "window_local_alignment_len",
            "notes",
            "fasta",
        ],
    )

    seed_tsv = os.path.join(OUT, "results", "Chr23997_spiny_window_seed_hits.tsv")
    write_tsv(
        seed_tsv,
        sorted(seed_rows, key=lambda r: (r["genome"], r["window"], -int(r["seed_hits"]))),
        ["genome", "window", "target_seqid", "predicted_context_start_bin", "seed_hits", "labels"],
    )

    summary_md = os.path.join(OUT, "summary", "Chr23997_spiny_intron_window_comparison.md")
    with open(summary_md, "w") as handle:
        handle.write("# Chr23997 spiny species intron-window comparison\n\n")
        handle.write("Reference: genome_Msa Chr4 around Chr23997 second intron.\n\n")
        handle.write("## Windows\n\n")
        for win in WINDOWS:
            handle.write(
                f"- {win['window']}: {CHROM}:{win['start']}-{win['end']}, length={win['end'] - win['start'] + 1} bp, {win['description']}\n"
            )
        handle.write("\n## Main homology table\n\n")
        handle.write(f"- `{homology_tsv}`\n")
        handle.write(f"- `{seed_tsv}`\n")
        handle.write(f"- `{all_fa}`\n")
        handle.write(f"- `{ref_fa}`\n\n")
        handle.write("## Extracted window overview\n\n")
        for row in rows:
            handle.write(
                f"- {row.get('genome')} {row.get('window')}: status={row.get('status')}, target={row.get('target_coord','')}, len={row.get('target_len','')}, local_identity={row.get('window_local_identity_vs_Msa','')}, notes={row.get('notes','')}\n"
            )

    print("OUT", OUT)
    print("REFERENCE_FASTA", ref_fa)
    print("HOMOLOGY_TSV", homology_tsv)
    print("SEED_TSV", seed_tsv)
    print("WINDOW_FASTA", all_fa)
    print("SUMMARY", summary_md)


if __name__ == "__main__":
    main()
