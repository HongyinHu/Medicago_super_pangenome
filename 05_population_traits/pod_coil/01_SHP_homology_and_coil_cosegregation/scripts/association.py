#!/usr/bin/env python3
import argparse
import csv
import math
from collections import Counter
from pathlib import Path


def comb(n, k):
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))


def fisher_two_sided(a, b, c, d):
    r1, r2, c1 = a + b, c + d, a + c
    n = r1 + r2
    lo = max(0, r1 - (n - c1))
    hi = min(r1, c1)
    den = comb(n, r1)
    def prob(x):
        return comb(c1, x) * comb(n - c1, r1 - x) / den
    observed = prob(a)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1) if prob(x) <= observed + 1e-15))


def bh(values):
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    out = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        idx = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[idx] * n / rank)
        out[idx] = min(1.0, running)
    return out


def parse_fasta_alignment(path):
    records = {}
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            elif name is not None:
                seq.append(line)
    if name is not None:
        records[name] = "".join(seq)
    return records


def cross_species(alignment, phenotype_file, output):
    phen = {}
    with open(phenotype_file) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2 and row[0] not in {"genome_A17", "genome_474b"}:
                if row[1] == "Spiral": phen[row[0]] = 1
                elif row[1] == "Non-Spiral": phen[row[0]] = 0
    seqs = parse_fasta_alignment(alignment)
    names = sorted(set(seqs) & set(phen))
    if not names:
        Path(output).write_text("alignment_column\talleles\tn\tfisher_p\tbest_match_fraction\torientation\n")
        return
    length = len(seqs[names[0]])
    reference_name = "genome_Msa" if "genome_Msa" in seqs else names[0]
    reference_pos = 0
    rows = []
    for pos in range(length):
        reference_residue = seqs[reference_name][pos]
        if reference_residue != "-":
            reference_pos += 1
        alleles = Counter(seqs[n][pos] for n in names if seqs[n][pos] not in {"-", "X", "N", "?"})
        if len(alleles) != 2:
            continue
        major = alleles.most_common(1)[0][0]
        usable = [n for n in names if seqs[n][pos] in alleles]
        if len(usable) < 6:
            continue
        a = sum(seqs[n][pos] != major and phen[n] == 1 for n in usable)
        b = sum(seqs[n][pos] != major and phen[n] == 0 for n in usable)
        c = sum(seqs[n][pos] == major and phen[n] == 1 for n in usable)
        d = sum(seqs[n][pos] == major and phen[n] == 0 for n in usable)
        p = fisher_two_sided(a, b, c, d)
        match_alt_spiral = (a + d) / len(usable)
        match_alt_nonspiral = (b + c) / len(usable)
        best = max(match_alt_spiral, match_alt_nonspiral)
        orient = "minor_allele_spiral" if match_alt_spiral >= match_alt_nonspiral else "minor_allele_nonspiral"
        species_alleles = ";".join(f"{n}:{seqs[n][pos]}:{'Spiral' if phen[n] else 'Non-Spiral'}" for n in usable)
        rows.append([pos + 1, reference_pos if reference_residue != "-" else "", reference_residue, ",".join(f"{k}:{v}" for k, v in sorted(alleles.items())), len(usable), a, b, c, d, p, best, orient, species_alleles])
    qvals = bh([r[9] for r in rows]) if rows else []
    with open(output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["alignment_column", "reference_aa_position", "reference_residue", "alleles", "n", "minor_spiral", "minor_nonspiral", "major_spiral", "major_nonspiral", "fisher_p", "BH_FDR", "best_match_fraction", "orientation", "species_alleles"])
        for r, q in zip(rows, qvals):
            writer.writerow(r[:10] + [q] + r[10:])


def parse_gt(gt):
    gt = gt.split(":", 1)[0]
    if gt in {".", "./.", ".|."} or "." in gt:
        return None
    alleles = gt.replace("|", "/").split("/")
    return 1 if any(x != "0" for x in alleles) else 0


def population(query_tsv, samples_file, phenotype_file, output, gene):
    samples = [x.strip() for x in Path(samples_file).read_text().splitlines() if x.strip()]
    phen = {}
    with open(phenotype_file) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            val = row.get("Pod_spiral_binary_strict_spiral1_nonspiral0_weakNA", "")
            if val in {"0", "1"}:
                phen[row["Sample_ID"]] = int(val)
    rows = []
    with open(query_tsv) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt = f[:4]
            gts = f[4:]
            use = []
            for sample, gt in zip(samples, gts):
                dosage = parse_gt(gt)
                if sample in phen and dosage is not None:
                    use.append((dosage, phen[sample]))
            if len(use) < 20:
                continue
            a = sum(g == 1 and p == 1 for g, p in use)
            b = sum(g == 1 and p == 0 for g, p in use)
            c = sum(g == 0 and p == 1 for g, p in use)
            d = sum(g == 0 and p == 0 for g, p in use)
            pval = fisher_two_sided(a, b, c, d)
            m1 = (a + d) / len(use)
            m2 = (b + c) / len(use)
            rows.append([gene, chrom, int(pos), ref, alt, len(use), a, b, c, d, pval, max(m1, m2), "ALT_spiral" if m1 >= m2 else "ALT_nonspiral"])
    qvals = bh([r[10] for r in rows]) if rows else []
    with open(output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["gene", "chrom", "pos", "ref", "alt", "n", "ALT_spiral", "ALT_nonspiral", "REF_spiral", "REF_nonspiral", "fisher_p", "BH_FDR", "best_match_fraction", "orientation"])
        for r, q in zip(rows, qvals):
            writer.writerow(r[:11] + [q] + r[11:])


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode")
    cross = sub.add_parser("cross")
    cross.add_argument("--alignment", required=True)
    cross.add_argument("--phenotype", required=True)
    cross.add_argument("--output", required=True)
    pop = sub.add_parser("population")
    pop.add_argument("--query", required=True)
    pop.add_argument("--samples", required=True)
    pop.add_argument("--phenotype", required=True)
    pop.add_argument("--output", required=True)
    pop.add_argument("--gene", required=True)
    args = parser.parse_args()
    if args.mode is None:
        parser.error("a mode is required: cross or population")
    if args.mode == "cross":
        cross_species(args.alignment, args.phenotype, args.output)
    else:
        population(args.query, args.samples, args.phenotype, args.output, args.gene)


if __name__ == "__main__":
    main()
