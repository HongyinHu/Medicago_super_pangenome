#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_fasta(path):
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name.split()[0], "".join(seq).upper()
                name = line[1:]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name.split()[0], "".join(seq).upper()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--chunk-size", type=int, default=2000)
    parser.add_argument(
        "--denominator",
        choices=["union", "both"],
        default="union",
        help="union treats gaps as mismatches; both ignores columns where either sequence has a gap.",
    )
    args = parser.parse_args()

    names = []
    seqs = []
    for name, seq in parse_fasta(args.alignment):
        names.append(name)
        seqs.append(seq)
    if not seqs:
        raise SystemExit("No sequences found.")
    lengths = {len(seq) for seq in seqs}
    if len(lengths) != 1:
        raise SystemExit("Alignment sequences have different lengths: {}".format(sorted(lengths)[:10]))

    n = len(seqs)
    aln_len = len(seqs[0])
    arr = np.frombuffer("".join(seqs).encode("ascii"), dtype="S1").reshape(n, aln_len)
    matches = np.zeros((n, n), dtype=np.uint32)
    comparable = np.zeros((n, n), dtype=np.uint32)
    observed = np.unique(arr)
    ignored = {b"-", b".", b"*"}
    alphabet = [char for char in observed if char not in ignored]

    for start in range(0, aln_len, args.chunk_size):
        block = arr[:, start:start + args.chunk_size]
        nongap = (block != b"-") & (block != b".") & (block != b"*")
        if args.denominator == "both":
            comp = nongap.astype(np.uint16).dot(nongap.astype(np.uint16).T)
        else:
            gap = (~nongap).astype(np.uint16)
            gap_both = gap.dot(gap.T)
            comp = np.uint32(block.shape[1]) - gap_both.astype(np.uint32)
        comparable += comp.astype(np.uint32)
        for base in alphabet:
            is_base = (block == base).astype(np.uint16)
            matches += is_base.dot(is_base.T).astype(np.uint32)

    with np.errstate(divide="ignore", invalid="ignore"):
        sim = matches / comparable
    sim[~np.isfinite(sim)] = 0.0
    np.fill_diagonal(sim, 1.0)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    matrix = pd.DataFrame(sim, index=names, columns=names)
    matrix.to_csv(str(out_prefix) + ".matrix_0to1.tsv", sep="\t", float_format="%.6f")
    (matrix * 100.0).to_csv(str(out_prefix) + ".matrix_0to100_TBtools.tsv", sep="\t", float_format="%.2f")

    tri = sim[np.tril_indices(n, -1)]
    print("sequences={}".format(n))
    print("alignment_length={}".format(aln_len))
    print("denominator={}".format(args.denominator))
    print("min={:.6f}\tmedian={:.6f}\tmean={:.6f}\tq90={:.6f}\tmax={:.6f}".format(
        float(np.min(tri)), float(np.median(tri)), float(np.mean(tri)), float(np.quantile(tri, 0.9)), float(np.max(tri))
    ))


if __name__ == "__main__":
    main()
