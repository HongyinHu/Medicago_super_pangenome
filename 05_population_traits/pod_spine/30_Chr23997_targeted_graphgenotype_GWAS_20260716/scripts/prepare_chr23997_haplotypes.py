#!/usr/bin/env python3
"""Create reference and exact-DEL local haplotypes for Chr23997."""

import argparse
import importlib.util
from pathlib import Path

import pysam


MODULE = Path(__file__).with_name("chr23997_haplotype_genotyper.py")
SPEC = importlib.util.spec_from_file_location("chr23997_haplotype_genotyper", MODULE)
GENOTYPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENOTYPER)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chrom", default="Chr4")
    parser.add_argument("--pos", type=int, default=88980831)
    parser.add_argument("--end", type=int, default=88981052)
    parser.add_argument("--window", type=int, default=2000)
    args = parser.parse_args()

    start0 = args.pos - args.window
    end0 = args.end + args.window
    with pysam.FastaFile(args.reference) as fasta:
        sequence = fasta.fetch(args.chrom, start0, end0).upper()
    alt, _ = GENOTYPER.build_alt_haplotype(
        sequence, args.window, args.window + args.end - args.pos
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        handle.write(f">Chr23997_REF|{args.chrom}:{start0 + 1}-{end0}\n{sequence}\n")
        handle.write(
            f">Chr23997_DEL221|{args.chrom}:{start0 + 1}-{end0}|DEL:{args.pos + 1}-{args.end}\n{alt}\n"
        )


if __name__ == "__main__":
    main()
