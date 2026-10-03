#!/usr/bin/env python3
import argparse


def strip_db(name):
    return name.split(".", 1)[1] if "." in name else name


def main():
    parser = argparse.ArgumentParser(description="Convert minimap2 PAF records to plotting block TSV.")
    parser.add_argument("--paf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    n = 0
    with open(args.paf) as inp, open(args.out, "w") as out:
        out.write(
            "\t".join(
                [
                    "query_chr",
                    "query_start",
                    "query_end",
                    "target_chr",
                    "target_start",
                    "target_end",
                    "strand",
                    "score",
                    "query_len",
                    "target_len",
                    "block_id",
                ]
            )
            + "\n"
        )
        for line in inp:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, _qlen, qstart, qend, strand, tname, _tlen, tstart, tend, nmatch, block_len, _mapq = fields[:12]
            qstart = int(qstart)
            qend = int(qend)
            tstart = int(tstart)
            tend = int(tend)
            out.write(
                "\t".join(
                    map(
                        str,
                        [
                            strip_db(qname),
                            qstart,
                            qend,
                            strip_db(tname),
                            tstart,
                            tend,
                            strand,
                            nmatch,
                            abs(qend - qstart),
                            abs(tend - tstart),
                            n,
                        ],
                    )
                )
                + "\n"
            )
            n += 1
    print(f"wrote {n} PAF blocks to {args.out}")


if __name__ == "__main__":
    main()
