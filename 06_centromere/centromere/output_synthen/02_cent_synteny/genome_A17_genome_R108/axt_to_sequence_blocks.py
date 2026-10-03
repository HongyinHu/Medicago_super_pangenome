#!/usr/bin/env python3
import argparse
import re


HEADER_RE = re.compile(r"^\d+\s+")


def strip_db(name):
    if "." in name:
        return name.split(".", 1)[1]
    return name


def main():
    parser = argparse.ArgumentParser(
        description="Extract coordinate blocks from UCSC AXT without loading alignment sequences."
    )
    parser.add_argument("--axt", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    n = 0
    with open(args.axt) as inp, open(args.out, "w") as out:
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
            if not HEADER_RE.match(line):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 9:
                continue

            block_id, target_name, target_start, target_end, query_name, query_start, query_end, strand, score = fields[:9]
            target_chr = strip_db(target_name)
            query_chr = strip_db(query_name)
            target_start = int(target_start)
            target_end = int(target_end)
            query_start = int(query_start)
            query_end = int(query_end)

            out.write(
                "\t".join(
                    map(
                        str,
                        [
                            query_chr,
                            query_start,
                            query_end,
                            target_chr,
                            target_start,
                            target_end,
                            strand,
                            score,
                            abs(query_end - query_start),
                            abs(target_end - target_start),
                            block_id,
                        ],
                    )
                )
                + "\n"
            )
            n += 1

    print(f"wrote {n} blocks to {args.out}")


if __name__ == "__main__":
    main()
