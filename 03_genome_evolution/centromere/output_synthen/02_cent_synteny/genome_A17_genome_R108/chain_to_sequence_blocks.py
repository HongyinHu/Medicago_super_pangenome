#!/usr/bin/env python3
import argparse


def strip_db(name):
    return name.split(".", 1)[1] if "." in name else name


def plus_query_interval(q_pos, size, q_size, q_strand):
    if q_strand == "+":
        return q_pos, q_pos + size
    return q_size - (q_pos + size), q_size - q_pos


def main():
    parser = argparse.ArgumentParser(description="Convert UCSC chain blocks to plotting block TSV.")
    parser.add_argument("--chain", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    n = 0
    current = None
    t_pos = None
    q_pos = None
    with open(args.chain) as inp, open(args.out, "w") as out:
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
        for raw in inp:
            line = raw.strip()
            if not line:
                current = None
                t_pos = None
                q_pos = None
                continue
            fields = line.split()
            if fields[0] == "chain":
                if len(fields) < 13:
                    continue
                current = {
                    "score": fields[1],
                    "target_chr": strip_db(fields[2]),
                    "query_chr": strip_db(fields[7]),
                    "query_size": int(fields[8]),
                    "query_strand": fields[9],
                }
                t_pos = int(fields[5])
                q_pos = int(fields[10])
                continue
            if current is None:
                continue

            size = int(fields[0])
            target_start = t_pos
            target_end = t_pos + size
            query_start, query_end = plus_query_interval(
                q_pos,
                size,
                current["query_size"],
                current["query_strand"],
            )
            out.write(
                "\t".join(
                    map(
                        str,
                        [
                            current["query_chr"],
                            query_start,
                            query_end,
                            current["target_chr"],
                            target_start,
                            target_end,
                            current["query_strand"],
                            current["score"],
                            abs(query_end - query_start),
                            abs(target_end - target_start),
                            n,
                        ],
                    )
                )
                + "\n"
            )
            n += 1

            if len(fields) >= 3:
                t_pos = target_end + int(fields[1])
                q_pos = q_pos + size + int(fields[2])

    print(f"wrote {n} chain blocks to {args.out}")


if __name__ == "__main__":
    main()
