#!/usr/bin/env python3
import re
import sys


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: normalize_paftools_maf.py in.maf out.maf")

    in_path, out_path = sys.argv[1], sys.argv[2]
    pat = re.compile(r"^a\s+([0-9.]+)\s*$")

    with open(in_path) as inp, open(out_path, "w") as out:
        for line in inp:
            m = pat.match(line.rstrip("\n"))
            if m:
                out.write(f"a score={m.group(1)}\n")
            else:
                out.write(line)


if __name__ == "__main__":
    main()
