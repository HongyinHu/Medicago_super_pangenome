#!/usr/bin/env python3
import argparse
from pathlib import Path


def read_fasta(path):
    records = {}
    header = None
    seq = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = (header, "".join(seq))
                header = line[1:]
                seq = []
            elif header is not None:
                seq.append(line.strip())
    if header is not None:
        records[header.split()[0]] = (header, "".join(seq))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--ids", required=True, help="Comma list or file with one ID per line")
    parser.add_argument("--output", required=True)
    parser.add_argument("--rename-prefix", default="")
    args = parser.parse_args()

    ids_path = Path(args.ids)
    if ids_path.exists():
        wanted = [x.strip().split()[0] for x in ids_path.read_text().splitlines() if x.strip()]
    else:
        wanted = [x.strip() for x in args.ids.split(",") if x.strip()]

    records = read_fasta(args.fasta)
    missing = []
    with open(args.output, "w", encoding="utf-8") as out:
        for ident in wanted:
            matches = [key for key in records if key == ident or key.startswith(ident + ".")]
            if not matches:
                missing.append(ident)
                continue
            key = sorted(matches)[0]
            name = f"{args.rename_prefix}{key}" if args.rename_prefix else key
            out.write(f">{name}\n")
            seq = records[key][1]
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    if missing:
        raise SystemExit("Missing FASTA IDs: " + ",".join(missing))


if __name__ == "__main__":
    main()
