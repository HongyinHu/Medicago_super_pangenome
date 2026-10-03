#!/usr/bin/env python3
import argparse
import os


def load_manifest(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            rows.append({name: fields[i] for name, i in idx.items()})
    return rows


def main():
    parser = argparse.ArgumentParser(description="Concatenate non-overlapping Sniffles2 chunk genotype VCFs.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--chunk-dir", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--ref-alias", required=True)
    parser.add_argument("--out-vcf", required=True)
    args = parser.parse_args()

    rows = load_manifest(args.manifest)
    wrote_header = False
    sample_columns = None
    record_count = 0
    os.makedirs(os.path.dirname(args.out_vcf), exist_ok=True)
    with open(args.out_vcf, "w", encoding="utf-8") as out:
        for row in rows:
            chunk_id = row["chunk_id"]
            chunk_vcf = os.path.join(args.chunk_dir, f"{chunk_id}.vcf")
            if not os.path.exists(chunk_vcf):
                raise FileNotFoundError(chunk_vcf)
            with open(chunk_vcf, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith("##"):
                        if not wrote_header:
                            out.write(line)
                    elif line.startswith("#CHROM"):
                        cols = line.rstrip("\n").split("\t")
                        if len(cols) < 10:
                            raise ValueError(f"No sample column in {chunk_vcf}")
                        this_sample = cols[9:]
                        if sample_columns is None:
                            sample_columns = this_sample
                        elif sample_columns != this_sample:
                            raise ValueError(f"Sample columns differ in {chunk_vcf}: {this_sample} != {sample_columns}")
                        if not wrote_header:
                            out.write(line)
                            wrote_header = True
                    else:
                        if not wrote_header:
                            raise ValueError(f"Missing VCF header before records in {chunk_vcf}")
                        out.write(line)
                        record_count += 1
    print(f"records={record_count}\tout={args.out_vcf}")


if __name__ == "__main__":
    main()
