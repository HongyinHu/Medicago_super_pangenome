#!/usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict


def read_fai(path):
    lengths = {}
    order = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            lengths[chrom] = int(fields[1])
            order.append(chrom)
    return order, lengths


def open_text(path):
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(info):
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def record_span(fields):
    chrom = fields[0]
    pos = int(fields[1])
    info = parse_info(fields[7] if len(fields) > 7 else "")
    end = pos
    if "END" in info:
        try:
            end = int(info["END"])
        except ValueError:
            end = pos
    elif "SVLEN" in info:
        m = re.match(r"^-?\d+", info["SVLEN"])
        if m:
            end = max(pos, pos + abs(int(m.group(0))))
    return chrom, min(pos, end), max(pos, end)


def flush_chunk(args, header, records, chunk_index, chrom_lengths, manifest):
    if not records:
        return chunk_index
    first = records[0].rstrip("\n").split("\t")
    chrom, start, end = record_span(first)
    min_pos = start
    max_end = end
    for rec in records[1:]:
        fields = rec.rstrip("\n").split("\t")
        _, rec_start, rec_end = record_span(fields)
        min_pos = min(min_pos, rec_start)
        max_end = max(max_end, rec_end)
    chunk_id = f"{args.ref_alias}.chunk{chunk_index:05d}.{chrom}.{min_pos}_{max_end}"
    chunk_vcf = os.path.join(args.out_dir, f"{chunk_id}.vcf")
    chunk_bed = os.path.join(args.out_dir, f"{chunk_id}.bed")
    region_start = max(0, min_pos - 1 - args.padding)
    region_end = min(chrom_lengths.get(chrom, max_end + args.padding), max_end + args.padding)
    if region_end <= region_start:
        region_end = region_start + 1
    with open(chunk_vcf, "w", encoding="utf-8") as out:
        out.writelines(header)
        out.writelines(records)
    with open(chunk_bed, "w", encoding="utf-8") as out:
        out.write(f"{chrom}\t{region_start}\t{region_end}\n")
    manifest.write(
        "\t".join(
            [
                chunk_id,
                chrom,
                str(region_start),
                str(region_end),
                str(len(records)),
                chunk_vcf,
                chunk_bed,
            ]
        )
        + "\n"
    )
    return chunk_index + 1


def main():
    parser = argparse.ArgumentParser(description="Split a sorted discovery VCF into Sniffles2 genotype chunks.")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--ref-alias", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--max-records", type=int, default=5000)
    parser.add_argument("--padding", type=int, default=20000)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    chrom_order, chrom_lengths = read_fai(args.fai)
    chrom_rank = defaultdict(lambda: 10**9)
    for i, chrom in enumerate(chrom_order):
        chrom_rank[chrom] = i

    header = []
    current = []
    current_chrom = None
    chunk_index = 1
    total_records = 0
    with open(args.manifest, "w", encoding="utf-8") as manifest:
        manifest.write("chunk_id\tchrom\tstart\tend\trecord_count\tvcf\tbed\n")
        with open_text(args.vcf) as handle:
            for line in handle:
                if line.startswith("#"):
                    header.append(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                chrom = fields[0]
                if current and (chrom != current_chrom or len(current) >= args.max_records):
                    chunk_index = flush_chunk(args, header, current, chunk_index, chrom_lengths, manifest)
                    current = []
                current_chrom = chrom
                current.append(line)
                total_records += 1
            if current:
                chunk_index = flush_chunk(args, header, current, chunk_index, chrom_lengths, manifest)
    print(f"chunks={chunk_index - 1}\trecords={total_records}\tmanifest={args.manifest}")


if __name__ == "__main__":
    main()
