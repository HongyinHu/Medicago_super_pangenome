#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a biallelic, QC-filtered cohort VCF to GEMMA BIMBAM input."
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--sample-order", required=True)
    parser.add_argument("--geno", required=True)
    parser.add_argument("--anno", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--summary", required=True)
    return parser.parse_args()


def read_order(path):
    samples = [line.strip() for line in open(path, encoding="utf-8") if line.strip()]
    if len(samples) != len(set(samples)):
        raise SystemExit("Duplicate sample IDs in expected sample order")
    return samples


def info_value(record, key, default="NA"):
    value = record.info.get(key)
    if value is None:
        return default
    if isinstance(value, tuple):
        return ",".join(map(str, value))
    return str(value)


def main():
    args = parse_args()
    expected = read_order(args.sample_order)
    vcf = pysam.VariantFile(args.vcf)
    observed = list(vcf.header.samples)
    if observed != expected:
        raise SystemExit("VCF sample order does not match SNP-derived 143-sample order")

    for path in (args.geno, args.anno, args.mapping, args.summary):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    n_records = 0
    n_missing = 0
    n_called = 0
    type_counts = {}

    with open(args.geno, "w", encoding="utf-8", newline="") as geno, \
         open(args.anno, "w", encoding="utf-8", newline="") as anno, \
         gzip.open(args.mapping, "wt", encoding="utf-8", newline="") as mapping:
        map_fields = [
            "marker_id", "chrom", "pos", "end", "original_id", "svtype",
            "svlen", "mateid", "event", "ref", "alt"
        ]
        map_writer = csv.DictWriter(mapping, fieldnames=map_fields, delimiter="\t")
        map_writer.writeheader()

        for record in vcf:
            if len(record.alts or ()) != 1:
                raise SystemExit(f"Non-biallelic record after QC: {record.contig}:{record.pos}")
            n_records += 1
            marker = f"SV{n_records:07d}"
            dosage = []
            for sample in expected:
                gt = record.samples[sample].get("GT")
                if gt is None or any(allele is None for allele in gt):
                    dosage.append("NA")
                    n_missing += 1
                    continue
                if any(allele not in (0, 1) for allele in gt):
                    raise SystemExit(f"Unexpected allele in {marker}, sample {sample}: {gt}")
                dosage.append(str(sum(gt)))
                n_called += 1

            geno.write(",".join([marker, "A", "C", *dosage]) + "\n")
            chrom = record.contig[3:] if record.contig.startswith("Chr") else record.contig
            anno.write(f"{marker},{record.pos},{chrom}\n")
            svtype = info_value(record, "SVTYPE")
            type_counts[svtype] = type_counts.get(svtype, 0) + 1
            end = record.stop if record.stop is not None else record.pos
            map_writer.writerow(
                {
                    "marker_id": marker,
                    "chrom": record.contig,
                    "pos": record.pos,
                    "end": end,
                    "original_id": record.id or ".",
                    "svtype": svtype,
                    "svlen": info_value(record, "SVLEN"),
                    "mateid": info_value(record, "MATEID"),
                    "event": info_value(record, "EVENT"),
                    "ref": record.ref,
                    "alt": record.alts[0],
                }
            )

    total = n_records * len(expected)
    if total != n_called + n_missing:
        raise SystemExit("Internal genotype count mismatch")
    with open(args.summary, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"samples\t{len(expected)}\n")
        handle.write(f"variants\t{n_records}\n")
        handle.write(f"called_genotypes\t{n_called}\n")
        handle.write(f"missing_genotypes\t{n_missing}\n")
        handle.write(f"overall_call_rate\t{(n_called / total if total else 0):.8f}\n")
        for svtype, count in sorted(type_counts.items()):
            handle.write(f"SVTYPE_{svtype}\t{count}\n")

    print(f"Converted {n_records} SVs across {len(expected)} samples")


if __name__ == "__main__":
    main()
