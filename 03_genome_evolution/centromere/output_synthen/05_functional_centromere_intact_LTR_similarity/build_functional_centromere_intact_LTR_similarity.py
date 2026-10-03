#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


CHR_RE = re.compile(r"^Chr[0-9]+$")
HEADER_RE = re.compile(r"^(?P<teid>[^|]+)\|(?P<chr>[^:]+):(?P<start>\d+)\.\.(?P<end>\d+)#(?P<classification>\S+)")


def read_bed(path):
    regions = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if not CHR_RE.match(chrom):
                continue
            regions.setdefault(chrom, []).append((start, end))
    return regions


def overlaps(regions, chrom, start0, end):
    for rstart, rend in regions.get(chrom, []):
        if start0 < rend and end > rstart:
            return True
    return False


def read_pass_list(path):
    rows = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            loc = row.get("#LTR_loc") or row.get("LTR_loc")
            if not loc:
                continue
            rows[loc] = row
    return rows


def read_fasta(path):
    header = None
    seq_parts = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts)


def clean_token(value):
    return re.sub(r"[^A-Za-z0-9_.:-]+", "_", str(value))


def extract_one(genome, label, pass_list, fasta, bed, out_fa, meta_writer):
    pass_rows = read_pass_list(pass_list)
    regions = read_bed(bed)
    kept = 0
    for header, seq in read_fasta(fasta):
        match = HEADER_RE.search(header)
        if not match:
            continue
        chrom = match.group("chr")
        if not CHR_RE.match(chrom):
            continue
        start = int(match.group("start"))
        end = int(match.group("end"))
        start0 = start - 1
        if not overlaps(regions, chrom, start0, end):
            continue
        ltr_loc = f"{chrom}:{start}..{end}"
        pass_row = pass_rows.get(ltr_loc, {})
        identity = pass_row.get("Identity", "")
        insertion_time = pass_row.get("Insertion_Time", "")
        insertion_mya = ""
        try:
            insertion_mya = float(insertion_time) / 1_000_000.0
        except Exception:
            pass
        classification = match.group("classification")
        teid = match.group("teid")
        seq_id = clean_token(f"{genome}__{teid}__{chrom}_{start}_{end}__{classification.replace('/', '_')}")
        out_fa.write(f">{seq_id}\n")
        for i in range(0, len(seq), 80):
            out_fa.write(seq[i : i + 80] + "\n")
        meta_writer.writerow(
            {
                "id": seq_id,
                "genome": genome,
                "label": label,
                "teid": teid,
                "chr": chrom,
                "start": start,
                "end": end,
                "length": len(seq),
                "classification": classification,
                "ltr_identity": identity,
                "insertion_time_years": insertion_time,
                "insertion_time_mya": insertion_mya,
                "original_header": header,
            }
        )
        kept += 1
    return kept


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_fasta = outdir / "functional_centromere.intact_LTRRT.fa"
    out_meta = outdir / "functional_centromere.intact_LTRRT.metadata.tsv"

    fieldnames = [
        "id",
        "genome",
        "label",
        "teid",
        "chr",
        "start",
        "end",
        "length",
        "classification",
        "ltr_identity",
        "insertion_time_years",
        "insertion_time_mya",
        "original_header",
    ]
    counts = []
    with open(out_fasta, "w") as fa, open(out_meta, "w", newline="") as meta:
        writer = csv.DictWriter(meta, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        with open(args.config) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                count = extract_one(
                    row["genome"],
                    row["label"],
                    row["pass_list"],
                    row["intact_fasta"],
                    row["centromere_bed"],
                    fa,
                    writer,
                )
                counts.append((row["genome"], row["label"], count))
                print(f"{row['genome']}: kept {count} functional-centromere intact LTR-RTs")

    with open(outdir / "functional_centromere.intact_LTRRT.counts.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome", "label", "count"])
        writer.writerows(counts)


if __name__ == "__main__":
    main()
