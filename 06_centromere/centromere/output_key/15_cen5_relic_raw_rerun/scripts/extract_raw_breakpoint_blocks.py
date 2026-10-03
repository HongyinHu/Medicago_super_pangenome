import csv
from pathlib import Path


RUN = Path("path/to/project/10.centromere_analysis/output_key/15_cen5_relic_raw_rerun")
PAF = RUN / "paf" / "Mpo_to_R108.raw_asm5.paf"
OUT = RUN / "results" / "raw_minimap_breakpoint_local_blocks.tsv"

WINDOWS = [
    ("Chr3", 18_000_000, 25_000_000, "Chr3_CEN5_to_CEN6_transition"),
    ("Chr3", 38_000_000, 46_500_000, "Chr3_CEN6_to_CEN3_transition"),
    ("Chr5", 18_000_000, 30_000_000, "Chr5_CEN5_to_CEN6_transition"),
    ("Chr5", 54_000_000, 59_500_000, "Chr5_CEN6_to_CEN3_transition"),
]
ANCESTRAL = {"Chr3", "Chr5", "Chr6"}


def overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


rows = []
with PAF.open() as handle:
    for line in handle:
        if not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq = p[:12]
        if qname not in {"Chr3", "Chr5"} or tname not in ANCESTRAL:
            continue
        qstart_i, qend_i = int(qstart), int(qend)
        alen_i = int(alen)
        if alen_i < 3000:
            continue
        ident = int(nmatch) / max(1, alen_i)
        if ident < 0.70:
            continue
        for chrom, w0, w1, label in WINDOWS:
            if qname != chrom:
                continue
            ov = overlap(qstart_i, qend_i, w0, w1)
            if ov <= 0:
                continue
            rows.append(
                {
                    "breakpoint_region": label,
                    "Mpo_chr": qname,
                    "Mpo_block_start": qstart_i,
                    "Mpo_block_end": qend_i,
                    "overlap_bp": ov,
                    "R108_chr": tname,
                    "R108_block_start": int(tstart),
                    "R108_block_end": int(tend),
                    "strand": strand,
                    "aligned_bp": alen_i,
                    "identity": round(ident, 4),
                    "mapq": int(mapq),
                }
            )

rows = sorted(rows, key=lambda r: (r["breakpoint_region"], r["Mpo_block_start"], -r["aligned_bp"]))
with OUT.open("w", newline="", encoding="utf-8") as handle:
    fields = [
        "breakpoint_region",
        "Mpo_chr",
        "Mpo_block_start",
        "Mpo_block_end",
        "overlap_bp",
        "R108_chr",
        "R108_block_start",
        "R108_block_end",
        "strand",
        "aligned_bp",
        "identity",
        "mapq",
    ]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

print(OUT)
