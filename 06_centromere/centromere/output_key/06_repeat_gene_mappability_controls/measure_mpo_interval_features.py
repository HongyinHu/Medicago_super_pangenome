#!/usr/bin/env python3
import csv
import os

ROOT = "path/to/project/10.centromere_analysis"
RUN_ROOT = os.path.join(ROOT, "output_key")
EDTA = os.path.join(ROOT, "output_key/data/03_repeat/genome_Mpo.EDTA/genome_Mpo.fa.mod.EDTA.TEanno.gff3")
TRASH = os.path.join(ROOT, "output_key/data/03_repeat/genome_Mpo.TRASH/genome_Mpo.TRASH_arrays.sorted.bed")
GENES = os.path.join(ROOT, "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.gene.coord.tsv")
SIGNAL = os.path.join(ROOT, "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.CENH3_vs_Input.q0.50k.log2.bedgraph")
MPO_CEN = os.path.join(ROOT, "output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/genome_Mpo.cen_core.bed")
OUT = os.path.join(RUN_ROOT, "06_repeat_gene_mappability_controls/Mpo_interval_repeat_gene_CENH3_features.tsv")


def ov(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def merge_intervals(items):
    if not items:
        return []
    items = sorted(items)
    merged = [list(items[0])]
    for s, e in items[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def overlap_fraction(chrom, start, end, intervals):
    overlaps = []
    for c, s, e, _label in intervals:
        if c != chrom:
            continue
        x = ov(start, end, s, e)
        if x:
            overlaps.append((max(start, s), min(end, e)))
    merged = merge_intervals(overlaps)
    return sum(e - s for s, e in merged) / max(1, end - start)


def read_edta():
    rows = []
    with open(EDTA) as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom = p[0]
            start = int(p[3]) - 1
            end = int(p[4])
            attr = p[8]
            label = "repeat"
            if "LTR/Gypsy" in attr:
                label = "LTR/Gypsy"
            elif "LTR/Copia" in attr:
                label = "LTR/Copia"
            elif "classification=LTR" in attr:
                label = "LTR"
            rows.append((chrom, start, end, label))
    return rows


def read_bed(path, label):
    rows = []
    with open(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            rows.append((p[0], int(float(p[1])), int(float(p[2])), label))
    return rows


def read_genes():
    rows = []
    with open(GENES) as handle:
        for raw in handle:
            if not raw.strip():
                continue
            p = raw.rstrip("\n").split("\t")
            rows.append((p[1], int(p[2]), int(p[3]), "gene"))
    return rows


def mean_signal(chrom, start, end):
    total = 0.0
    covered = 0
    with open(SIGNAL) as handle:
        for raw in handle:
            if not raw.strip():
                continue
            p = raw.rstrip("\n").split("\t")
            if p[0] != chrom:
                continue
            s, e = int(float(p[1])), int(float(p[2]))
            x = ov(start, end, s, e)
            if x:
                total += float(p[3]) * x
                covered += x
    return total / covered if covered else ""


def intervals_to_measure():
    intervals = [
        ("CEN5_projected_cluster", "Chr4", 18316439, 18432188),
        ("CEN5_projected_cluster_pm500kb", "Chr4", 17816439, 18932188),
        ("CEN5_projected_cluster_pm1Mb", "Chr4", 17316439, 19432188),
    ]
    with open(MPO_CEN) as handle:
        for raw in handle:
            if not raw.strip():
                continue
            p = raw.rstrip("\n").split("\t")
            intervals.append((f"Mpo_active_CEN_{p[0]}", p[0], int(p[1]), int(p[2])))
    return intervals


def main():
    edta = read_edta()
    trash = read_bed(TRASH, "TRASH")
    genes = read_genes()
    all_ltr = [x for x in edta if x[3].startswith("LTR")]
    gypsy = [x for x in edta if x[3] == "LTR/Gypsy"]
    copia = [x for x in edta if x[3] == "LTR/Copia"]
    rows = []
    for name, chrom, start, end in intervals_to_measure():
        size = end - start
        gene_count = sum(1 for c, s, e, _ in genes if c == chrom and ov(start, end, s, e) > 0)
        rows.append({
            "interval_id": name,
            "chr": chrom,
            "start": start,
            "end": end,
            "size_bp": size,
            "CENH3_log2_mean": mean_signal(chrom, start, end),
            "EDTA_repeat_fraction": overlap_fraction(chrom, start, end, edta),
            "LTR_fraction": overlap_fraction(chrom, start, end, all_ltr),
            "Gypsy_fraction": overlap_fraction(chrom, start, end, gypsy),
            "Copia_fraction": overlap_fraction(chrom, start, end, copia),
            "TRASH_fraction": overlap_fraction(chrom, start, end, trash),
            "gene_count": gene_count,
            "gene_density_per_Mb": gene_count / (size / 1e6),
        })
    fields = ["interval_id", "chr", "start", "end", "size_bp", "CENH3_log2_mean", "EDTA_repeat_fraction", "LTR_fraction", "Gypsy_fraction", "Copia_fraction", "TRASH_fraction", "gene_count", "gene_density_per_Mb"]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(OUT)


if __name__ == "__main__":
    main()
