#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare WGDI-like inputs:
1) PREFIX.gff  : chr_num  gene_id  start  end  strand  order  genename
2) PREFIX.lens : chr_num  chr_length  gene_count_on_chr
3) PREFIX.cds  : renamed CDS fasta
4) PREFIX.pep  : renamed PEP fasta (remove '.')

Key idea:
- Robustly map contig/chrom names to chr numbers (1..8) using regex,
  instead of exact string match.
- Use chr_num ("1".."8") as unified chromosome key across fasta and gff.

Usage:
  python wgdi_prepare.py genome.fa genome.gff3 genome.cds.fa genome.pep.fa PREFIX

Notes:
- It assumes genes are on chromosomes 1..8. Others (contigs/scaffolds) are skipped.
- For gene name, we try multiple fields; if none found, write "NA".
"""

import re
import sys
from collections import defaultdict
from Bio import SeqIO


# ----------------------------
# helpers
# ----------------------------
def parse_attr(attr: str) -> dict:
    """Parse GFF3 attributes into dict; tolerant to missing/odd fields."""
    d = {}
    if not attr:
        return d
    for part in attr.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def get_chr_num(name: str, min_chr=1, max_chr=8):
    """
    Extract chromosome number from seqid.
    Accepts patterns like:
      Chr1, chr1, MtrunA17Chr1, ... Chr01 ...
    Returns '1'..'8' as string, or None if not matched / out of range.
    """
    if not name:
        return None
    # FASTA header sometimes has spaces after ID, keep the first token
    token = name.strip().split()[0]

    # Common patterns: ...Chr1... / ...chr1... / ...MtrunA17Chr1...
    m = re.search(r"(?:Chr|chr|GWHBEDJ|MtrunA17Chr|LG)(0*\d+)", token)
    if not m:
        return None
    n = int(m.group(1))
    if min_chr <= n <= max_chr:
        return str(n)
    return None


def guess_gene_name(attrs: dict) -> str:
    """
    Try common fields for a human-readable gene name.
    Priority can be adjusted.
    """
    for key in ("gene", "gene_name", "Name", "Parent", "ID"):
        if key in attrs and attrs[key]:
            return attrs[key].split(",")[0]
    return "NA"


# ----------------------------
# main
# ----------------------------
def main(genome_fa, genome_gff, genome_cds, genome_pep, prefix, min_chr=1, max_chr=8):
    # 1) chr lengths from fasta
    chr_len = {}
    for rec in SeqIO.parse(genome_fa, "fasta"):
        chr_num = get_chr_num(rec.id, min_chr=min_chr, max_chr=max_chr)
        if chr_num is not None:
            chr_len[chr_num] = len(rec.seq)

    if not chr_len:
        sys.stderr.write(
            "[ERROR] No chromosome lengths were captured from FASTA.\n"
            "        Check FASTA headers contain Chr1..Chr8 (or chr/MtrunA17Chr).\n"
        )
        sys.exit(1)

    # 2) parse GFF, write PREFIX.gff, build rename map for CDS/PEP
    # gene_name_change: old mRNA/transcript ID -> new gene ID
    gene_name_change = {}

    # track gene order per chr (1..8)
    gene_order_counter = defaultdict(int)

    # gene count per chr for lens
    gene_count = defaultdict(int)

    out_gff = f"{prefix}.gff"
    out_lens = f"{prefix}.lens"
    out_cds = f"{prefix}.cds"
    out_pep = f"{prefix}.pep"

    with open(genome_gff, "r") as inf, open(out_gff, "w") as ouf:
        for line in inf:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attr = cols

            chr_num = get_chr_num(seqid, min_chr=min_chr, max_chr=max_chr)
            if chr_num is None:
                continue

            # 只处理mRNA（与你原脚本一致）；如果你想支持transcript，也可加条件
            if ftype != "mRNA":
                continue

            attrs = parse_attr(attr)
            old_id = attrs.get("ID")
            if not old_id:
                continue

            gene_order_counter[chr_num] += 1
            order = gene_order_counter[chr_num]

            new_id = f"{prefix}_{int(chr_num)}g{order:05d}"
            gene_name_change[old_id] = new_id

            genename = guess_gene_name(attrs)

            # 输出列保持你原始结构：
            # chr_num, new_id, start, end, strand, order, genename
            ouf.write("\t".join(map(str, [chr_num, new_id, start, end, strand, order, genename])) + "\n")

            gene_count[chr_num] += 1

    # 3) write lens (always sorted by chr num)
    with open(out_lens, "w") as ouf:
        for chr_num in sorted(chr_len.keys(), key=lambda x: int(x)):
            ouf.write("\t".join(map(str, [chr_num, chr_len[chr_num], gene_count.get(chr_num, 0)])) + "\n")

    # 4) rename CDS
    with open(out_cds, "w") as ouf:
        for rec in SeqIO.parse(genome_cds, "fasta"):
            if rec.id in gene_name_change:
                ouf.write(f">{gene_name_change[rec.id]}\n{str(rec.seq)}\n")

    # 5) rename PEP (remove '.' as in your script)
    with open(out_pep, "w") as ouf:
        for rec in SeqIO.parse(genome_pep, "fasta"):
            if rec.id in gene_name_change:
                seq = str(rec.seq).replace(".", "")
                ouf.write(f">{gene_name_change[rec.id]}\n{seq}\n")

    sys.stderr.write(
        "[DONE]\n"
        f"  {out_gff}\n"
        f"  {out_lens}\n"
        f"  {out_cds}\n"
        f"  {out_pep}\n"
        f"  Chromosomes captured from FASTA: {', '.join(sorted(chr_len.keys(), key=lambda x: int(x)))}\n"
    )


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.stderr.write(
            "Usage:\n"
            f"  python {sys.argv[0]} genome_fa genome_gff genome_cds genome_pep prefix\n"
        )
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
