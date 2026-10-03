#窗口甲基化统计脚本

#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

parser = argparse.ArgumentParser()
parser.add_argument("--cx", required=True, help="Bismark CX_report.txt or .gz")
parser.add_argument("--chrom_sizes", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--window", type=int, default=100000)
parser.add_argument("--min_cov", type=int, default=5)
parser.add_argument("--context", default="ALL", choices=["ALL", "CG", "CHG", "CHH"])
args = parser.parse_args()

chrom_sizes = {}
with open(args.chrom_sizes) as f:
    for line in f:
        chrom, size = line.strip().split()[:2]
        chrom_sizes[chrom] = int(size)

meth = defaultdict(int)
total = defaultdict(int)

with open_maybe_gz(args.cx) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        arr = line.rstrip("\n").split()
        if len(arr) < 6:
            continue

        chrom = arr[0]
        if chrom not in chrom_sizes:
            continue

        pos = int(arr[1]) - 1
        m = int(arr[3])
        u = int(arr[4])
        ctx = arr[5]

        cov = m + u
        if cov < args.min_cov:
            continue

        if args.context != "ALL" and ctx != args.context:
            continue

        bin_start = (pos // args.window) * args.window
        key = (chrom, bin_start)

        meth[key] += m
        total[key] += cov

with open(args.out, "w") as out:
    for chrom in chrom_sizes:
        size = chrom_sizes[chrom]
        for start in range(0, size, args.window):
            end = min(start + args.window, size)
            key = (chrom, start)

            if total[key] == 0:
                continue

            value = meth[key] / total[key]
            out.write(f"{chrom}\t{start}\t{end}\t{value:.6f}\n")