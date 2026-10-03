#!/usr/bin/env python3
import argparse
import gzip
import os
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def read_list(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                raise SystemExit("Bad list row: %s" % line)
            rows.append((parts[0], parts[1]))
    return rows


class VcfReader:
    def __init__(self, sample, path):
        self.sample = sample
        self.path = path
        self.fh = open_text(path)
        self.headers = []
        self.first = None
        self._init()

    def _init(self):
        for line in self.fh:
            if line.startswith("#CHROM"):
                continue
            if line.startswith("#"):
                self.headers.append(line)
                continue
            if line.strip():
                self.first = line.rstrip("\n")
                return
        self.first = None

    def next_record(self):
        if self.first is not None:
            rec = self.first
            self.first = None
            return rec
        for line in self.fh:
            if line.strip() and not line.startswith("#"):
                return line.rstrip("\n")
        return None


def sample_dict(format_s, sample_s):
    keys = format_s.split(":") if format_s else []
    vals = sample_s.split(":") if sample_s else []
    return {k: vals[i] if i < len(vals) else "." for i, k in enumerate(keys)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True, help="Two columns: sample and vcf")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = read_list(args.list)
    if not rows:
        raise SystemExit("No VCFs in list")
    readers = [VcfReader(sample, path) for sample, path in rows]
    samples = [x[0] for x in rows]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as out:
        for h in readers[0].headers:
            if h.startswith("##fileformat") or h.startswith("##FORMAT") or not h.startswith("##"):
                out.write(h)
            elif h.startswith("##"):
                out.write(h)
        out.write("##panSV_merge=custom_stream_merge_by_discovery_record_order;FORMAT_fields_are_union_per_record\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(samples))

        n = 0
        while True:
            recs = [r.next_record() for r in readers]
            if all(r is None for r in recs):
                break
            if any(r is None for r in recs):
                raise SystemExit("ERROR: genotype VCF record counts differ near record %d" % (n + 1))
            split = [r.split("\t") for r in recs]
            ids = [s[2] for s in split]
            if len(set(ids)) != 1:
                raise SystemExit("ERROR: genotype VCF order/IDs differ at record %d: %s" %
                                 (n + 1, ",".join(ids[:5])))
            fixed = split[0][:8]
            fmt_order = []
            sample_maps = []
            for s in split:
                fmt = s[8] if len(s) > 8 else "GT"
                smp = s[9] if len(s) > 9 else "./."
                for k in fmt.split(":"):
                    if k and k not in fmt_order:
                        fmt_order.append(k)
                sample_maps.append(sample_dict(fmt, smp))
            if "GT" in fmt_order:
                fmt_order = ["GT"] + [x for x in fmt_order if x != "GT"]
            sample_values = []
            for d in sample_maps:
                sample_values.append(":".join(d.get(k, ".") for k in fmt_order))
            out.write("\t".join(fixed + [":".join(fmt_order)] + sample_values) + "\n")
            n += 1
    sys.stderr.write("merged_genotyped_vcfs\trecords=%d\tout=%s\n" % (n, args.out))


if __name__ == "__main__":
    main()
