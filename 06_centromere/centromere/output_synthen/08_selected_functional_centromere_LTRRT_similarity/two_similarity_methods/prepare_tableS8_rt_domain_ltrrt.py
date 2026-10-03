#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path

import pandas as pd


def parse_fasta(path):
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name = line[1:]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq).upper()


def header_value(header, key):
    m = re.search(r"(?:^|[ ;]){}=([^; ]+)".format(re.escape(key)), header)
    return m.group(1) if m else ""


def parse_evalue(value):
    try:
        return float(value)
    except Exception:
        return math.inf


def parse_float(value, default=0.0):
    try:
        return float(value)
    except Exception:
        return default


def choose_better(current, candidate):
    if current is None:
        return candidate
    return candidate["sort_key"] < current["sort_key"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--detail", required=True)
    parser.add_argument("--tesorter-dir", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    detail = pd.read_csv(args.detail, sep="\t")
    selected = detail[
        (detail["in_functional_centromere"].astype(str).str.lower() == "true")
        & (detail["has_rt"].astype(str).str.lower() == "true")
    ].copy()
    selected["insertion_time_mya"] = pd.to_numeric(selected["insertion_time_years"], errors="coerce") / 1_000_000.0
    selected_ids = set(selected["seq_id"].astype(str))

    rt_by_locus = {}
    tesorter_dir = Path(args.tesorter_dir)
    for faa in sorted(tesorter_dir.glob("*/*.dom.faa")):
        for header, seq in parse_fasta(faa):
            token = header.split()[0]
            locus = token.split("|", 1)[0]
            if locus not in selected_ids:
                continue
            gene = header_value(header, "gene")
            if gene != "RT" and not token.endswith("-RT"):
                continue
            evalue = parse_evalue(header_value(header, "evalue"))
            coverage = parse_float(header_value(header, "coverage"))
            seq = seq.replace(".", "").replace("-", "")
            candidate = {
                "header": header,
                "seq": seq,
                "evalue": evalue,
                "coverage": coverage,
                "aa_length": len(seq.replace("*", "")),
                "sort_key": (evalue, -coverage, -len(seq.replace("*", "")), header),
            }
            if choose_better(rt_by_locus.get(locus), candidate):
                rt_by_locus[locus] = candidate

    found = [seq_id for seq_id in selected["seq_id"].astype(str) if seq_id in rt_by_locus]
    missing = sorted(selected_ids - set(found))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_path = outdir / "TableS8_functional_centromere_RT_domain.pep.fa"
    meta_path = outdir / "TableS8_functional_centromere_RT_domain.metadata.tsv"
    missing_path = outdir / "TableS8_functional_centromere_RT_domain.missing_ids.txt"

    with open(fasta_path, "w") as out:
        for seq_id in found:
            row = rt_by_locus[seq_id]
            out.write(">{}\n".format(seq_id))
            seq = row["seq"]
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")

    meta = selected.set_index("seq_id").loc[found].reset_index()
    meta["rt_evalue"] = [rt_by_locus[seq_id]["evalue"] for seq_id in found]
    meta["rt_coverage"] = [rt_by_locus[seq_id]["coverage"] for seq_id in found]
    meta["rt_aa_length"] = [rt_by_locus[seq_id]["aa_length"] for seq_id in found]
    meta.to_csv(meta_path, sep="\t", index=False)
    missing_path.write_text("\n".join(missing) + ("\n" if missing else ""))

    print("selected_records={}".format(len(selected)))
    print("rt_domain_records={}".format(len(found)))
    print("missing_rt_domain_records={}".format(len(missing)))
    print(meta.groupby(["label", "superfamily", "clade"]).size().reset_index(name="count").to_string(index=False))
    print("fasta={}".format(fasta_path))
    print("metadata={}".format(meta_path))
    if missing:
        print("missing_ids={}".format(missing_path))


if __name__ == "__main__":
    main()
