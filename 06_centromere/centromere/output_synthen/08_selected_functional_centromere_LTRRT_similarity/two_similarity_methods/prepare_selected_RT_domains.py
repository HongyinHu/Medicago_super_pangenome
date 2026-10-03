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
                    yield name, "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq)


def parse_header(header):
    first = header.split()[0]
    seq_id = first.split("|", 1)[0]
    attrs = {}
    for item in header.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.strip()] = v.strip()
    gene = attrs.get("gene", "")
    clade = attrs.get("clade", "")
    try:
        coverage = float(attrs.get("coverage", "nan"))
    except ValueError:
        coverage = math.nan
    try:
        evalue = float(attrs.get("evalue", "inf"))
    except ValueError:
        evalue = math.inf
    return seq_id, gene, clade, coverage, evalue


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--dom-files", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.metadata, sep="\t")
    meta["seq_id"] = meta["seq_id"].astype(str)
    meta["label"] = meta["label"].astype(str).str.strip()
    prefix_to_selected = {}
    for _, row in meta.iterrows():
        prefix = (
            f"{row['label']}__{row['chr']}__"
            f"LTR{int(row['full_start_1based'])}_{int(row['full_end_1based'])}__"
            f"INT{int(row['internal_start'])}_{int(row['internal_end'])}"
        )
        prefix_to_selected[prefix] = row["seq_id"]

    best = {}
    for dom_file in args.dom_files:
        for header, seq in parse_fasta(dom_file):
            original_seq_id, gene, clade, coverage, evalue = parse_header(header)
            prefix = re.sub(r"__N\d+$", "", original_seq_id)
            if prefix not in prefix_to_selected or gene != "RT":
                continue
            seq_id = prefix_to_selected[prefix]
            score = (0 if math.isnan(coverage) else coverage, -evalue, len(seq))
            if seq_id not in best or score > best[seq_id]["score"]:
                best[seq_id] = {
                    "seq_id": seq_id,
                    "rt_clade": clade,
                    "rt_coverage": coverage,
                    "rt_evalue": evalue,
                    "rt_length_aa": len(seq),
                    "rt_sequence": seq,
                    "score": score,
                }

    rows = []
    fasta_path = outdir / "selected_functional_centromere_RT_domains.fa"
    with open(fasta_path, "w") as out:
        for seq_id in meta["seq_id"]:
            if seq_id not in best:
                continue
            rec = best[seq_id]
            out.write(f">{seq_id}\n")
            seq = rec["rt_sequence"]
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
            rows.append({k: v for k, v in rec.items() if k not in {"rt_sequence", "score"}})

    rt_meta = meta.merge(pd.DataFrame(rows), on="seq_id", how="inner")
    rt_meta.to_csv(outdir / "selected_functional_centromere_RT_domains.metadata.tsv", sep="\t", index=False)
    print(f"selected_loci\t{len(meta)}")
    print(f"with_RT_domain\t{len(rt_meta)}")
    print(rt_meta.groupby(["label", "selected_type"]).size().reset_index(name="count").to_string(index=False))


if __name__ == "__main__":
    main()
