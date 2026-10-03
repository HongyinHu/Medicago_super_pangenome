#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


QUERIES = ("Chr14120.1", "Chr24229.1")


def first_hits(path):
    hits = {}
    if not Path(path).exists():
        return hits
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            hits.setdefault(fields[0], fields)
    return hits


def read_fasta(path):
    records = {}
    header = None
    seq = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = "".join(seq)
                header = line[1:]
                seq = []
            elif header is not None:
                seq.append(line.strip())
    if header is not None:
        records[header.split()[0]] = "".join(seq)
    return records


def write_fasta_record(out, name, seq):
    out.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        out.write(seq[i:i + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--stage", choices=("candidates", "final"), required=True)
    args = parser.parse_args()
    work = Path(args.work)
    species = [x.strip() for x in Path(args.species).read_text().splitlines() if x.strip()]

    if args.stage == "candidates":
        manifest = []
        for sp in species:
            hits = first_hits(work / "blast" / f"{sp}.forward.tsv")
            pep = read_fasta(work / "proteomes" / f"{sp}.pep")
            cds = read_fasta(work / "proteomes" / f"{sp}.cds")
            with open(work / "rbh" / f"{sp}.candidate.pep", "w") as pout, open(work / "rbh" / f"{sp}.candidate.cds", "w") as cout:
                for query in QUERIES:
                    row = hits.get(query)
                    if row is None:
                        manifest.append([sp, query, "", "no_forward_hit"])
                        continue
                    target = row[1]
                    pep_key = target if target in pep else next((k for k in pep if k.startswith(target + ".")), "")
                    cds_key = target if target in cds else next((k for k in cds if k.startswith(target + ".")), "")
                    if not pep_key or not cds_key:
                        manifest.append([sp, query, target, "sequence_missing"])
                        continue
                    label = f"{sp}|{query}|{target}"
                    write_fasta_record(pout, label, pep[pep_key])
                    write_fasta_record(cout, label, cds[cds_key])
                    manifest.append([sp, query, target, "candidate_ready"])
        with open(work / "tables" / "forward_candidates.tsv", "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["species", "msa_query", "candidate_id", "status"])
            writer.writerows(manifest)
        return

    forward_rows = list(csv.DictReader(open(work / "tables" / "forward_candidates.tsv"), delimiter="\t"))
    final_rows = []
    combined_pep = {q: open(work / "rbh" / f"strict_{q}.pep", "w") for q in QUERIES}
    combined_cds = {q: open(work / "rbh" / f"strict_{q}.cds", "w") for q in QUERIES}
    hq_pep = {q: open(work / "rbh" / f"strict_HQ_{q}.pep", "w") for q in QUERIES}
    hq_cds = {q: open(work / "rbh" / f"strict_HQ_{q}.cds", "w") for q in QUERIES}
    try:
        for sp in species:
            reverse = first_hits(work / "blast" / f"{sp}.reverse.tsv")
            pep = read_fasta(work / "rbh" / f"{sp}.candidate.pep")
            cds = read_fasta(work / "rbh" / f"{sp}.candidate.cds")
            for row in [r for r in forward_rows if r["species"] == sp]:
                query = row["msa_query"]
                target = row["candidate_id"]
                forward = first_hits(work / "blast" / f"{sp}.forward.tsv").get(query, [])
                pident = float(forward[2]) if len(forward) > 4 else float("nan")
                qcov = float(forward[4]) if len(forward) > 4 else float("nan")
                label = f"{sp}|{query}|{target}"
                reverse_row = reverse.get(label)
                reverse_best = reverse_row[1] if reverse_row else ""
                strict = reverse_best == query
                status = "strict_RBH" if strict else "non_RBH"
                high_quality = strict and pident >= 90 and qcov >= 90
                final_rows.append([sp, query, target, reverse_best, pident, qcov, int(strict), int(high_quality), status])
                if strict and label in pep and label in cds:
                    write_fasta_record(combined_pep[query], sp, pep[label])
                    write_fasta_record(combined_cds[query], sp, cds[label])
                if high_quality and label in pep and label in cds:
                    write_fasta_record(hq_pep[query], sp, pep[label])
                    write_fasta_record(hq_cds[query], sp, cds[label])
    finally:
        for handle in combined_pep.values():
            handle.close()
        for handle in combined_cds.values():
            handle.close()
        for handle in hq_pep.values():
            handle.close()
        for handle in hq_cds.values():
            handle.close()
    with open(work / "tables" / "strict_RBH_orthologs.tsv", "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["species", "msa_query", "candidate_id", "reverse_best_Msa", "pident", "query_coverage", "strict_RBH", "strict_HQ", "status"])
        writer.writerows(final_rows)


if __name__ == "__main__":
    main()
