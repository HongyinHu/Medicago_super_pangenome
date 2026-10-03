#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import subprocess
from pathlib import Path


def xopen(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def fasta_records(path):
    name = None
    chunks = []
    with xopen(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).replace("*", "")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).replace("*", "")


def write_fasta(records, path):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def load_fasta_map(path):
    return {name: seq for name, seq in fasta_records(path)}


def run(cmd, log_path):
    with open(log_path, "a") as log:
        log.write(" ".join(map(str, cmd)) + "\n")
    subprocess.run(cmd, check=True)


def load_summary(path, species):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [r for r in rows if r["Species"] == species]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    data = Path(args.data)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for d in ("Msa", "R108"):
        (out / d).mkdir(exist_ok=True)

    mafft = "path/to/home/anaconda3/bin/mafft"
    iqtree2 = "path/to/home/anaconda3/envs/biosofeware/bin/iqtree2"

    target_records = list(fasta_records(args.target))
    summary_path = Path(args.summary)

    genome_pep_maps = {
        "Msa": load_fasta_map(data / "genome_Msa.pep"),
        "R108": load_fasta_map(data / "genome_R108.pep"),
    }

    with open(out / "README.txt", "w") as fh:
        fh.write("Separate SPL trees built from the same Target.SPL.pep.fa reference set.\n")
        fh.write("Each species tree includes all reference target SPL proteins plus only the homologs from that genome.\n")

    for species in ("Msa", "R108"):
        homolog_rows = load_summary(summary_path, species)
        species_dir = out / species
        log_path = species_dir / "commands_log.txt"
        tree_fa = species_dir / f"Target_{species}_SPL_for_tree.fa"
        mafft_fa = species_dir / f"Target_{species}_SPL_for_tree.mafft.fa"
        iq_prefix = f"Target_{species}_SPL_for_tree.iqtree"

        records = []
        seen = set()
        for name, seq in target_records:
            label = f"Ref|{name}"
            if label not in seen:
                seen.add(label)
                records.append((label, seq))
        for r in sorted(homolog_rows, key=lambda x: (x["Best_reference_hit"], x["Gene_ID"])):
            pid = r["Protein_ID"]
            seq = genome_pep_maps[species].get(pid)
            if not seq:
                raise SystemExit(f"Missing protein sequence for {species} protein ID: {pid}")
            label = f"{species}|{r['Gene_ID']}|{pid}|best={r['Best_reference_hit']}"
            if label in seen:
                continue
            seen.add(label)
            records.append((label, seq))

        write_fasta(records, tree_fa)

        with open(log_path, "a") as log:
            log.write(f"MAFFT {tree_fa}\n")
        with open(mafft_fa, "w") as out_mafft, open(species_dir / "mafft.log", "w") as mafft_log:
            subprocess.run([mafft, "--auto", "--thread", str(args.threads), str(tree_fa)], check=True, stdout=out_mafft, stderr=mafft_log)

        # Clean any stale FastTree placeholders from earlier interrupted runs.
        for stale in (species_dir / f"Target_{species}_SPL_for_tree.fasttree.nwk", species_dir / "fasttree.log"):
            if stale.exists():
                stale.unlink()

        with open(log_path, "a") as log:
            log.write(f"IQ-TREE2 {mafft_fa}\n")
        with open(species_dir / "iqtree.stdout.log", "w") as iq_out, open(species_dir / "iqtree.stderr.log", "w") as iq_err:
            subprocess.run([
                iqtree2, "-s", str(mafft_fa), "-m", "LG+G4", "-B", "1000",
                "-T", str(args.threads), "--prefix", iq_prefix
            ], check=True, cwd=str(species_dir), stdout=iq_out, stderr=iq_err)

        species_summary = species_dir / f"Target_{species}_SPL_for_tree.summary.tsv"
        with open(species_summary, "w", newline="") as fh_out:
            if homolog_rows:
                fields = list(homolog_rows[0].keys())
            else:
                fields = []
            if fields:
                w = csv.DictWriter(fh_out, fieldnames=fields, delimiter="\t")
                w.writeheader()
                for r in homolog_rows:
                    w.writerow(r)

        print(f"{species}: homolog genes = {len(homolog_rows)}")
        print(f"{species}: tree FASTA = {tree_fa}")
        print(f"{species}: IQ-TREE prefix = {iq_prefix}")


if __name__ == "__main__":
    main()
