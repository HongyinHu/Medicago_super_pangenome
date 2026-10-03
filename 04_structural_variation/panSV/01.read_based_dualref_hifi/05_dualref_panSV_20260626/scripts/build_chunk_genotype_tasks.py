#!/usr/bin/env python3
import argparse
from pathlib import Path


def read_species(path: Path):
    return [x.strip() for x in path.read_text().splitlines() if x.strip() and not x.startswith("#")]


def read_refs(path: Path):
    refs = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            refs.append((row[idx["ref_alias"]], row[idx["ref_source"]], row[idx["ref_fa"]]))
    return refs


def read_chunks(path: Path):
    chunks = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            chunks.append(
                (
                    row[idx["chunk_id"]],
                    row[idx["vcf"]],
                    row[idx["bed"]],
                )
            )
    return chunks


def main():
    ap = argparse.ArgumentParser(description="Build interleaved Slurm-array task list for chunked Sniffles2 genotype.")
    ap.add_argument("--step", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    step = Path(args.step)
    species = read_species(step / "config" / "species.list")
    refs = read_refs(step / "config" / "refs.tsv")
    chunks_by_ref = {
        ref_alias: read_chunks(step / "results" / ref_alias / "genotype_chunks" / "chunks.tsv")
        for ref_alias, _, _ in refs
    }
    max_chunks = max(len(v) for v in chunks_by_ref.values())

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with out.open("w", newline="\n") as oh:
        for sp in species:
            for i in range(max_chunks):
                for ref_alias, ref_source, ref_fa in refs:
                    chunks = chunks_by_ref[ref_alias]
                    if i >= len(chunks):
                        continue
                    chunk_id, chunk_vcf, chunk_bed = chunks[i]
                    out_dir = step / "results" / ref_alias / "genotype_chunks" / sp
                    out_vcf = out_dir / f"{chunk_id}.vcf"
                    done = out_dir / f"{chunk_id}.done"
                    if out_vcf.is_file() and out_vcf.stat().st_size > 0 and done.is_file() and done.stat().st_size >= 0:
                        continue
                    oh.write("\t".join([ref_alias, ref_source, ref_fa, sp, chunk_id, chunk_vcf, chunk_bed]) + "\n")
                    written += 1
    print(f"tasks={written}")


if __name__ == "__main__":
    main()
