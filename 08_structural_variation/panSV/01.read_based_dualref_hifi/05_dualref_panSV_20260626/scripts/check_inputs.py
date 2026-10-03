#!/usr/bin/env python3
import argparse
import os
import sys

CALLERS = {
    "pbsv": ("pbsv", "pbsv"),
    "sniffles2": ("sniffles2", "sniffles2"),
    "cutesv": ("cutesv", "cutesv"),
}


def read_species(path):
    species = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                species.append(line)
    return species


def read_refs(path):
    refs = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            refs.append((parts[idx["ref_alias"]], parts[idx["ref_source"]], parts[idx["ref_fa"]]))
    return refs


def check_file(path, min_size=1):
    if not os.path.exists(path):
        return "missing"
    if os.path.getsize(path) < min_size:
        return "empty"
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    species_path = os.path.join(args.step, "config", "species.list")
    refs_path = os.path.join(args.step, "config", "refs.tsv")
    species = read_species(species_path)
    refs = read_refs(refs_path)

    if "genome_A17" in species:
        sys.stderr.write("ERROR: genome_A17 is present in species.list, but this workflow must exclude it.\n")
        return 2
    if len(species) != 18:
        sys.stderr.write("ERROR: species.list has %d entries; expected 18 after excluding genome_A17.\n" % len(species))
        return 2

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    summary_dir = os.path.join(args.step, "summary")
    os.makedirs(summary_dir, exist_ok=True)
    check_tsv = os.path.join(summary_dir, "input_check.tsv")

    missing = []
    with open(args.out, "w") as manifest, open(check_tsv, "w") as check:
        manifest.write("ref_alias\tref_source\tspecies\tcaller\tvcf\tbam\tbai\tref_fa\n")
        check.write("kind\tref_alias\tref_source\tspecies\tcaller\tpath\tstatus\n")
        for ref_alias, ref_source, ref_fa in refs:
            for path, kind in [(ref_fa, "ref_fa"), (ref_fa + ".fai", "ref_fai")]:
                status = check_file(path)
                check.write("%s\t%s\t%s\t.\t.\t%s\t%s\n" % (kind, ref_alias, ref_source, path, status))
                if status != "ok":
                    missing.append((kind, ref_alias, ".", ".", path, status))
            for sp in species:
                bam = os.path.join(args.run_dir, "03_per_sample", ref_source, "bam", sp + ".sorted.bam")
                bai = bam + ".bai"
                for path, kind in [(bam, "bam"), (bai, "bai")]:
                    status = check_file(path)
                    check.write("%s\t%s\t%s\t%s\t.\t%s\t%s\n" % (kind, ref_alias, ref_source, sp, path, status))
                    if status != "ok":
                        missing.append((kind, ref_alias, sp, ".", path, status))
                for caller, (subdir, suffix) in CALLERS.items():
                    vcf = os.path.join(args.run_dir, "03_per_sample", ref_source, subdir, "%s.%s.vcf.gz" % (sp, suffix))
                    status = check_file(vcf)
                    check.write("vcf\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ref_alias, ref_source, sp, caller, vcf, status))
                    if status != "ok":
                        missing.append(("vcf", ref_alias, sp, caller, vcf, status))
                    manifest.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                        ref_alias, ref_source, sp, caller, vcf, bam, bai, ref_fa))

    if missing:
        sys.stderr.write("ERROR: input check failed; see %s\n" % check_tsv)
        for row in missing[:50]:
            sys.stderr.write("  %s\n" % ("\t".join(row),))
        if len(missing) > 50:
            sys.stderr.write("  ... %d more missing/empty entries\n" % (len(missing) - 50))
        return 2

    done = os.path.join(summary_dir, "input_check.done")
    with open(done, "w") as fh:
        fh.write("ok\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
