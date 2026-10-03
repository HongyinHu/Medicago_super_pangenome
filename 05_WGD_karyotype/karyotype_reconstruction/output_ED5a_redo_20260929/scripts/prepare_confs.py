#!/usr/bin/env python3
"""Stage per-species WGDI runs for the filtered ED5a karyotype mapping.

Re-uses each species' existing WGDI -icl/-ks/-bi outputs under output_Kary/
(parameters identical across species) and writes new configs for
  -c  : ortholog-only blocks (Ks within the speciation peak), >=20 genes
  -km : karyotype mapping from the filtered, genome-wide one-to-one blocks
Mfis uses the new 7-chromosome assembly genome_M46_2 (NOT genome_M46).
"""
import os
import sys

import pandas as pd

ROOT = "path/to/project/37.karyotype_reconstruction"
OUT = os.path.join(ROOT, "output_ED5a_redo_20260929")
BLOCK_LENGTH = 20      # minimum genes per collinear block (-c block_length)
LIMIT_LENGTH = 20      # minimum consecutive genes per mapped segment (-km limit_length)
KS_HIT = 0.5           # fraction of block gene pairs that must fall inside ks_area
# 0.05 (first run, 02_wgdi) removed long, clearly orthologous blocks in gene-sparse
# regions; 0.2 matches the -icl threshold (run 04_wgdi_p0.2).
PVALUE = float(os.environ.get("PVALUE", "0.05"))
RUN = os.environ.get("RUN", "02_wgdi")


def parse_section(path, wanted):
    sec, out = None, {}
    for raw in open(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            sec = line[1:-1].strip()
            continue
        if sec == wanted and "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def resolve(d, name):
    """Same rule as extract_ed5_panel_a_from_wgdi.py: fall back to data/<g>/old
    when the output-dir symlink dangles (Mcar, Mtru_R108 re-annotated in June)."""
    p = name if os.path.isabs(name) else os.path.join(d, name)
    p = os.path.abspath(p)
    if os.path.exists(p):
        return os.path.realpath(p)
    if os.path.islink(p):
        t = os.readlink(p)
        if not os.path.isabs(t):
            t = os.path.normpath(os.path.join(os.path.dirname(p), t))
        old = os.path.join(os.path.dirname(t), "old", os.path.basename(t))
        if os.path.exists(old):
            return os.path.realpath(old)
    raise FileNotFoundError(p)


def main():
    species = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t")
    peaks = pd.read_csv(os.path.join(OUT, "01_ks_survey", "ks_peaks_by_species.tsv"), sep="\t")
    species = species.merge(peaks[["species", "ortholog_peak_ks", "wgd_peak_ks", "trough_ks"]],
                            left_on="label", right_on="species", how="left")
    assert species["trough_ks"].notnull().all()
    assert (species.loc[species.label == "Mfis", "subdir"] == "genome_M46_2").all()
    rows = []
    for _, s in species.iterrows():
        src = os.path.join(ROOT, "output_Kary", s.subdir)
        conf = os.path.join(src, "total.conf")
        bi = parse_section(conf, "blockinfo")
        km = parse_section(conf, "karyotype_mapping")
        files = {
            "gff1": resolve(src, km["gff1"]),
            "gff2": resolve(src, km["gff2"]),
            "lens1": resolve(src, km["the_other_lens"]),
            "lens2": resolve(src, bi["lens2"]),
            "blast": resolve(src, km["blast"]),
            "ancestor": resolve(src, km["ancestor_top"]),
            "blockinfo_raw": resolve(src, bi["savefile"]),
        }
        d = os.path.join(OUT, RUN, s.label)
        os.makedirs(d, exist_ok=True)
        for k, v in files.items():
            link = os.path.join(d, "in." + k)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(v, link)
        ks_cut = float(s.trough_ks)
        with open(os.path.join(d, "c.conf"), "w") as fh:
            fh.write(
                "[correspondence]\n"
                "blockinfo = in.blockinfo_raw\n"
                "lens1 = in.lens1\n"
                "lens2 = in.lens2\n"
                "tandem = true\n"
                "tandem_length = 200\n"
                "pvalue = {p}\n"
                "block_length = {bl}\n"
                "multiple = 1\n"
                "homo = 0,1\n"
                "ks_area = 0,{ks}\n"
                "ks_hit = {kh}\n"
                "savefile = {lab}_aak.ortholog_blocks.csv\n".format(
                    p=PVALUE, bl=BLOCK_LENGTH, ks=ks_cut, kh=KS_HIT, lab=s.label))
        with open(os.path.join(d, "km.conf"), "w") as fh:
            fh.write(
                "[karyotype_mapping]\n"
                "blast = in.blast\n"
                "blast_reverse = false\n"
                "gff1 = in.gff1\n"
                "gff2 = in.gff2\n"
                "score = 100\n"
                "evalue = 1e-5\n"
                "repeat_number = 5\n"
                "ancestor_top = in.ancestor\n"
                "the_other_lens = in.lens1\n"
                "blockinfo = {lab}_aak.ortholog_blocks.1to1.csv\n"
                "blockinfo_reverse = false\n"
                "limit_length = {ll}\n"
                "the_other_ancestor_file = km_result.txt\n".format(lab=s.label, ll=LIMIT_LENGTH))
        with open(os.path.join(d, "k.conf"), "w") as fh:
            fh.write("[karyotype]\nancestor = km_result.txt\nwidth = 0.5\n"
                     "figsize = 10,6.18\nsavefig = karyotype.pdf\n")
        row = {"label": s.label, "subdir": s.subdir, "ks_cut": ks_cut, "pvalue": PVALUE,
               "ortholog_peak_ks": s.ortholog_peak_ks, "wgd_peak_ks": s.wgd_peak_ks}
        row.update(files)
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, RUN, "run_manifest.tsv"), sep="\t", index=False)
    print("staged", len(rows), "species")


if __name__ == "__main__":
    sys.exit(main())
