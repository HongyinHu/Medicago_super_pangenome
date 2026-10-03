#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import subprocess
from pathlib import Path

INFO_DEFS = [
    '##INFO=<ID=FINAL_SV_ID,Number=1,Type=String,Description="Final species-level panSV identifier">',
    '##INFO=<ID=READ_SV_ID,Number=1,Type=String,Description="Original read-mapping/Jasmine discovery SV identifier">',
    '##INFO=<ID=SOURCE_METHODS,Number=1,Type=String,Description="Evidence source category: read_based or read_based+SVGAP">',
    '##INFO=<ID=SVGAP_SUPPORT_COUNT,Number=1,Type=Integer,Description="Number of overlapping SVGAP/genome-synteny SV records supporting this read-based event">',
    '##INFO=<ID=SVGAP_SUPPORT_SAMPLES,Number=.,Type=String,Description="Species with overlapping SVGAP/genome-synteny evidence after sample-name normalization">',
    '##INFO=<ID=FINAL_CONFIDENCE,Number=1,Type=String,Description="Final confidence label for species-level panSV event">',
]


def clean_info_value(value):
    if value is None or value == '' or value == '.':
        return '.'
    return str(value).replace(' ', '_').replace(';', ',')


def link_or_copy(src, dst):
    src = Path(src)
    dst = Path(dst)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.link(src, dst)
    except OSError:
        # Fallback to symlink to avoid duplicating very large TSVs on shared storage.
        os.symlink(src.resolve(), dst)


def read_events(path):
    by_read = {}
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            read_id = row.get('read_SV_ID')
            if read_id:
                by_read[read_id] = row
    return by_read


def annotate_vcf(in_vcf_gz, events_by_read, out_vcf_gz, bgzip, tabix):
    out_vcf_gz = Path(out_vcf_gz)
    tmp = out_vcf_gz.with_suffix(out_vcf_gz.suffix + '.inprogress')
    if tmp.exists():
        tmp.unlink()
    if out_vcf_gz.exists():
        out_vcf_gz.unlink()
    tbi = Path(str(out_vcf_gz) + '.tbi')
    if tbi.exists():
        tbi.unlink()

    with open(tmp, 'wb') as out_bin:
        proc = subprocess.Popen([bgzip, '-c'], stdin=subprocess.PIPE, stdout=out_bin)
        assert proc.stdin is not None
        with gzip.open(in_vcf_gz, 'rt', encoding='utf-8', errors='replace') as fh:
            inserted = False
            for line in fh:
                if line.startswith('##'):
                    proc.stdin.write(line.encode())
                    continue
                if line.startswith('#CHROM'):
                    if not inserted:
                        for info in INFO_DEFS:
                            proc.stdin.write((info + '\n').encode())
                        inserted = True
                    proc.stdin.write(line.encode())
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 8:
                    continue
                read_id = parts[2]
                ev = events_by_read.get(read_id)
                if ev:
                    final_id = ev['FinalSV_ID']
                    parts[2] = final_id
                    additions = [
                        ('FINAL_SV_ID', final_id),
                        ('READ_SV_ID', read_id),
                        ('SOURCE_METHODS', ev.get('source_methods', '.')),
                        ('SVGAP_SUPPORT_COUNT', ev.get('svgap_support_count', '0')),
                        ('SVGAP_SUPPORT_SAMPLES', ev.get('svgap_support_samples', '.')),
                        ('FINAL_CONFIDENCE', ev.get('final_confidence', '.')),
                    ]
                    extra = ';'.join('%s=%s' % (k, clean_info_value(v)) for k, v in additions)
                    if parts[7] in ('', '.'):
                        parts[7] = extra
                    else:
                        parts[7] = parts[7] + ';' + extra
                proc.stdin.write(('\t'.join(parts) + '\n').encode())
        proc.stdin.close()
        ret = proc.wait()
        if ret != 0:
            raise SystemExit('bgzip failed with exit code %s' % ret)
    tmp.rename(out_vcf_gz)
    subprocess.check_call([tabix, '-f', '-p', 'vcf', str(out_vcf_gz)])


def count_data_lines(path):
    n = 0
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line.startswith('#'):
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pan-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--bgzip', required=True)
    ap.add_argument('--tabix', required=True)
    args = ap.parse_args()

    pan = Path(args.pan_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    discovery = pan / 'panSV.Msa.flankQC.discovery.vcf.gz'
    events = pan / 'panSV.Msa.flankQC.read_primary.events.tsv'
    pav = pan / 'panSV.Msa.flankQC.read_primary.PAV.matrix.tsv'
    evidence = pan / 'panSV.Msa.flankQC.svgap_overlap_evidence.tsv'
    stats = pan / 'panSV.Msa.flankQC.read_primary.stats.tsv'
    read_mapping_stats = pan / 'panSV.Msa.flankQC.read_mapping.stats.tsv'

    for path in [discovery, events, pav, evidence, stats, read_mapping_stats]:
        if not path.exists() or path.stat().st_size == 0:
            raise SystemExit('missing or empty input: %s' % path)

    events_by_read = read_events(events)
    final_vcf = out / 'species_level_panSV.Msa.read_mapping_flankQC.svgap_evidence.vcf.gz'
    annotate_vcf(discovery, events_by_read, final_vcf, args.bgzip, args.tabix)

    link_or_copy(events, out / 'species_level_panSV.Msa.events.tsv')
    link_or_copy(pav, out / 'species_level_panSV.Msa.PAV.matrix.tsv')
    link_or_copy(evidence, out / 'species_level_panSV.Msa.svgap_overlap_evidence.tsv')
    link_or_copy(stats, out / 'species_level_panSV.Msa.stats.tsv')
    link_or_copy(read_mapping_stats, out / 'species_level_panSV.Msa.read_mapping_only.stats.tsv')

    method = out / 'species_level_panSV.Msa.method.tsv'
    with open(method, 'w') as w:
        w.write('item\tvalue\n')
        w.write('reference\tgenome_Msa\n')
        w.write('species_count\t18\n')
        w.write('excluded_samples\tgenome_A17;genome_474b\n')
        w.write('genome_474_policy\tuse_genome_474/read_mapping_and_genome_474a_for_assembly_context\n')
        w.write('read_mapping_callers\tpbsv;sniffles2;cuteSV\n')
        w.write('within_species_merge\tJasmine min_support=2 max_dist=500 nonlinear_dist ignore_strand normalize_type\n')
        w.write('flank_coverage_filter\tleft/right flank 500bp; min_each_flank_depth>=5; mean_flank_depth>=8; mosdepth\n')
        w.write('panSV_merge\tJasmine min_support=1 max_dist=1000 nonlinear_dist ignore_strand normalize_type\n')
        w.write('svgap_policy\tSVGAP/genome-synteny SVs are overlap evidence only; SVGAP-only SVs are not retained\n')
        w.write('svtypes\tDEL;INS;DUP;INV;TRA\n')

    manifest = out / 'MANIFEST.tsv'
    with open(manifest, 'w') as w:
        w.write('file\tdescription\trecords_or_size\n')
        w.write('%s\tFinal annotated VCF; read-mapping+flankQC panSV with SVGAP evidence annotations; no SVGAP-only records\t%d variants\n' % (final_vcf.name, count_data_lines(final_vcf)))
        w.write('species_level_panSV.Msa.events.tsv\tFinal event table with source_methods and SVGAP support fields\t%d events\n' % (sum(1 for _ in open(events)) - 1))
        w.write('species_level_panSV.Msa.PAV.matrix.tsv\tFinal species-level PAV matrix for 18 species\t%d events\n' % (sum(1 for _ in open(pav)) - 1))
        w.write('species_level_panSV.Msa.svgap_overlap_evidence.tsv\tOverlap evidence table linking SVGAP records to read-based events\t%d evidence rows\n' % (sum(1 for _ in open(evidence)) - 1))
        w.write('species_level_panSV.Msa.stats.tsv\tFinal counts and policy summary\t%s\n' % stats.stat().st_size)
        w.write('species_level_panSV.Msa.method.tsv\tMethod/provenance summary\t%s\n' % method.stat().st_size)

    print('FINAL_SPECIES_LEVEL_PANSV_DIR=%s' % out)
    print('FINAL_VCF=%s' % final_vcf)
    print('FINAL_VCF_VARIANTS=%d' % count_data_lines(final_vcf))


if __name__ == '__main__':
    main()
