#!/usr/bin/env python3
from __future__ import annotations
import csv, json, os, re, subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

RUN = Path('path/to/project/N_4.pod_spiny/03_candidate_gene_sv_cosegregation_Chr23997')
BASE = Path('path/to/project/N_4.pod_spiny')
N3 = Path('path/to/project/N_3.call_SV/01.read_based_dualref_hifi')
CREATE_REPORT = Path('path/to/home/anaconda3/bin/create_report')
REF_FASTA = BASE / '02_IGV_validation_29SV/references/Msa.candidate_chroms.fa'
GFF = BASE / '00_data/4.reference_anno/genome_Msa.gff'
MATRIX = RUN / 'results/Chr23997.region_sv_genotype_matrix.tsv'
OUT = RUN / 'igv_reports_html'
WORK = OUT / 'report_inputs'
BAM_ROOT = N3 / '03_per_sample/Msa/bam'
GENE_CHROM = 'Chr4'
GENE_START = 88979818
GENE_END = 88982331
GENE_ID = 'Chr23997'

SAMPLES = [
    ('01', 'spiny', 'genome_410', 'main_spiny', 'rgb(205,40,40)'),
    ('02', 'spiny', 'genome_474', 'main_spiny', 'rgb(205,40,40)'),
    ('03', 'spiny', 'genome_Mpo', 'main_spiny', 'rgb(205,40,40)'),
    ('04', 'spiny', 'genome_R108', 'main_spiny', 'rgb(205,40,40)'),
    ('05', 'weak_spiny', 'genome_436', 'ambiguous_weak_spiny', 'rgb(230,145,30)'),
    ('06', 'spiny_lowQ', 'genome_457', 'low_quality_spiny', 'rgb(230,145,30)'),
    ('07', 'spineless', 'genome_395', 'main_spineless', 'rgb(65,110,205)'),
    ('08', 'spineless_hairy', 'genome_454', 'main_spineless_hairy', 'rgb(65,110,205)'),
    ('09', 'spineless', 'genome_461', 'main_spineless', 'rgb(65,110,205)'),
    ('10', 'spineless', 'genome_468', 'main_spineless', 'rgb(65,110,205)'),
    ('11', 'spineless', 'genome_472', 'main_spineless', 'rgb(65,110,205)'),
    ('12', 'spineless', 'genome_482', 'main_spineless', 'rgb(65,110,205)'),
    ('13', 'spineless', 'genome_M22', 'main_spineless', 'rgb(65,110,205)'),
    ('14', 'spineless', 'genome_Mar', 'main_spineless', 'rgb(65,110,205)'),
    ('15', 'spineless', 'genome_Mru', 'main_spineless', 'rgb(65,110,205)'),
    ('16', 'spineless_ref', 'genome_Msa', 'main_spineless_ref', 'rgb(35,145,90)'),
    ('17', 'spineless', 'genome_ZM4', 'main_spineless', 'rgb(65,110,205)'),
    ('18', 'spineless_lowQ', 'genome_M46', 'low_quality_spineless', 'rgb(120,120,120)'),
]

TARGET_IDS = {
    '1_0_pbsv.INS.137428': 'GB01_gene_body_intronic_INS',
    '0_0_pbsv.DEL.130205': 'GB02_gene_body_intronic_DEL',
    '2_0_pbsv.DEL.106554': 'GB03_gene_body_intronic_DEL',
    '2_0_pbsv.DEL.106551': 'GB04_gene_body_intronic_DEL',
    '3_0_pbsv.DEL.106472': 'GB05_gene_body_CDS_DEL',
    '10_0_pbsv.DEL.136629': 'SP01_downstream_spiny_enriched_DEL',
    '1_0_pbsv.INS.137368': 'SP02_downstream_spiny_enriched_INS',
    '1_0_pbsv.DEL.137490': 'SP03_upstream_spiny_enriched_DEL',
    '1_0_pbsv.INS.137487': 'SP04_upstream_spiny_enriched_INS',
    '3_1_Sniffles2.DEL.5949S3': 'SL01_upstream_spineless_enriched_DEL',
    '2_1_Sniffles2.INS.22F9S3': 'SL02_downstream_spineless_enriched_INS',
}


def safe(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s).strip('_')


def read_fai(path: Path):
    d = {}
    with path.open() as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                d[p[0]] = int(p[1])
    return d


def write_vcf(path: Path, contigs: dict[str,int], r: dict[str,str]):
    chrom, pos, end = r['chrom'], int(r['pos']), int(r['end'])
    svtype, svlen = r['svtype'], r['svlen']
    alt = f'<{svtype}>'
    with path.open('w') as out:
        out.write('##fileformat=VCFv4.2\n')
        for c,l in contigs.items(): out.write(f'##contig=<ID={c},length={l}>\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        vid = safe(f"{TARGET_IDS[r['sv_id']]}|{r['sv_id']}|{svtype}|{chrom}:{pos}-{end}")
        out.write(f'{chrom}\t{pos}\t{vid}\tN\t{alt}\t.\tPASS\tEND={end};SVTYPE={svtype};SVLEN={svlen}\n')


def write_marker(path: Path, r: dict[str,str]):
    chrom, pos, end = r['chrom'], int(r['pos']), int(r['end'])
    left, right = min(pos, end), max(pos, end)
    if right <= left: right = left + 1
    svtype, svlen = r['svtype'], r['svlen']
    with path.open('w') as out:
        out.write(f'{chrom}\t{max(left-1,0)}\t{right}\t{svtype}_len{svlen}_Area\t0\t.\t{max(left-1,0)}\t{right}\t170,170,170\n')
        out.write(f'{chrom}\t{max(pos-1,0)}\t{max(pos,1)}\tSTART\t1000\t.\t{max(pos-1,0)}\t{max(pos,1)}\t255,0,0\n')
        if abs(end-pos) > 5:
            out.write(f'{chrom}\t{max(end-1,0)}\t{max(end,1)}\tEND\t1000\t.\t{max(end-1,0)}\t{max(end,1)}\t0,0,255\n')
        # gene body marker in each report, useful when nearby SV is outside gene
        out.write(f'{GENE_CHROM}\t{GENE_START-1}\t{GENE_END}\t{GENE_ID}_gene_body\t500\t-\t{GENE_START-1}\t{GENE_END}\t80,80,80\n')


def write_gff_subset(path: Path, chrom: str, start: int, end: int):
    start0 = max(1, start)
    end0 = end
    with GFF.open() as inp, path.open('w') as out:
        out.write('##gff-version 3\n')
        for line in inp:
            if not line.strip() or line.startswith('#'): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[0] != chrom: continue
            fs, fe = int(p[3]), int(p[4])
            if min(fe, end0) >= max(fs, start0):
                out.write(line)


def link_bams(cand_dir: Path):
    cand_dir.mkdir(parents=True, exist_ok=True)
    for order, phen, sample, group, color in SAMPLES:
        src = BAM_ROOT / f'{sample}.sorted.bam'
        bai = Path(str(src) + '.bai')
        name = f'{order}_{phen}_{sample}'
        dst = cand_dir / f'{name}.bam'
        dbi = cand_dir / f'{name}.bam.bai'
        if not dst.exists(): dst.symlink_to(src)
        if not dbi.exists(): dbi.symlink_to(bai)


def write_tracks(path: Path, cand_dir: Path, gff_subset: Path):
    tracks = [{
        'name': 'Msa gene model', 'url': str(gff_subset), 'format': 'gff3', 'type': 'annotation',
        'displayMode': 'EXPANDED', 'height': 150, 'color': 'rgb(40,40,40)'
    }]
    for order, phen, sample, group, color in SAMPLES:
        name = f'{order}_{phen}_{sample}'
        tracks.append({
            'name': f'{name} ({group})', 'url': str(cand_dir / f'{name}.bam'),
            'indexURL': str(cand_dir / f'{name}.bam.bai'), 'format': 'bam', 'type': 'alignment',
            'height': 90, 'color': color, 'showCoverage': True
        })
    with path.open('w') as out: json.dump(tracks, out, indent=2)


def html_escape(s):
    return str(s).replace('&','&amp;').replace('<','&lt;').replace('>','&gt;').replace('"','&quot;')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    rows=[]
    with MATRIX.open() as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['sv_id'] in TARGET_IDS:
                rows.append(r)
    rows.sort(key=lambda r: TARGET_IDS[r['sv_id']])
    contigs = read_fai(Path(str(REF_FASTA)+'.fai'))
    commands=[]; index=[]
    manifest = OUT / 'Chr23997_igv_reports_manifest.tsv'
    with manifest.open('w') as man:
        man.write('label\tsv_id\tlocus\tsvtype\tsvlen\tcontext\tstatus\tspiny_present\tspineless_present\thtml\n')
        for r in rows:
            label = TARGET_IDS[r['sv_id']]
            tag = safe(f"{label}_{r['sv_id']}_{r['svtype']}_{r['chrom']}_{r['pos']}_{r['end']}")
            w = WORK / tag
            w.mkdir(parents=True, exist_ok=True)
            cand_dir = w / 'bams'
            link_bams(cand_dir)
            vcf = w / f'{tag}.vcf'
            marker = w / f'{tag}.mark.bed'
            gff_subset = w / f'{tag}.gene_model.gff3'
            track_json = w / f'{tag}.tracks.json'
            write_vcf(vcf, contigs, r)
            write_marker(marker, r)
            pos, end = int(r['pos']), int(r['end'])
            # detailed breakpoint view; keep enough flank for local reads, not the whole +/-100 kb unless the SV itself is long
            flank = 5000 if abs(end-pos) < 5000 else 10000
            write_gff_subset(gff_subset, r['chrom'], min(pos,end)-flank, max(pos,end)+flank)
            write_tracks(track_json, cand_dir, gff_subset)
            html = OUT / f'{tag}.html'
            title = f"{label} {r['sv_id']} {r['svtype']} len={r['svlen']} {r['chrom']}:{r['pos']}-{r['end']} {r['status']} {r['context']}"
            cmd = [str(CREATE_REPORT), str(vcf), '--fasta', str(REF_FASTA), '--flanking', str(flank), '--tracks', str(marker), '--roi', str(marker), '--track-config', str(track_json), '--standalone', '--title', title, '--output', str(html)]
            commands.append(cmd)
            locus=f"{r['chrom']}:{r['pos']}-{r['end']}"
            man.write('\t'.join([label,r['sv_id'],locus,r['svtype'],r['svlen'],r['context'],r['status'],r['spiny_present'],r['spineless_present'],str(html)])+'\n')
            index.append((label,r,html,flank))
    cmdlist = OUT / 'Chr23997_igv_report_commands.list'
    with cmdlist.open('w') as out:
        for cmd in commands:
            out.write(' '.join(subprocess.list2cmdline([x]) for x in cmd)+'\n')
    failures=[]
    def run(cmd):
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return cmd,res.returncode,res.stdout,res.stderr
    with ThreadPoolExecutor(max_workers=2) as ex:
        futs=[ex.submit(run,cmd) for cmd in commands]
        for fut in as_completed(futs):
            cmd,code,so,se=fut.result()
            log = OUT / (Path(cmd[-1]).stem + '.create_report.log')
            log.write_text(so + '\n--- STDERR ---\n' + se)
            if code != 0:
                failures.append((cmd,code,log))
            print(('OK' if code==0 else f'FAIL={code}'), cmd[-1], flush=True)
    # index after reports
    index_path = OUT / 'index.html'
    with index_path.open('w') as out:
        out.write("<!doctype html><html><head><meta charset='utf-8'><title>Chr23997 IGV SV reports</title>")
        out.write("<style>body{font-family:Arial,sans-serif;margin:24px}table{border-collapse:collapse;font-size:13px}th,td{border:1px solid #ddd;padding:6px 8px;vertical-align:top}th{background:#f3f5f7}.sp{background:#fff4e6}.gb{background:#eef7ff}.sl{background:#f0f0f0}</style></head><body>")
        out.write('<h1>Chr23997 SV IGV reports</h1>')
        out.write('<p>Tracks: 4 main spiny, 2 weak/low-quality spiny flags, 11 main spineless, 1 low-quality spineless flag. Red/orange tracks are spiny-side materials; blue/green/gray are spineless-side materials.</p>')
        out.write(f'<p>Gene: {GENE_ID}, Msa {GENE_CHROM}:{GENE_START}-{GENE_END}(-). Reports are local breakpoint views with gene model if overlapping the view.</p>')
        out.write('<table><thead><tr><th>report</th><th>label</th><th>SV ID</th><th>locus</th><th>type</th><th>len</th><th>context</th><th>status</th><th>main spiny present</th><th>main spineless present</th><th>flank bp</th></tr></thead><tbody>')
        for label,r,html,flank in index:
            cls = 'gb' if label.startswith('GB') else 'sp' if label.startswith('SP') else 'sl'
            rel=os.path.relpath(html, OUT)
            out.write(f"<tr class='{cls}'><td><a href='{html_escape(rel)}'>open</a></td><td>{html_escape(label)}</td><td>{html_escape(r['sv_id'])}</td><td>{html_escape(r['chrom']+':'+r['pos']+'-'+r['end'])}</td><td>{html_escape(r['svtype'])}</td><td>{html_escape(r['svlen'])}</td><td>{html_escape(r['context'])}</td><td>{html_escape(r['status'])}</td><td>{html_escape(r['spiny_present'])}/4</td><td>{html_escape(r['spineless_present'])}/11</td><td>{flank}</td></tr>")
        out.write('</tbody></table></body></html>')
    if failures:
        for cmd,code,log in failures:
            print(f'FAIL {code}: {cmd[-1]} log={log}', flush=True)
        raise SystemExit(1)
    print(f'Generated {len(index)} IGV reports')
    print(f'Index: {index_path}')
    print(f'Manifest: {manifest}')

if __name__ == '__main__':
    main()
