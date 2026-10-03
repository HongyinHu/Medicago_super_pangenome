#!/usr/bin/env python3
from pathlib import Path
base = Path('path/to/project/N_1.EDTA_single/02_EDTA_single')
out = base/'LAI_recalc_single_EDTA_20260625'

def exists_file(p):
    try:
        return p.exists() and p.is_file() and p.stat().st_size >= 0
    except OSError:
        return False

def list_glob(d, pat):
    try:
        return sorted([p for p in d.glob(pat) if exists_file(p)])
    except OSError:
        return []

def choose_first(paths):
    for p in paths:
        if p and exists_file(p):
            return p
    return None

rows = []
for d in sorted([p for p in base.iterdir() if p.is_dir() and p.name.startswith('genome_')]):
    sample = d.name
    top_fa = []
    for p in d.iterdir():
        if p.name.endswith(('.fa.mod', '.fasta.mod')) and exists_file(p) and '.EDTA.' not in p.name:
            top_fa.append(p)
    top_fa = sorted(top_fa, key=lambda p: (len(p.name), p.name))
    fa = top_fa[0] if top_fa else None
    prefix = fa.name if fa else ''

    pass_candidates = []
    rm_candidates = []
    if prefix:
        raw_ltr = d/(prefix + '.EDTA.raw')/'LTR'
        anno = d/(prefix + '.EDTA.anno')
        pass_candidates.append(raw_ltr/(prefix + '.pass.list'))
        pass_candidates.extend([p for p in list_glob(raw_ltr, '*.pass.list') if 'nmtf.pass.list' not in p.name])
        rm_candidates.append(anno/(prefix + '.out'))
        rm_candidates.extend(list_glob(anno, '*.mod.out'))
        rm_candidates.extend(list_glob(anno, '*.out'))
        # Existing cal_LAI dirs are used as fallback only; EDTA.raw + EDTA.anno above are preferred.
        for cname in ('cal_LAI', 'cal_LAI2'):
            c = d/cname
            pass_candidates.append(c/(prefix + '.pass.list'))
            pass_candidates.extend([p for p in list_glob(c, '*.pass.list') if 'nmtf.pass.list' not in p.name])
            rm_candidates.append(c/(prefix + '.out'))
            rm_candidates.extend([p for p in list_glob(c, '*.out') if not any(x in p.name for x in ('.LAI', '.cat.gz', '.masked', '.tbl'))])
    pass_candidates = [p for p in pass_candidates if p and 'nmtf.pass.list' not in p.name]
    rm_candidates = [p for p in rm_candidates if p and not any(x in p.name for x in ('.LAI', '.cat.gz', '.masked', '.tbl'))]
    pl = choose_first(pass_candidates)
    rm = choose_first(rm_candidates)
    missing = []
    if not fa: missing.append('fa')
    if not pl: missing.append('pass_list')
    if not rm: missing.append('rmout')
    status = 'ready' if not missing else 'missing_' + ','.join(missing)
    rows.append((sample, str(fa or ''), str(pl or ''), str(rm or ''), status))

manifest = out/'manifest_single_EDTA_LAI.tsv'
with manifest.open('w') as f:
    f.write('sample\tfa\tpass_list\trmout\tstatus\n')
    for r in rows:
        f.write('\t'.join(r)+'\n')
print(f'manifest={manifest}')
print(f'total={len(rows)} ready={sum(r[4]=="ready" for r in rows)} not_ready={sum(r[4]!="ready" for r in rows)}')
for r in rows:
    if r[4] != 'ready':
        print('NOT_READY\t' + '\t'.join(r))
