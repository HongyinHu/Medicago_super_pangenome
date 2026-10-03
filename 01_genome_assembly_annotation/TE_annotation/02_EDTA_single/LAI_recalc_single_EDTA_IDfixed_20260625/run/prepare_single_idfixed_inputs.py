#!/usr/bin/env python3
from pathlib import Path
import csv, sys
BASE=Path('path/to/project/N_1.EDTA_single')
IN=BASE/'02_EDTA_single'
OUT=IN/'LAI_recalc_single_EDTA_IDfixed_20260625'

def exists(p):
    try: return p.exists() and p.stat().st_size > 0
    except OSError: return False

def headers(path):
    hs=[]
    with path.open('r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('>'):
                hs.append(line[1:].strip().split()[0])
    return hs

def choose_mod(d):
    c=[]
    for p in d.iterdir():
        if '.EDTA.' in p.name: continue
        if p.name.endswith(('.fa.mod','.fasta.mod')) and exists(p):
            c.append(p)
    return sorted(c, key=lambda p:(len(p.name),p.name))[0] if c else None

def first_non_nmtf(globbed):
    arr=sorted([p for p in globbed if exists(p) and 'nmtf.pass.list' not in p.name])
    return arr[0] if arr else None

def rewrite_out(rmout, out_path, mapping):
    total=mapped=unchanged=unmapped=0
    unmapped_ids={}
    with rmout.open('r', encoding='utf-8', errors='replace') as inp, out_path.open('w', encoding='utf-8') as out:
        for line in inp:
            stripped=line.strip()
            if not stripped:
                out.write(line); continue
            parts=stripped.split()
            if len(parts) >= 5 and parts[0].lstrip('-').isdigit():
                total += 1
                q=parts[4]
                if q in mapping:
                    parts[4]=mapping[q]
                    mapped += 1
                else:
                    unchanged += 1
                    # record suspicious IDs that are not in the original header set after rewrite
                    if q.startswith('_J'):
                        unmapped += 1
                        unmapped_ids[q]=unmapped_ids.get(q,0)+1
                out.write('\t'.join(parts)+'\n')
            else:
                out.write(line)
    return total,mapped,unchanged,unmapped,','.join(list(unmapped_ids)[:20])

rows=[]
queue=[]
for d in sorted([p for p in IN.iterdir() if p.is_dir() and p.name.startswith('genome_')]):
    sample=d.name
    work=OUT/sample
    work.mkdir(parents=True, exist_ok=True)
    mod=choose_mod(d)
    prefix=mod.name if mod else ''
    orig=d/prefix[:-4] if prefix.endswith('.mod') else None
    if not (orig and exists(orig)):
        # fallback to sample.fa only if the EDTA original name is unavailable
        fa_candidates=[p for p in d.iterdir() if p.name.endswith(('.fa','.fasta')) and exists(p) and '.mod' not in p.name]
        orig=sorted(fa_candidates, key=lambda p:(p.name!=f'{sample}.fa', len(p.name), p.name))[0] if fa_candidates else None
    raw_ltr=d/(prefix+'.EDTA.raw')/'LTR' if prefix else Path('')
    anno=d/(prefix+'.EDTA.anno') if prefix else Path('')
    pass_list=raw_ltr/(prefix+'.pass.list') if prefix else None
    if not (pass_list and exists(pass_list)):
        pass_list=first_non_nmtf(raw_ltr.glob('*.pass.list')) if prefix and raw_ltr.exists() else None
    if not (pass_list and exists(pass_list)):
        for cname in ('cal_LAI','cal_LAI2'):
            c=d/cname
            pass_list=first_non_nmtf(c.glob('*.pass.list')) if c.exists() else None
            if pass_list: break
    rmout=anno/(prefix+'.out') if prefix else None
    if not (rmout and exists(rmout)):
        cand=sorted([p for p in anno.glob('*.mod.out') if exists(p)]) if prefix and anno.exists() else []
        rmout=cand[0] if cand else None
    if not (rmout and exists(rmout)):
        for cname in ('cal_LAI','cal_LAI2'):
            c=d/cname
            cand=sorted([p for p in c.glob('*.out') if exists(p) and not any(x in p.name for x in ('.LAI','.cat.gz','.masked','.tbl'))]) if c.exists() else []
            if cand:
                rmout=cand[0]; break
    missing=[]
    if not (orig and exists(orig)): missing.append('original_genome')
    if not (mod and exists(mod)): missing.append('mod_genome')
    if not (pass_list and exists(pass_list)): missing.append('pass_list')
    if not (rmout and exists(rmout)): missing.append('rmout')
    idfixed=work/f'{sample}.singleEDTA.IDfixed.out'
    genome_link=work/f'{sample}.singleEDTA.genome.fa'
    pass_link=work/f'{sample}.singleEDTA.pass.list'
    status='ok'
    seq_count=mapped=unchanged=unmapped=total=0
    unmapped_ids=''
    note=''
    if missing:
        status='missing_' + ','.join(missing)
        note=status
    else:
        orig_h=headers(orig); mod_h=headers(mod)
        if not orig_h or not mod_h or len(orig_h)!=len(mod_h):
            status='bad_header_count'
            note=f'orig_headers={len(orig_h)} mod_headers={len(mod_h)}'
        else:
            mapping=dict(zip(mod_h, orig_h))
            total,mapped,unchanged,unmapped,unmapped_ids=rewrite_out(rmout, idfixed, mapping)
            seq_count=len(mapping)
            genome_link.unlink(missing_ok=True); genome_link.symlink_to(orig)
            pass_link.unlink(missing_ok=True); pass_link.symlink_to(pass_list)
            status='ok' if unmapped==0 else 'ok_with_unmapped_mod_ids'
            note=f'total_records={total}'
            queue.append(sample)
    rows.append(dict(sample=sample, original_genome=str(orig or ''), mod_genome=str(mod or ''), pass_list=str(pass_list or ''), rmout=str(rmout or ''), workdir=str(work), genome_link=str(genome_link), pass_link=str(pass_link), idfixed_out=str(idfixed), seq_count=seq_count, mapped_records=mapped, unchanged_records=unchanged, unmapped_records=unmapped, status=status, note=note, unmapped_ids=unmapped_ids))

manifest=OUT/'prepare_manifest_single_idfixed.tsv'
fields=['sample','original_genome','mod_genome','pass_list','rmout','workdir','genome_link','pass_link','idfixed_out','seq_count','mapped_records','unchanged_records','unmapped_records','status','note','unmapped_ids']
with manifest.open('w', newline='') as f:
    w=csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)
with (OUT/'sample.queue').open('w') as f:
    for s in queue:
        f.write(s+'\n')
print(f'manifest={manifest}')
print(f'total={len(rows)} ok={sum(r["status"].startswith("ok") for r in rows)} bad={sum(not r["status"].startswith("ok") for r in rows)} queue={len(queue)}')
for r in rows:
    if not r['status'].startswith('ok'):
        print('BAD\t' + '\t'.join(str(r[k]) for k in fields))
