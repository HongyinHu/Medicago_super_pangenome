#!/usr/bin/env bash
set -euo pipefail
BASE=${BASE:-path/to/project/N_3.call_SV}
READ_RUN=${READ_RUN:-$BASE/01.read_based_dualref_hifi}
OUT=${OUT:-$BASE/02.assembly_svgap_dualref}
SRC=${SRC:-$BASE/00_data/1.reference_genome}
PYTHON=${PYTHON:-path/to/home/anaconda3/envs/biosofeware/bin/python}
SVGAP=${SVGAP:-$OUT/software/SVGAP}
mkdir -p "$OUT" "$OUT/01_metadata" "$OUT/genome" "$OUT/logs" "$OUT/logs/tasks" "$OUT/logs/workers" "$OUT/status/wga" "$OUT/locks/wga" "$OUT/tmp"
mkdir -p "$OUT/Msa_ref/alignment" "$OUT/Msa_ref/chainnet" "$OUT/R108_ref/alignment" "$OUT/R108_ref/chainnet"
"$PYTHON" - <<'PY'
import csv, gzip, os, re, sys
base=os.environ.get('BASE','path/to/project/N_3.call_SV')
read_run=os.environ.get('READ_RUN',base+'/01.read_based_dualref_hifi')
out=os.environ.get('OUT',base+'/02.assembly_svgap_dualref')
src=os.environ.get('SRC',base+'/00_data/1.reference_genome')
meta=os.path.join(read_run,'01_metadata','assembly_samples.tsv')
if not os.path.exists(meta):
    raise SystemExit(f'missing metadata: {meta}')

def internal_id(asm):
    if asm == 'genome_Msa': return 'Msa'
    if asm == 'genome_R108': return 'R108'
    x = asm
    if x.startswith('genome_'):
        x = x[len('genome_'):]
    x = re.sub(r'[^A-Za-z0-9]+','',x)
    if not x:
        x = re.sub(r'[^A-Za-z0-9]+','',asm)
    if x and x[0].isdigit():
        x = 'G' + x
    return x

def clean_chr(h):
    h = h.strip().split()[0]
    h = re.sub(r'^>','',h)
    h = re.sub(r'[^A-Za-z0-9]+','x',h)
    return h or 'ctg'

def opener(path):
    return gzip.open(path,'rt') if path.endswith('.gz') else open(path,'r')

rows=[]
with open(meta, newline='') as fh:
    reader=csv.DictReader(fh, delimiter='\t')
    for r in reader:
        asm=r['assembly_id']
        iid=internal_id(asm)
        fasta=os.path.join(src, asm + '.fa')
        if not os.path.exists(fasta):
            for ext in ('.fasta','.fna','.fa.gz','.fasta.gz','.fna.gz'):
                cand=os.path.join(src, asm + ext)
                if os.path.exists(cand):
                    fasta=cand; break
        if not os.path.exists(fasta):
            raise SystemExit(f'missing fasta for {asm}: {fasta}')
        r=dict(r)
        r['internal_id']=iid
        r['input_fasta']=fasta
        r['prepared_fasta']=os.path.join(out,'genome',iid)
        rows.append(r)
ids=[r['internal_id'] for r in rows]
if len(ids) != len(set(ids)):
    dup=[x for x in ids if ids.count(x)>1]
    raise SystemExit('duplicate internal IDs: '+','.join(sorted(set(dup))))
map_path=os.path.join(out,'01_metadata','assembly_id_map.tsv')
header_path=os.path.join(out,'01_metadata','header_map.tsv')
with open(map_path,'w',newline='') as mf, open(header_path,'w',newline='') as hf:
    fields=['assembly_id','internal_id','source_group','role','pod_coiling','pod_spine','assembly_gb','independent_trait_sample','input_fasta','prepared_fasta']
    mf.write('\t'.join(fields)+'\n')
    hf.write('assembly_id\tinternal_id\toriginal_header\tsvgap_header\tlength\n')
    for r in rows:
        mf.write('\t'.join(str(r.get(k,'')) for k in fields)+'\n')
        outfa=r['prepared_fasta']
        tmp=outfa+'.tmp'
        sizes=outfa+'.sizes'
        seen={}
        cur_orig=None; cur_new=None; length=0
        seq_count=0
        with opener(r['input_fasta']) as inp, open(tmp,'w') as of, open(sizes+'.tmp','w') as sf:
            def flush():

                if cur_new is not None:
                    sf.write(f'{cur_new}\t{length}\n')
                    hf.write(f"{r['assembly_id']}\t{r['internal_id']}\t{cur_orig}\t{cur_new}\t{length}\n")
            for line in inp:
                if line.startswith('>'):
                    flush()
                    cur_orig=line[1:].strip().split()[0]
                    c=clean_chr(cur_orig)
                    n=seen.get(c,0)+1; seen[c]=n
                    if n>1:
                        c=f'{c}x{n}'
                    cur_new=f"{r['internal_id']}.{c}"
                    of.write('>'+cur_new+'\n')
                    length=0; seq_count += 1
                else:
                    s=line.strip()
                    if not s: continue
                    of.write(s+'\n')
                    length += len(s)
            flush()
        os.replace(tmp,outfa)
        os.replace(sizes+'.tmp',sizes)
        if seq_count == 0:
            raise SystemExit(f'no sequences in {r["input_fasta"]}')
with open(os.path.join(out,'genome','g.lst'),'w') as gl:
    for r in rows:
        gl.write(r['internal_id']+'\n')
for ref in ('Msa','R108'):
    refdir=os.path.join(out,ref+'_ref')
    os.makedirs(os.path.join(refdir,'alignment'),exist_ok=True)
    os.makedirs(os.path.join(refdir,'chainnet'),exist_ok=True)
    with open(os.path.join(refdir,'query.lst'),'w') as qh:
        for r in rows:
            if r['internal_id'] != ref:
                qh.write(r['internal_id']+'\n')
tasks=[]
for ref in ('Msa','R108'):
    for r in rows:
        q=r['internal_id']
        if q == ref: continue
        task=f'{ref}_vs_{q}'
        refdir=os.path.join(out,ref+'_ref')
        tasks.append({
            'task_id':task,
            'ref_run':ref+'_ref',
            'ref_id':ref,
            'query_id':q,
            'query_assembly_id':r['assembly_id'],
            'source_group':r.get('source_group',''),
            'role':r.get('role',''),
            'pod_coiling':r.get('pod_coiling',''),
            'pod_spine':r.get('pod_spine',''),
            'independent_trait_sample':r.get('independent_trait_sample',''),
            'ref_fasta':os.path.join(out,'genome',ref),
            'query_fasta':os.path.join(out,'genome',q),
            'paf':os.path.join(refdir,'alignment',f'{ref}vs{q}.chr.paf')
        })
task_path=os.path.join(out,'01_metadata','svgap_wga_tasks.tsv')
fields=['task_id','ref_run','ref_id','query_id','query_assembly_id','source_group','role','pod_coiling','pod_spine','independent_trait_sample','ref_fasta','query_fasta','paf']
with open(task_path,'w',newline='') as th:
    th.write('\t'.join(fields)+'\n')
    for t in tasks:
        th.write('\t'.join(str(t[k]) for k in fields)+'\n')
print(f'prepared_genomes={len(rows)}')
print(f'wga_tasks={len(tasks)}')
print(f'map={map_path}')
print(f'tasks={task_path}')
PY
for f in "$OUT"/genome/*; do
  case "$f" in *.sizes|*.2bit|*.tmp) continue;; esac
  [ -f "$f" ] || continue
  if [ ! -s "$f.sizes" ]; then awk '/^>/{if(n){print name"\t"n}; name=substr($0,2); n=0; next}{n+=length($0)}END{if(name){print name"\t"n}}' "$f" > "$f.sizes"; fi
done
printf 'DONE prepare_svgap_inputs %s\n' "$(date '+%F %T %Z')"
