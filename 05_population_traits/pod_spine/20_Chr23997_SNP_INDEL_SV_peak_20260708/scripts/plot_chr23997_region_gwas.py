#!/usr/bin/env python3
import os, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_plink(path, start=88000000, end=90000000):
    pts=[]
    if not os.path.exists(path): return pts
    with open(path) as f:
        header=f.readline().split()
        idx={h:i for i,h in enumerate(header)}
        for line in f:
            if not line.strip(): continue
            p=line.split()
            try:
                pos=int(float(p[idx['POS']]))
                pv=p[idx['P']]
                if pv=='NA': continue
                pv=float(pv)
            except Exception:
                continue
            if start <= pos <= end and pv>0:
                pts.append((pos/1e6, -math.log10(pv)))
    return pts
RUN='path/to/project/N_4.pod_spiny/20_Chr23997_SNP_INDEL_SV_peak_20260708'
GWAS='path/to/project/N_4.pod_spiny/08_GWAS_pod_spine_SNP_SV_20260706'
files={
 'No covariates': [
   ('SNP', GWAS+'/results/sativa182_pod_spine_snp_gwas.corrected.no_covar.pod_spine.glm.logistic', '#6a51a3'),
   ('InDel', RUN+'/results/indel.Chr4_88_90Mb.no_covar.pod_spine.glm.logistic', '#e6550d'),
   ('Delly DEL', RUN+'/results/dellyDEL.Chr4_88_90Mb.no_covar.pod_spine.glm.logistic', '#1b9e77'),
 ],
 'PC1-5 corrected': [
   ('SNP', GWAS+'/results/sativa182_pod_spine_snp_gwas.corrected.PC5.pod_spine.glm.logistic', '#6a51a3'),
   ('InDel', RUN+'/results/indel.Chr4_88_90Mb.PC5.pod_spine.glm.logistic', '#e6550d'),
   ('Delly DEL', RUN+'/results/dellyDEL.Chr4_88_90Mb.PC5.pod_spine.glm.logistic', '#1b9e77'),
 ]
}
fig, axes=plt.subplots(2,1,figsize=(9,5.8),sharex=True)
gene_start=88979818/1e6; gene_end=88982331/1e6
for ax,(title,arr) in zip(axes,files.items()):
    ymax=0
    for name,path,color in arr:
        pts=read_plink(path)
        if pts:
            xs,ys=zip(*pts)
            ymax=max(ymax,max(ys))
            ax.scatter(xs,ys,s=9,alpha=0.7,label=f'{name} (n={len(pts)})',color=color,edgecolors='none')
        else:
            ax.scatter([],[],label=f'{name} (n=0)',color=color)
    ax.axvspan(gene_start,gene_end,color='red',alpha=0.15,label='Chr23997' if title=='No covariates' else None)
    ax.axvline(gene_start,color='red',lw=0.8,alpha=0.5)
    ax.axvline(gene_end,color='red',lw=0.8,alpha=0.5)
    ax.set_title(title,fontsize=11)
    ax.set_ylabel(r'$-log_{10}(P)$')
    ax.set_ylim(bottom=0, top=max(4, ymax*1.12 if ymax else 4))
    ax.grid(axis='y',color='#dddddd',lw=0.5)
axes[-1].set_xlabel('Chr4 position (Mb)')
axes[-1].set_xlim(88,90)
axes[0].legend(frameon=False,ncol=4,fontsize=9,loc='upper right')
fig.suptitle('Pod spine association around Chr23997 (Chr4:88-90 Mb)', y=0.99, fontsize=12)
fig.tight_layout(rect=[0,0,1,0.96])
outbase=RUN+'/figures/Chr23997_SNP_INDEL_DELLY_Chr4_88_90Mb_zoom'
fig.savefig(outbase+'.png',dpi=300)
fig.savefig(outbase+'.pdf')
print(outbase+'.png')
