#!/usr/bin/env bash
set -eo pipefail

BASE="path/to/project/38.medicago_resequence/test_population"
VCF="$BASE/DP_7-50_miss_0.2_sativa182.SNP.vcf.gz"
OUT="$BASE/01.population_analysis"
THREADS="${THREADS:-16}"
KMIN="${KMIN:-2}"
KMAX="${KMAX:-10}"

mkdir -p "$OUT"/{00_input,01_plink,02_pca,03_admixture,04_tree,05_plots,06_summary,99_scripts,logs}
LOG="$OUT/logs/run.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

source ~/.bashrc 2>/dev/null || true
conda activate biosofeware
export OMP_NUM_THREADS="$THREADS"

cd "$OUT"
ln -sf "$VCF" 00_input/input.vcf.gz

printf '[%s] Start population analysis on %s\n' "$(date '+%F %T')" "$(hostname)"
printf 'VCF=%s\nOUT=%s\nTHREADS=%s\nK=%s..%s\n' "$VCF" "$OUT" "$THREADS" "$KMIN" "$KMAX"

python - <<'PY'
import gzip, pathlib
vcf = pathlib.Path('path/to/project/38.medicago_resequence/test_population/DP_7-50_miss_0.2_sativa182.SNP.vcf.gz')
out = pathlib.Path('path/to/project/38.medicago_resequence/test_population/01.population_analysis/06_summary')
out.mkdir(parents=True, exist_ok=True)
with gzip.open(vcf, 'rt', errors='replace') as fh:
    for line in fh:
        if line.startswith('#CHROM'):
            samples = line.rstrip('\n').split('\t')[9:]
            break
    else:
        raise SystemExit('No #CHROM header found')
(out / 'samples.txt').write_text('\n'.join(samples) + '\n')
(out / 'sample_count.txt').write_text(str(len(samples)) + '\n')
print(f'Samples: {len(samples)}')
PY

if [[ ! -s 01_plink/sativa182.filtered.bed ]]; then
  printf '[%s] Convert/filter VCF to PLINK BED\n' "$(date '+%F %T')"
  plink \
    --vcf "$VCF" \
    --double-id \
    --allow-extra-chr \
    --keep-allele-order \
    --vcf-half-call m \
    --biallelic-only strict \
    --snps-only just-acgt \
    --geno 0.2 \
    --maf 0.05 \
    --set-missing-var-ids @:# \
    --make-bed \
    --out 01_plink/sativa182.filtered
fi

if [[ ! -s 01_plink/sativa182.ld.prune.in ]]; then
  printf '[%s] LD pruning\n' "$(date '+%F %T')"
  plink --bfile 01_plink/sativa182.filtered --allow-extra-chr --indep-pairwise 50 10 0.2 --out 01_plink/sativa182.ld
fi

if [[ ! -s 01_plink/sativa182.pruned.bed ]]; then
  printf '[%s] Make pruned PLINK BED\n' "$(date '+%F %T')"
  plink --bfile 01_plink/sativa182.filtered --allow-extra-chr --extract 01_plink/sativa182.ld.prune.in --make-bed --out 01_plink/sativa182.pruned
fi

if [[ ! -s 02_pca/sativa182.pca.eigenvec ]]; then
  printf '[%s] PCA\n' "$(date '+%F %T')"
  plink --bfile 01_plink/sativa182.pruned --allow-extra-chr --pca 20 --out 02_pca/sativa182.pca
fi

python - <<'PY'
from pathlib import Path
p = Path('path/to/project/38.medicago_resequence/test_population/01.population_analysis')
eig = p/'02_pca/sativa182.pca.eigenvec'
out = p/'02_pca/sativa182.pca.tsv'
if eig.exists():
    with eig.open() as fi, out.open('w') as fo:
        header = ['FID','IID'] + [f'PC{i}' for i in range(1,21)]
        fo.write('\t'.join(header)+'\n')
        for line in fi:
            cols = line.split()
            fo.write('\t'.join(cols)+'\n')
PY

cat > 99_scripts/plot_pca.R <<'RS'
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
eigval <- args[2]
outprefix <- args[3]
d <- read.table(infile, header=TRUE, sep="\t", check.names=FALSE)
varlab <- c("PC1", "PC2")
if (file.exists(eigval)) {
  ev <- scan(eigval, quiet=TRUE)
  if (length(ev) >= 2 && sum(ev) > 0) {
    varlab <- paste0(c("PC1", "PC2"), " (", round(ev[1:2] / sum(ev) * 100, 2), "%)")
  }
}
pdf(paste0(outprefix, ".pdf"), width=6.5, height=5.5)
plot(d$PC1, d$PC2, pch=21, bg="#2c7fb8", col="#08306b", xlab=varlab[1], ylab=varlab[2], main="PCA")
dev.off()
png(paste0(outprefix, ".png"), width=1800, height=1500, res=250)
plot(d$PC1, d$PC2, pch=21, bg="#2c7fb8", col="#08306b", xlab=varlab[1], ylab=varlab[2], main="PCA")
dev.off()
pdf(paste0(outprefix, ".labels.pdf"), width=8, height=7)
plot(d$PC1, d$PC2, pch=21, bg="#2c7fb8", col="#08306b", xlab=varlab[1], ylab=varlab[2], main="PCA with sample IDs")
text(d$PC1, d$PC2, labels=d$IID, cex=0.45, pos=3)
dev.off()
RS
Rscript 99_scripts/plot_pca.R 02_pca/sativa182.pca.tsv 02_pca/sativa182.pca.eigenval 02_pca/sativa182.pca
cp -f 02_pca/sativa182.pca.pdf 05_plots/ 2>/dev/null || true
cp -f 02_pca/sativa182.pca.png 05_plots/ 2>/dev/null || true
cp -f 02_pca/sativa182.pca.labels.pdf 05_plots/ 2>/dev/null || true

if [[ ! -s 03_admixture/sativa182.admix.bed ]]; then
  printf '[%s] Prepare ADMIXTURE BED with numeric chromosome IDs\n' "$(date '+%F %T')"
  cp 01_plink/sativa182.pruned.bed 03_admixture/sativa182.admix.bed
  cp 01_plink/sativa182.pruned.fam 03_admixture/sativa182.admix.fam
  awk 'BEGIN{OFS="\t"}{$1=1; print}' 01_plink/sativa182.pruned.bim > 03_admixture/sativa182.admix.bim
fi

cd 03_admixture
for K in $(seq "$KMIN" "$KMAX"); do
  if [[ ! -s sativa182.admix.${K}.Q ]]; then
    printf '[%s] ADMIXTURE K=%s\n' "$(date '+%F %T')" "$K"
    admixture --cv=10 -j"$THREADS" sativa182.admix.bed "$K" > "K${K}.log" 2>&1
    mv -f sativa182.admix.${K}.P "sativa182.admix.${K}.P" 2>/dev/null || true
    mv -f sativa182.admix.${K}.Q "sativa182.admix.${K}.Q" 2>/dev/null || true
  fi
done
cd "$OUT"

grep -H "CV error" 03_admixture/K*.log 2>/dev/null | sed -E 's#03_admixture/K([0-9]+).log:CV error \(K=[0-9]+\): #K\t#' > 03_admixture/admixture_cv_errors.tsv || true

cat > 99_scripts/plot_admixture.R <<'RS'
args <- commandArgs(trailingOnly=TRUE)
fam <- args[1]
outdir <- args[2]
kmin <- as.integer(args[3]); kmax <- as.integer(args[4])
f <- read.table(fam, stringsAsFactors=FALSE)
ids <- f[,2]
cv <- data.frame(K=integer(), CV=numeric())
for (K in kmin:kmax) {
  qfile <- file.path(outdir, paste0("sativa182.admix.", K, ".Q"))
  if (!file.exists(qfile)) next
  q <- as.matrix(read.table(qfile))
  ord <- order(max.col(q, ties.method="first"), ids)
  pdf(file.path(outdir, paste0("admixture.K", K, ".pdf")), width=10, height=3.8)
  par(mar=c(5,4,2,1))
  barplot(t(q[ord,,drop=FALSE]), col=rainbow(K), border=NA, space=0, names.arg=ids[ord], las=2, cex.names=0.32, ylab="Ancestry proportion", main=paste0("ADMIXTURE K=", K))
  dev.off()
  png(file.path(outdir, paste0("admixture.K", K, ".png")), width=2600, height=900, res=220)
  par(mar=c(5,4,2,1))
  barplot(t(q[ord,,drop=FALSE]), col=rainbow(K), border=NA, space=0, names.arg=ids[ord], las=2, cex.names=0.32, ylab="Ancestry proportion", main=paste0("ADMIXTURE K=", K))
  dev.off()
}
logs <- list.files(outdir, pattern="^K[0-9]+\\.log$", full.names=TRUE)
if (length(logs)) {
  for (lf in logs) {
    txt <- readLines(lf, warn=FALSE)
    hit <- grep("CV error", txt, value=TRUE)
    if (length(hit)) {
      K <- as.integer(sub("^K([0-9]+)\\.log$", "\\1", basename(lf)))
      val <- as.numeric(sub(".*: ", "", hit[length(hit)]))
      cv <- rbind(cv, data.frame(K=K, CV=val))
    }
  }
  if (nrow(cv)) {
    cv <- cv[order(cv$K),]
    write.table(cv, file=file.path(outdir, "admixture_cv_errors.clean.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    pdf(file.path(outdir, "admixture_cv_errors.pdf"), width=5.5, height=4.5)
    plot(cv$K, cv$CV, type="b", pch=19, xlab="K", ylab="Cross-validation error", main="ADMIXTURE CV error")
    dev.off()
    png(file.path(outdir, "admixture_cv_errors.png"), width=1300, height=1000, res=220)
    plot(cv$K, cv$CV, type="b", pch=19, xlab="K", ylab="Cross-validation error", main="ADMIXTURE CV error")
    dev.off()
  }
}
RS
Rscript 99_scripts/plot_admixture.R 03_admixture/sativa182.admix.fam 03_admixture "$KMIN" "$KMAX"
cp -f 03_admixture/admixture*.pdf 03_admixture/admixture*.png 05_plots/ 2>/dev/null || true

if [[ ! -s 04_tree/sativa182.pruned.vcf ]]; then
  printf '[%s] Export pruned VCF for tree\n' "$(date '+%F %T')"
  plink --bfile 01_plink/sativa182.pruned --allow-extra-chr --recode vcf-iid --out 04_tree/sativa182.pruned
fi

cat > 99_scripts/vcf_to_phy.py <<'PY'
#!/usr/bin/env python3
from pathlib import Path
import gzip
vcf = Path('path/to/project/38.medicago_resequence/test_population/01.population_analysis/04_tree/sativa182.pruned.vcf')
out_phy = Path('path/to/project/38.medicago_resequence/test_population/01.population_analysis/04_tree/sativa182.pruned.phy')
out_fa = Path('path/to/project/38.medicago_resequence/test_population/01.population_analysis/04_tree/sativa182.pruned.fasta')
iupac = {frozenset('AG'):'R', frozenset('CT'):'Y', frozenset('CG'):'S', frozenset('AT'):'W', frozenset('GT'):'K', frozenset('AC'):'M'}
openfn = gzip.open if vcf.suffix == '.gz' else open
samples = []
seqs = []
nvar = 0
with openfn(vcf, 'rt', errors='replace') as fh:
    for line in fh:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            samples = line.rstrip('\n').split('\t')[9:]
            seqs = [[] for _ in samples]
            continue
        if not samples:
            continue
        cols = line.rstrip('\n').split('\t')
        ref, alt = cols[3].upper(), cols[4].upper()
        if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT':
            continue
        alleles = [ref, alt]
        for i, gtfield in enumerate(cols[9:]):
            gt = gtfield.split(':', 1)[0].replace('|', '/')
            if gt in ('0/0', '0'):
                base = ref
            elif gt in ('1/1', '1'):
                base = alt
            elif gt in ('0/1', '1/0'):
                base = iupac.get(frozenset((ref, alt)), 'N')
            else:
                base = 'N'
            seqs[i].append(base)
        nvar += 1
if not samples or nvar == 0:
    raise SystemExit('No usable SNPs found for tree alignment')
with out_phy.open('w') as fo:
    fo.write(f'{len(samples)} {nvar}\n')
    for name, seq in zip(samples, seqs):
        safe = ''.join(c if c.isalnum() or c in '_.-' else '_' for c in name)[:30]
        fo.write(f'{safe:<32} {"".join(seq)}\n')
with out_fa.open('w') as fo:
    for name, seq in zip(samples, seqs):
        safe = ''.join(c if c.isalnum() or c in '_.-' else '_' for c in name)
        fo.write(f'>{safe}\n')
        s = ''.join(seq)
        for j in range(0, len(s), 80):
            fo.write(s[j:j+80] + '\n')
print(f'Wrote {len(samples)} samples x {nvar} SNPs')
PY
chmod +x 99_scripts/vcf_to_phy.py
if [[ ! -s 04_tree/sativa182.pruned.phy ]]; then
  printf '[%s] Convert VCF to PHYLIP/FASTA\n' "$(date '+%F %T')"
  python 99_scripts/vcf_to_phy.py
fi

if [[ ! -s 04_tree/iqtree/sativa182_pruned.treefile ]]; then
  printf '[%s] IQ-TREE ML tree\n' "$(date '+%F %T')"
  mkdir -p 04_tree/iqtree
  iqtree2 -s 04_tree/sativa182.pruned.phy -m GTR+ASC -bb 1000 -nt "$THREADS" -pre 04_tree/iqtree/sativa182_pruned
fi

cat > 06_summary/methods.txt <<EOF
Input VCF: $VCF
Conda env: biosofeware
Filtering: PLINK --biallelic-only strict --snps-only just-acgt --geno 0.2 --maf 0.05
LD pruning: PLINK --indep-pairwise 50 10 0.2
PCA: PLINK --pca 20 on LD-pruned SNPs
Population structure: ADMIXTURE K=$KMIN..$KMAX, --cv=10, on LD-pruned SNPs with numeric chromosome IDs
ML tree: IQ-TREE2, SNP-only PHYLIP from LD-pruned SNPs, model GTR+ASC, ultrafast bootstrap 1000
EOF

{
  echo -e "item\tvalue"
  echo -e "sample_count\t$(cat 06_summary/sample_count.txt)"
  [[ -s 01_plink/sativa182.filtered.bim ]] && echo -e "filtered_snp_count\t$(wc -l < 01_plink/sativa182.filtered.bim)"
  [[ -s 01_plink/sativa182.ld.prune.in ]] && echo -e "ld_pruned_snp_count\t$(wc -l < 01_plink/sativa182.ld.prune.in)"
  [[ -s 04_tree/sativa182.pruned.phy ]] && awk 'NR==1{print "tree_alignment\t"$1" samples x "$2" SNPs"}' 04_tree/sativa182.pruned.phy
  [[ -s 03_admixture/admixture_cv_errors.clean.tsv ]] && awk 'NR>1{if(min=="" || $2<min){min=$2;k=$1}} END{if(k!="") print "best_K_by_CV\t"k" (CV="min")"}' 03_admixture/admixture_cv_errors.clean.tsv
} > 06_summary/summary.tsv

printf '[%s] Done. Summary: %s\n' "$(date '+%F %T')" "$OUT/06_summary/summary.tsv"
