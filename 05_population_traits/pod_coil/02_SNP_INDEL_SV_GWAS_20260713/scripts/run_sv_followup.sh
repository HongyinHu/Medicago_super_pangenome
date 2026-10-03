#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/02_SNP_INDEL_SV_GWAS_20260713
SOURCE=path/to/project/N_4.pod_spiny/28_GWAS_spiny
SVSRC="$SOURCE/05_sv_genotyping"
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
TABIX=path/to/home/anaconda3/envs/panpop/bin/tabix
SMOOVE=path/to/home/anaconda3/envs/smoove/bin/smoove
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python

export PATH=path/to/home/anaconda3/envs/smoove/bin:path/to/home/anaconda3/envs/panpop/bin:path/to/home/anaconda3/envs/biosofeware/bin:/usr/bin:/bin

mkdir -p "$ROOT/03_sv_gwas"/{joint,plink,gemma,results,summary} "$ROOT/logs" "$ROOT/tmp"
LOG="$ROOT/logs/run_sv_followup.log"
exec > >(tee -a "$LOG") 2>&1

echo "[$(date '+%F %T %Z')] SV follow-up waiting for joint candidate discovery"
while [[ ! -s "$SVSRC/summary/sv_site_discovery.done" ]]; do
  if grep -Eqi 'No space left|Disk quota exceeded|Traceback|fatal|killed' "$SVSRC/logs/run_sv_site_discovery.log"; then
    echo "Fatal pattern detected in SV discovery log" >&2
    exit 2
  fi
  sleep 300
done

count_done() {
  find "$SVSRC/joint_genotype_per_sample" -mindepth 2 -maxdepth 2 -name '*.done' -type f 2>/dev/null | wc -l
}

retries=0
while (( $(count_done) < 143 )); do
  active=$(squeue -h -u $USER -n spiny_sv_gt 2>/dev/null | wc -l)
  if (( active == 0 )); then
    if (( retries >= 3 )); then
      echo "SV joint genotyping incomplete after three submissions: $(count_done)/143" >&2
      exit 3
    fi
    job=$(bash "$SVSRC/scripts/submit_sv_genotype_array.sh")
    retries=$((retries + 1))
    printf '%s\t%s\t%s\n' "$(date '+%F %T %Z')" "$job" "$retries" >> "$ROOT/03_sv_gwas/summary/genotype_submissions.tsv"
    echo "Submitted SV genotype array $job (attempt $retries)"
  fi
  echo "[$(date '+%F %T')] SV genotype complete $(count_done)/143"
  sleep 300
done

MANIFEST="$SVSRC/summary/sv_analysis_manifest.tsv"
VCF_LIST="$ROOT/tmp/joint_genotype_vcfs.list"
awk -F '\t' 'NR>1 {print $1}' "$MANIFEST" | while read -r sid; do
  vcf="$SVSRC/joint_genotype_per_sample/$sid/$sid-smoove.genotyped.vcf.gz"
  [[ -s "$vcf" ]] || { echo "Missing genotype VCF: $vcf" >&2; exit 1; }
  printf '%s\n' "$vcf"
done > "$VCF_LIST"
[[ $(wc -l < "$VCF_LIST") -eq 143 ]]

JOINT_DIR="$ROOT/03_sv_gwas/joint"
JOINT="$JOINT_DIR/medicago143.smoove.square.vcf.gz"
if [[ ! -s "$JOINT" ]]; then
  mapfile -t vcfs < "$VCF_LIST"
  "$SMOOVE" paste --name medicago143 --outdir "$JOINT_DIR" "${vcfs[@]}"
fi
[[ -s "$JOINT" ]]
if [[ ! -s "$JOINT.csi" && ! -s "$JOINT.tbi" ]]; then
  "$BCFTOOLS" index -f "$JOINT"
fi
[[ $("$BCFTOOLS" query -l "$JOINT" | wc -l) -eq 143 ]]
"$BCFTOOLS" index -n "$JOINT" >/dev/null

MAP="$ROOT/tmp/chrom_rename.tsv"
printf 'Chr1\t1\nChr2\t2\nChr3\t3\nChr4\t4\nChr5\t5\nChr6\t6\nChr7\t7\nChr8\t8\n' > "$MAP"
STRICT_VCF="$JOINT_DIR/coil125.SV.GQ20.vcf.gz"
if [[ ! -s "$STRICT_VCF" ]]; then
  "$BCFTOOLS" view -S "$ROOT/00_input/strict_sample_ids.txt" \
    -r Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8 -Ou "$JOINT" \
    | "$BCFTOOLS" filter -S . -e 'FMT/GQ<20' -Ou \
    | "$BCFTOOLS" annotate --rename-chrs "$MAP" --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o "$STRICT_VCF"
  "$TABIX" -f -p vcf "$STRICT_VCF"
fi

SV_PLINK_STAGE1="$ROOT/03_sv_gwas/plink/coil.SV.GQ20.geno20"
if [[ ! -s "$SV_PLINK_STAGE1.bed" ]]; then
  "$PLINK" --vcf "$STRICT_VCF" --double-id --chr-set 8 no-xy \
    --biallelic-only strict --geno 0.20 \
    --make-bed --out "$SV_PLINK_STAGE1"
fi

SV_PLINK="$ROOT/03_sv_gwas/plink/coil.SV.GQ20.twostage"
if [[ ! -s "$SV_PLINK.bed" ]]; then
  "$PLINK" --bfile "$SV_PLINK_STAGE1" --chr-set 8 no-xy \
    --mind 0.20 --maf 0.05 \
    --make-bed --out "$SV_PLINK"
fi
SV_N=$(wc -l < "$SV_PLINK.fam")
(( SV_N >= 100 )) || {
  echo "SV two-stage QC retained too few samples: $SV_N" >&2
  exit 4
}

"$PYTHON" "$ROOT/scripts/reorder_kinship_for_sv.py" \
  "$ROOT/01_structure/plink/coil.SNP.fam" \
  "$ROOT/01_structure/gemma/output/coil_SNP_kinship.cXX.txt" \
  "$SV_PLINK.fam" "$ROOT/00_input/coil_phenotype_strict.tsv" \
  "$ROOT/01_structure/plink/coil.SNP.pca.eigenvec" \
  "$ROOT/03_sv_gwas/summary" 5

SV_ASSOC="$ROOT/03_sv_gwas/gemma/output/SV.coil.pc5.gq20.twostage.assoc.txt"
if [[ ! -s "$SV_ASSOC" ]]; then
  (
    cd "$ROOT/03_sv_gwas/gemma"
    "$GEMMA" -bfile "$SV_PLINK" \
      -k "$ROOT/03_sv_gwas/summary/SV_kinship.reordered.cXX.txt" \
      -c "$ROOT/03_sv_gwas/summary/SV.covariates.pc5.txt" \
      -p "$ROOT/03_sv_gwas/summary/SV.phenotype.0_1.txt" \
      -lmm 4 -o SV.coil.pc5.gq20.twostage
  )
fi
test -s "$SV_ASSOC"

"$PYTHON" "$ROOT/scripts/summarize_gwas.py" \
  --assoc "SV=$SV_ASSOC" --outdir "$ROOT/03_sv_gwas/results" --prefix coil_SV

"$PYTHON" "$ROOT/scripts/summarize_gwas.py" \
  --assoc "SNP=$ROOT/02_snp_indel_gwas/gemma/output/SNP.coil.pc5.assoc.txt" \
  --assoc "INDEL=$ROOT/02_snp_indel_gwas/gemma/output/INDEL.coil.pc5.assoc.txt" \
  --assoc "SV=$SV_ASSOC" --outdir "$ROOT/04_candidate_gene_peaks" --prefix coil_SNP_INDEL_SV

{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'samples_after_sv_qc\t%s\n' "$SV_N"
  printf 'joint_source_samples\t143\n'
  printf 'phenotype\t0=non_spiral,1=spiral; weak_spiral_excluded\n'
  printf 'genotype_qc\tGQ>=20; two-stage marker missingness<=0.20 then sample missingness<=0.20; MAF>=0.05\n'
  printf 'model\tGEMMA_LMM_Wald; SNP_kinship; PC1-PC5\n'
} > "$ROOT/summary/sv_gwas.done"

echo "[$(date '+%F %T %Z')] SNP/INDEL/SV coil GWAS completed"
