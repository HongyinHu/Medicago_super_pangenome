#!/usr/bin/env bash
set -u

for tool in \
  bcftools bgzip tabix samtools bedtools \
  plink plink2 gemma gcta64 Rscript python3 \
  delly smoove svtyper paragraph genotype_sv.py \
  survivor SURVIVOR java; do
  path=$(command -v "$tool" 2>/dev/null || true)
  printf '%s\t%s\n' "$tool" "${path:-NOT_FOUND}"
done

printf 'SEARCHED_ENV_EXECUTABLES\n'
for executable in \
  gemma gcta64 delly paragraph multigrmpy.py graphtyper \
  plink plink2 Rscript; do
  find path/to/home/anaconda3/envs -maxdepth 3 -type f \
    -name "$executable" -printf '%f\t%p\n' 2>/dev/null || true
done

for env_bin in \
  path/to/home/anaconda3/envs/panpop/bin \
  path/to/home/anaconda3/envs/sv/bin \
  path/to/home/anaconda3/envs/smoove/bin \
  path/to/home/anaconda3/envs/call_snp/bin; do
  if [[ -d "$env_bin" ]]; then
    printf 'ENV\t%s\n' "$env_bin"
    find "$env_bin" -maxdepth 1 -type f \
      \( -name 'plink*' -o -name 'gemma*' -o -name 'gcta*' -o \
         -name 'paragraph*' -o -name 'delly*' -o -name 'svtyper*' -o \
         -name 'bcftools' -o -name 'samtools' -o -name 'bedtools' \) \
      -printf '%f\t%p\n' 2>/dev/null || true
  fi
done
