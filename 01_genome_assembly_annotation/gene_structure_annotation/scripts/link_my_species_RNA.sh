#!/usr/bin/env bash
set -euo pipefail
ROOT=path/to/project/pan_genome_new
TARGET=path/to/project/N_1.coding_gene_anno/00_inputs/my_species_RNA
mkdir -p "$TARGET"
case "$TARGET" in
  path/to/project/N_1.coding_gene_anno/00_inputs/my_species_RNA) ;;
  *) echo "Unsafe TARGET: $TARGET" >&2; exit 2;;
esac
# Clean only symlinks from previous organization attempts; never delete raw files or real files.
find "$TARGET" -type l -delete 2>/dev/null || true
manifest="$TARGET/manifest.tsv"
summary="$TARGET/summary.tsv"
missing="$TARGET/missing_samples.tsv"
tmp_manifest="$manifest.tmp"
tmp_summary="$summary.tmp"
tmp_missing="$missing.tmp"
printf 'sample\tevidence_type\tlink_path\ttarget_path\ttarget_size\ttarget_mtime\traw_species_dir\tnote\n' > "$tmp_manifest"
printf 'sample\traw_species_dir\tillumina_files\tiso_raw_files\tstatus\tnote\n' > "$tmp_summary"
printf 'sample\treason\tnote\n' > "$tmp_missing"

safe_rel_name(){
  base_dir="$1"; file="$2"
  rel=${file#"$base_dir"/}
  rel=${rel//\//__}
  printf '%s' "$rel"
}

link_one(){
  sample="$1"; type="$2"; rawdir="$3"; base_dir="$4"; file="$5"; note="$6"; subdir="$TARGET/$sample/$type"
  [ -s "$file" ] || return 0
  mkdir -p "$subdir"
  link_name=$(safe_rel_name "$base_dir" "$file")
  link="$subdir/$link_name"
  if [ -e "$link" ] && [ ! -L "$link" ]; then
    printf '%s\tconflict_non_symlink\t%s\n' "$sample" "$link" >> "$tmp_missing"
    return 0
  fi
  ln -sfn "$file" "$link"
  size=$(stat -L -c %s "$file")
  mtime=$(stat -L -c '%y' "$file" | cut -d. -f1)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$type" "$link" "$file" "$size" "$mtime" "$rawdir" "$note" >> "$tmp_manifest"
}

add_links(){
  sample="$1"; rawdir="$2"; note="${3:-}"
  mkdir -p "$TARGET/$sample/illumina-seq" "$TARGET/$sample/iso-seq"
  illum_count=0
  iso_count=0
  status="ok"
  if [ ! -d "$rawdir" ]; then
    status="missing_rawdir"
    printf '%s\t%s\t%s\n' "$sample" "$status" "$rawdir" >> "$tmp_missing"
    printf '%s\t%s\t0\t0\t%s\t%s\n' "$sample" "$rawdir" "$status" "$note" >> "$tmp_summary"
    return 0
  fi

  while IFS= read -r d; do
    while IFS= read -r f; do
      link_one "$sample" illumina-seq "$rawdir" "$d" "$f" "$note"
      illum_count=$((illum_count+1))
    done < <(find -L "$d" -maxdepth 3 -type f \( -iname '*.fq.gz' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fastq' \) | sort)
  done < <(find -L "$rawdir" -maxdepth 2 -type d \( -iname '*illumina*RNA*' -o -iname '*RNA*illumina*' \) | sort)

  # For Iso-Seq, keep only raw inputs needed to reprocess from scratch.
  # Do not link 4_iso/output/* products such as ccs/flnc/transcripts/lima/demux intermediates.
  while IFS= read -r d; do
    while IFS= read -r f; do
      b=$(basename "$f")
      case "$b" in pb.fasta|pb.fasta.gz) continue ;; esac
      link_one "$sample" iso-seq "$rawdir" "$d" "$f" "$note"
      iso_count=$((iso_count+1))
    done < <(find -L "$d" -maxdepth 1 -type f \( -iname '*.bam' -o -iname '*.pbi' -o -iname '*.bai' -o -iname 'primer.fa' -o -iname '*.fa' -o -iname '*.fasta' -o -iname '*.fa.gz' -o -iname '*.fasta.gz' \) | sort)
  done < <(find -L "$rawdir" -maxdepth 2 -type d \( -iname '*iso*' -o -iname '*Iso*' \) | sort)

  if [ "$illum_count" -eq 0 ] && [ "$iso_count" -eq 0 ]; then
    status="no_rna_files_found"
    printf '%s\t%s\t%s\n' "$sample" "$status" "$rawdir" >> "$tmp_missing"
  elif [ "$illum_count" -eq 0 ]; then
    status="iso_only"
  elif [ "$iso_count" -eq 0 ]; then
    status="illumina_only"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$rawdir" "$illum_count" "$iso_count" "$status" "$note" >> "$tmp_summary"
}

add_links genome_395  "$ROOT/X_genome_TLMX_Mu_395_1/0.genome_raw_data" "mapped_from_X_genome_TLMX_Mu_395_1"
add_links genome_410  "$ROOT/N_genome_410_0/0.genome_raw_data" "mapped_from_N_genome_410_0"
add_links genome_436  "$ROOT/X_genome_436_0/0.genome_raw_data" "mapped_from_X_genome_436_0"
add_links genome_454  "$ROOT/X_genome_MS454_mianmaomuxu_0/0.genome_raw_data" "mapped_from_X_genome_MS454_mianmaomuxu_0"
add_links genome_457  "$ROOT/N_genome_Ms457_0/0.genome_raw_data" "mapped_from_N_genome_Ms457_0"
add_links genome_461  "$ROOT/N_genome_Ms461_0/0.genome_raw_data" "mapped_from_N_genome_Ms461_0"
add_links genome_468  "$ROOT/X_genome_468_1/0.genome_raw_data" "mapped_from_X_genome_468_1"
add_links genome_472  "$ROOT/N_genome_472_0/0.genome_raw_data" "mapped_from_N_genome_472_0"
add_links genome_474a "$ROOT/N_genome_474_0/0.genome_raw_data" "mapped_from_N_genome_474_0_shared_with_474b"
add_links genome_474b "$ROOT/N_genome_474_0/0.genome_raw_data" "mapped_from_N_genome_474_0_shared_with_474a"
add_links genome_482  "$ROOT/X_genome_MS482_0/0.genome_raw_data" "mapped_from_X_genome_MS482_0"
add_links genome_M22  "$ROOT/X_genome_M22_0/0.genome_raw_data" "mapped_from_X_genome_M22_0"
add_links genome_M46  "$ROOT/X_genome_M46_1/0.genome_raw_data" "mapped_from_X_genome_M46_1"
add_links genome_Mar  "$ROOT/P_genome_Medicago_archiducis-nicolai/0.genome_raw_data" "mapped_from_P_genome_Medicago_archiducis-nicolai"
add_links genome_Mru  "$ROOT/P_genome_Medicago_ruthenica/0.genome_raw_data" "mapped_from_P_genome_Medicago_ruthenica"

for s in genome_A17 genome_Mpo genome_Msa1 genome_Msa2 genome_R108 genome_ZM4; do
  mkdir -p "$TARGET/$s/illumina-seq" "$TARGET/$s/iso-seq"
  printf '%s\tno_0.genome_raw_data_found\tsearch_root=%s\n' "$s" "$ROOT" >> "$tmp_missing"
  printf '%s\tNA\t0\t0\tmissing_rawdir\tno matching 0.genome_raw_data found under pan_genome_new\n' "$s" >> "$tmp_summary"
done

mv "$tmp_manifest" "$manifest"
mv "$tmp_summary" "$summary"
mv "$tmp_missing" "$missing"
echo "manifest=$manifest"
echo "summary=$summary"
echo "missing=$missing"
