#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/N_1.coding_gene_anno
TARGET="$BASE/00_inputs/my_species_RNA"
ROOT=path/to/project/0.raw_data
EXTRA_MANIFEST="$TARGET/manifest_project_T2T_extra.tsv"
EXTRA_SUMMARY="$TARGET/summary_project_T2T_extra.tsv"
EXTRA_MISSING="$TARGET/missing_project_T2T_extra.tsv"
case "$TARGET" in
  path/to/project/N_1.coding_gene_anno/00_inputs/my_species_RNA) ;;
  *) echo "Unsafe TARGET: $TARGET" >&2; exit 2;;
esac
mkdir -p "$TARGET"
# Remove only previous project_T2T supplemental symlinks; preserve pan_genome_new links.
find "$TARGET" -type l -name 'project_T2T__*' -delete 2>/dev/null || true
printf 'sample\tevidence_type\tlink_path\ttarget_path\ttarget_size\ttarget_mtime\traw_species_dir\tnote\n' > "$EXTRA_MANIFEST.tmp"
printf 'sample\traw_species_dir\tillumina_or_shortRNA_files\tlongread_transcript_raw_files\tstatus\tnote\n' > "$EXTRA_SUMMARY.tmp"
printf 'sample\treason\tnote\n' > "$EXTRA_MISSING.tmp"

safe_name(){
  base="$1"; file="$2"
  rel=${file#"$base"/}
  rel=${rel//\//__}
  printf 'project_T2T__%s' "$rel"
}
link_one(){
  sample="$1"; evidence="$2"; rawdir="$3"; base="$4"; file="$5"; note="$6"
  [ -s "$file" ] || return 0
  subdir="$TARGET/$sample/$evidence"
  mkdir -p "$subdir"
  link="$subdir/$(safe_name "$base" "$file")"
  ln -sfn "$file" "$link"
  size=$(stat -L -c %s "$file")
  mtime=$(stat -L -c '%y' "$file" | cut -d. -f1)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$evidence" "$link" "$file" "$size" "$mtime" "$rawdir" "$note" >> "$EXTRA_MANIFEST.tmp"
}

add_t2t(){
  sample="$1"; rawdir="$2"; note="$3"
  mkdir -p "$TARGET/$sample/illumina-seq" "$TARGET/$sample/iso-seq"
  short_count=0
  long_count=0
  if [ ! -d "$rawdir" ]; then
    printf '%s\tmissing_rawdir\t%s\n' "$sample" "$rawdir" >> "$EXTRA_MISSING.tmp"
    printf '%s\t%s\t0\t0\tmissing_rawdir\t%s\n' "$sample" "$rawdir" "$note" >> "$EXTRA_SUMMARY.tmp"
    return 0
  fi

  # Short RNA / public RNA-seq evidence. Only public RNA-seq fastq is treated as short-read here;
  # public RNA BAMs are kept because no raw fastq is present for that sample.
  for d in "$rawdir"/data_RNA_public; do
    [ -d "$d" ] || continue
    while IFS= read -r f; do
      link_one "$sample" illumina-seq "$rawdir" "$d" "$f" "$note;short_RNA_or_public_RNA"
      short_count=$((short_count+1))
    done < <(find -L "$d" -maxdepth 2 -type f \( -iname '*.fq.gz' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fastq' -o -iname '*.bam' -o -iname '*.bai' \) | sort)
  done

  # Long-read transcript evidence: PacBio RNA fasta, ONT-cDNA fastq, raw transcript BAM/fasta.
  # data_RNA is long-read transcript input in this project_T2T raw tree.
  for d in "$rawdir"/data_RNA "$rawdir"/data_RNA_public_PacBio "$rawdir"/data_ONT-cDNA_public; do
    [ -d "$d" ] || continue
    while IFS= read -r f; do
      b=$(basename "$f")
      case "$b" in pb.fasta|pb.fasta.gz) continue ;; esac
      link_one "$sample" iso-seq "$rawdir" "$d" "$f" "$note;long_read_transcript_raw"
      long_count=$((long_count+1))
    done < <(find -L "$d" -maxdepth 2 -type f \( -iname '*.bam' -o -iname '*.pbi' -o -iname '*.bai' -o -iname '*.fa' -o -iname '*.fasta' -o -iname '*.fa.gz' -o -iname '*.fasta.gz' -o -iname '*.fq.gz' -o -iname '*.fastq.gz' \) | sort)
  done

  status=ok
  if [ "$short_count" -eq 0 ] && [ "$long_count" -eq 0 ]; then
    status=no_rna_files_found
    printf '%s\t%s\t%s\n' "$sample" "$status" "$rawdir" >> "$EXTRA_MISSING.tmp"
  elif [ "$short_count" -eq 0 ]; then
    status=longread_only
  elif [ "$long_count" -eq 0 ]; then
    status=shortRNA_only
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$rawdir" "$short_count" "$long_count" "$status" "$note" >> "$EXTRA_SUMMARY.tmp"
}

add_t2t genome_A17  "$ROOT/data_A17"  "mapped_from_project_T2T_data_A17"
add_t2t genome_R108 "$ROOT/data_R108" "mapped_from_project_T2T_data_R108"
add_t2t genome_Mpo  "$ROOT/data_Mpo"  "mapped_from_project_T2T_data_Mpo"
add_t2t genome_Msa1 "$ROOT/data_Msa"  "mapped_from_project_T2T_data_Msa_shared_with_Msa2"
add_t2t genome_Msa2 "$ROOT/data_Msa"  "mapped_from_project_T2T_data_Msa_shared_with_Msa1"
add_t2t genome_474a "$ROOT/data_474"  "mapped_from_project_T2T_data_474_shared_with_474b"
add_t2t genome_474b "$ROOT/data_474"  "mapped_from_project_T2T_data_474_shared_with_474a"

mv "$EXTRA_MANIFEST.tmp" "$EXTRA_MANIFEST"
mv "$EXTRA_SUMMARY.tmp" "$EXTRA_SUMMARY"
mv "$EXTRA_MISSING.tmp" "$EXTRA_MISSING"

# Rebuild combined manifest by taking the pan_genome_new manifest rows plus project_T2T extra rows.
main="$TARGET/manifest.tsv"
base_tmp="$main.base_no_project_T2T.tmp"
if [ -f "$main" ]; then
  awk 'NR==1 || $3 !~ /project_T2T__/' "$main" > "$base_tmp"
else
  printf 'sample\tevidence_type\tlink_path\ttarget_path\ttarget_size\ttarget_mtime\traw_species_dir\tnote\n' > "$base_tmp"
fi
tail -n +2 "$EXTRA_MANIFEST" >> "$base_tmp"
mv "$base_tmp" "$main"

# Build combined summary from symlinks by sample/evidence type.
summary="$TARGET/summary.tsv"
printf 'sample\tillumina_links\tiso_links\tstatus\tnote\n' > "$summary.tmp"
for d in "$TARGET"/genome_*; do
  [ -d "$d" ] || continue
  sample=$(basename "$d")
  illum=$(find "$d/illumina-seq" -maxdepth 1 -type l 2>/dev/null | wc -l)
  iso=$(find "$d/iso-seq" -maxdepth 1 -type l 2>/dev/null | wc -l)
  if [ "$illum" -eq 0 ] && [ "$iso" -eq 0 ]; then status=missing; elif [ "$illum" -eq 0 ]; then status=iso_only; elif [ "$iso" -eq 0 ]; then status=illumina_only; else status=ok; fi
  note="combined_pan_genome_new_plus_project_T2T"
  printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$illum" "$iso" "$status" "$note" >> "$summary.tmp"
done
mv "$summary.tmp" "$summary"

echo "extra_manifest=$EXTRA_MANIFEST"
echo "extra_summary=$EXTRA_SUMMARY"
echo "extra_missing=$EXTRA_MISSING"
echo "combined_manifest=$main"
echo "combined_summary=$summary"
