#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project
RUN="$ROOT/N_5.pod_coil/01_SHP_homology_and_coil_cosegregation"
DATA18="$ROOT/N_4.pod_spiny/00_data"
DATA="$ROOT/N_5.pod_coil/00_data"
ENV=path/to/home/anaconda3/envs/biosofeware/bin
BCF=path/to/home/anaconda3/envs/panpop/bin/bcftools
PYTHON=python3

mkdir -p "$RUN"/{scripts,proteomes,blast,db,rbh,alignments,trees,population,tables,logs,tmp}
cp "$0" "$RUN/scripts/run_shp_analysis.sh" 2>/dev/null || true

cat > "$RUN/species.list" <<'EOF'
genome_395
genome_410
genome_436
genome_454
genome_457
genome_461
genome_468
genome_472
genome_474a
genome_482
genome_M22
genome_M46
genome_Mar
genome_Mpo
genome_Mru
genome_Msa
genome_R108
genome_ZM4
EOF

echo "[$(date '+%F %T')] build proteomes" | tee "$RUN/logs/run.log"
while read -r sp; do
  genome="$DATA18/1.reference_genome/${sp}.fa"
  gff="$DATA18/4.reference_anno/${sp}.gff"
  if [[ ! -s "$genome" ]]; then
    genome="$ROOT/N_1.coding_gene_anno/01_genome_versions/unmask/${sp}.unmasked.fa"
  fi
  if [[ -s "$RUN/proteomes/${sp}.pep" && -s "$RUN/proteomes/${sp}.cds" ]]; then
    continue
  fi
  "$ENV/gffread" "$gff" -g "$genome" -y "$RUN/proteomes/${sp}.pep" -x "$RUN/proteomes/${sp}.cds" \
    >"$RUN/logs/${sp}.gffread.out" 2>"$RUN/logs/${sp}.gffread.err"
done < "$RUN/species.list"

# gffread may encode terminal stops as '.' for a small number of models.
# BLAST accepts X as an unknown amino acid but rejects '.'.
while read -r sp; do
  perl -pi -e 'unless(/^>/){s/[^A-Za-z\r\n]/X/g}' "$RUN/proteomes/${sp}.pep"
done < "$RUN/species.list"

$PYTHON "$RUN/scripts/fasta_tools.py" --fasta "$RUN/proteomes/genome_Msa.pep" \
  --ids Chr14120.1,Chr24229.1 --output "$RUN/rbh/Msa_queries.pep"
"$ENV/makeblastdb" -in "$RUN/proteomes/genome_Msa.pep" -dbtype prot -out "$RUN/db/genome_Msa" >/dev/null

while read -r sp; do
  "$ENV/makeblastdb" -in "$RUN/proteomes/${sp}.pep" -dbtype prot -out "$RUN/db/${sp}" >/dev/null
  "$ENV/blastp" -query "$RUN/rbh/Msa_queries.pep" -db "$RUN/db/${sp}" -max_target_seqs 5 \
    -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' \
    -out "$RUN/blast/${sp}.forward.tsv"
done < "$RUN/species.list"

$PYTHON "$RUN/scripts/build_rbh_tables.py" --work "$RUN" --species "$RUN/species.list" --stage candidates
while read -r sp; do
  "$ENV/blastp" -query "$RUN/rbh/${sp}.candidate.pep" -db "$RUN/db/genome_Msa" -max_target_seqs 5 \
    -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' \
    -out "$RUN/blast/${sp}.reverse.tsv"
done < "$RUN/species.list"
$PYTHON "$RUN/scripts/build_rbh_tables.py" --work "$RUN" --species "$RUN/species.list" --stage final

for gene in Chr14120.1 Chr24229.1; do
  "$ENV/mafft" --auto "$RUN/rbh/strict_${gene}.pep" > "$RUN/alignments/strict_${gene}.pep.aln.fa" 2> "$RUN/logs/${gene}.pep.mafft.log"
  "$ENV/mafft" --auto "$RUN/rbh/strict_${gene}.cds" > "$RUN/alignments/strict_${gene}.cds.aln.fa" 2> "$RUN/logs/${gene}.cds.mafft.log"
  $PYTHON "$RUN/scripts/association.py" cross --alignment "$RUN/alignments/strict_${gene}.pep.aln.fa" \
    --phenotype "$DATA18/5.phenotye_info/phenotype.tsv" --output "$RUN/tables/${gene}.protein_site_coil_association.tsv"
  $PYTHON "$RUN/scripts/association.py" cross --alignment "$RUN/alignments/strict_${gene}.cds.aln.fa" \
    --phenotype "$DATA18/5.phenotye_info/phenotype.tsv" --output "$RUN/tables/${gene}.cds_site_coil_association.tsv"
  (cd "$RUN/trees" && "$ENV/iqtree2" -s "$RUN/alignments/strict_${gene}.pep.aln.fa" -m MFP -B 1000 -alrt 1000 \
     -T 4 --prefix "strict_${gene}" >"$RUN/logs/${gene}.iqtree.log" 2>&1)
done

# Reference-based SHP clade validation with curated Arabidopsis proteins.
ARAPORT="$ROOT/16.T2T_ref_function_anno/output/blast_tair/tair_model/Araport11_pep_20220914"
cat > "$RUN/tmp/arabidopsis_mads.ids" <<'EOF'
AT3G58780.1
AT2G42830.1
AT4G18960.1
AT4G09960.1
AT5G60910.1
AT1G69120.1
AT1G24260.1
EOF
$PYTHON "$RUN/scripts/fasta_tools.py" --fasta "$ARAPORT" --ids "$RUN/tmp/arabidopsis_mads.ids" --output "$RUN/rbh/Arabidopsis_MADS_refs.pep"
cat "$RUN/rbh/Msa_queries.pep" "$RUN/rbh/Arabidopsis_MADS_refs.pep" > "$RUN/rbh/SHP_reference_tree.pep"
"$ENV/mafft" --auto "$RUN/rbh/SHP_reference_tree.pep" > "$RUN/alignments/SHP_reference_tree.pep.aln.fa" 2> "$RUN/logs/SHP_reference.mafft.log"
(cd "$RUN/trees" && "$ENV/iqtree2" -s "$RUN/alignments/SHP_reference_tree.pep.aln.fa" -m MFP -B 1000 -alrt 1000 \
  -T 4 --prefix SHP_reference >"$RUN/logs/SHP_reference.iqtree.log" 2>&1)

# Population SNP tests in gene bodies plus 2 kb flanks.
VCF="$DATA/SV_SNP/DP_6-85_miss_0.2.all_samples.biallelic.SNP.vcf.gz"
PHENO="$DATA/144_phenotype/medicago144_pod_spiral_spine_traits.tsv"
"$BCF" query -l "$VCF" > "$RUN/population/vcf.samples"
for spec in 'Chr14120 Chr3:53194386-53205316' 'Chr24229 Chr4:91171307-91180034'; do
  set -- $spec
  gene=$1
  region=$2
  "$BCF" query -r "$region" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$VCF" > "$RUN/population/${gene}.SNP.genotypes.tsv"
  $PYTHON "$RUN/scripts/association.py" population --query "$RUN/population/${gene}.SNP.genotypes.tsv" \
    --samples "$RUN/population/vcf.samples" --phenotype "$PHENO" --gene "$gene" \
    --output "$RUN/tables/${gene}.population_SNP_coil_association.tsv"
done

{
  echo -e 'item\tvalue'
  echo -e "analysis_date\t$(date '+%F %T')"
  echo -e "reference\tgenome_Msa"
  echo -e "species_count\t$(wc -l < "$RUN/species.list")"
  echo -e "population_phenotype_rows\t$(( $(wc -l < "$PHENO") - 1 ))"
  echo -e "population_vcf_samples\t$(wc -l < "$RUN/population/vcf.samples")"
  echo -e "Chr14120_coordinates\tChr3:53196386-53203316(-)"
  echo -e "Chr24229_coordinates\tChr4:91173307-91178034(+)"
} > "$RUN/tables/run_summary.tsv"
touch "$RUN/SHP_analysis.done"
echo "[$(date '+%F %T')] done" | tee -a "$RUN/logs/run.log"
