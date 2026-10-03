#!/usr/bin/env bash
set -euo pipefail
WORK_DIR=path/to/project/7.function_anno_genes/function_annotation_stats
OUTPUT_DIR=$WORK_DIR/output
SPECIES_PEP_DIR=path/to/project/7.function_anno_genes/data/species_pep
SAMPLE=genome_M46_2_new
SOURCE_PEP=path/to/projects_all/pan_genome/X_genome_M46_2/6.genome_annotation_gene_predict/output/4.evm_combind_3_finally/finally.genome.anno.pep
CONDA=path/to/home/anaconda3/bin/conda
PFAM_SCAN=path/to/home/anaconda3/envs/pfam/bin/pfam_scan.pl
INTERPROSCAN=path/to/home/sofeware/my_interproscan/interproscan-5.30-69.0/interproscan.sh
EGGNOG_DB_DIR=path/to/home/db/eggnog
PFAM_DB=path/to/data/database/pfam
NR_DB=path/to/project/7.function_anno_genes/data/nr.gz.dmnd
SWISS_DB=path/to/project/7.function_anno_genes/data/uniprot_sprot.fasta.gz.dmnd
TREMBL_DB=path/to/project/7.function_anno_genes/data/uniprot_trembl.fasta.dmnd
KOG_DB=path/to/data/database/KOG_database/KOG_dataset.fa.dmnd
CMD_DIR=$WORK_DIR/02_cmds
LOG_DIR=$WORK_DIR/01_logs
LOCK_DIR=$WORK_DIR/04_locks
DONE_DIR=$WORK_DIR/05_done
mkdir -p "$SPECIES_PEP_DIR" "$CMD_DIR" "$LOG_DIR" "$LOCK_DIR" "$DONE_DIR"
for d in pfam_anno InterproScan_anno NR_anno SwissPort_anno KOG_anno TrEMBL_anno eggNOG_anno; do
  mkdir -p "$OUTPUT_DIR/$d/$SAMPLE"
done
ln -sfn "$SOURCE_PEP" "$SPECIES_PEP_DIR/${SAMPLE}.source.pep"
awk '
  /^>/ {
    if (seq!="") { gsub(/[^A-Z]/,"X",seq); print ">"id; print seq }
    id=substr($0,2); sub(/[[:space:]].*$/, "", id); seq=""; next
  }
  { seq=seq toupper($0) }
  END { if (seq!="") { gsub(/[^A-Z]/,"X",seq); print ">"id; print seq } }
' "$SOURCE_PEP" > "$OUTPUT_DIR/pfam_anno/$SAMPLE/$SAMPLE.pep.format"
for d in InterproScan_anno NR_anno SwissPort_anno KOG_anno TrEMBL_anno eggNOG_anno; do
  ln -sfn "$OUTPUT_DIR/pfam_anno/$SAMPLE/$SAMPLE.pep.format" "$OUTPUT_DIR/$d/$SAMPLE/$SAMPLE.pep.format"
done
rm -f "$CMD_DIR/m46_2_new_pfam.cmd" "$CMD_DIR/m46_2_new_interproscan.cmd" "$CMD_DIR/m46_2_new_diamond.cmd" "$CMD_DIR/m46_2_new_eggnog.cmd"
rm -f "$DONE_DIR"/{pfam,interproscan,NR,SwissPort,KOG,TrEMBL,eggNOG}/$SAMPLE.done
cat > "$CMD_DIR/m46_2_new_pfam.cmd" <<EOF
$WORK_DIR/run_function_task_with_lock.sh $DONE_DIR/pfam/$SAMPLE.done $LOCK_DIR/pfam/$SAMPLE.lock bash -lc 'cd $OUTPUT_DIR/pfam_anno && rm -f $SAMPLE/$SAMPLE.pep.format.pfam.out && $CONDA run -n pfam $PFAM_SCAN -fasta $SAMPLE/$SAMPLE.pep.format -dir $PFAM_DB -outfile $SAMPLE/$SAMPLE.pep.format.pfam.out -cpu 30'
EOF
cat > "$CMD_DIR/m46_2_new_interproscan.cmd" <<EOF
$WORK_DIR/run_function_task_with_lock.sh $DONE_DIR/interproscan/$SAMPLE.done $LOCK_DIR/interproscan/$SAMPLE.lock bash -lc 'cd $OUTPUT_DIR/InterproScan_anno && rm -f $SAMPLE/$SAMPLE.pep.format.interproscan.out && $CONDA run -n java $INTERPROSCAN -f tsv -i $SAMPLE/$SAMPLE.pep.format -cpu 30 --highmem -o $SAMPLE/$SAMPLE.pep.format.interproscan.out -iprlookup -goterms -pa -td $SAMPLE/temp'
EOF
for spec in "NR:$NR_DB" "SwissPort:$SWISS_DB" "KOG:$KOG_DB" "TrEMBL:$TREMBL_DB"; do
  name=${spec%%:*}; db=${spec#*:}
  cat >> "$CMD_DIR/m46_2_new_diamond.cmd" <<EOF
$WORK_DIR/run_function_task_with_lock.sh $DONE_DIR/$name/$SAMPLE.done $LOCK_DIR/$name/$SAMPLE.lock bash -lc 'cd $OUTPUT_DIR/${name}_anno && rm -f $SAMPLE/$SAMPLE.pep.format_blastp.tab && $CONDA run -n biosofeware diamond blastp --threads 30 --db $db --query $SAMPLE/$SAMPLE.pep.format --out $SAMPLE/$SAMPLE.pep.format_blastp.tab --outfmt 6 --sensitive --evalue 1e-5 --quiet'
EOF
done
cat > "$CMD_DIR/m46_2_new_eggnog.cmd" <<EOF
$WORK_DIR/run_function_task_with_lock.sh $DONE_DIR/eggNOG/$SAMPLE.done $LOCK_DIR/eggNOG/$SAMPLE.lock bash -lc 'cd $OUTPUT_DIR/eggNOG_anno && rm -f $SAMPLE/$SAMPLE.emapper.annotations && $CONDA run -n eggnog_env emapper.py -i $SAMPLE/$SAMPLE.pep.format --itype proteins -m diamond --data_dir $EGGNOG_DB_DIR --cpu 30 --override --output $SAMPLE --output_dir $SAMPLE'
EOF
start_queue() {
  local node=$1 jobs=$2 name=$3 file=$4
  ssh "$node" "cd $WORK_DIR && nohup bash run_command_queue_no_parallel.sh $jobs $file $LOG_DIR/${name}.${node}.stdout $LOG_DIR/${name}.${node}.stderr >$LOG_DIR/${name}.${node}.queue.log 2>&1 & echo \$! > $WORK_DIR/${name}.${node}.pid"
}
start_queue lz31 1 m46_2_new_pfam "$CMD_DIR/m46_2_new_pfam.cmd"
start_queue lz31 1 m46_2_new_interproscan "$CMD_DIR/m46_2_new_interproscan.cmd"
start_queue lz31 2 m46_2_new_diamond "$CMD_DIR/m46_2_new_diamond.cmd"
start_queue lz31 1 m46_2_new_eggnog "$CMD_DIR/m46_2_new_eggnog.cmd"
echo "prepared and started $SAMPLE"
echo "input_records=$(grep -c '^>' "$OUTPUT_DIR/pfam_anno/$SAMPLE/$SAMPLE.pep.format")"
