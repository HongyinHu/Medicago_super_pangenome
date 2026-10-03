# gene_function_annotation

Functional annotation: DIAMOND (NR/Swiss-Prot/TrEMBL), eggNOG, InterProScan, Pfam, KOG

| Script | Description |
|---|---|
| [`function_annotation_stats/00_run_sh/00_run_command_queue_no_parallel.sh`](function_annotation_stats/00_run_sh/00_run_command_queue_no_parallel.sh) | Run command queue no parallel |
| [`function_annotation_stats/00_run_sh/00_run_function_task_with_lock.sh`](function_annotation_stats/00_run_sh/00_run_function_task_with_lock.sh) | Run function task with lock |
| [`function_annotation_stats/00_run_sh/03_add_three_species_function_annotation.sh`](function_annotation_stats/00_run_sh/03_add_three_species_function_annotation.sh) | Command: `WORK_DIR="${WORK_DIR:-function_annotation_stats}"` |
| [`function_annotation_stats/00_run_sh/04_add_M46_2_new_function_annotation.sh`](function_annotation_stats/00_run_sh/04_add_M46_2_new_function_annotation.sh) | Command: `WORK_DIR=function_annotation_stats` |
| [`function_annotation_stats/00_run_sh/05_finalize_M46_2_new_function_annotation.sh`](function_annotation_stats/00_run_sh/05_finalize_M46_2_new_function_annotation.sh) | Finalize M46 2 new function annotation |
| [`function_annotation_stats/00_run_sh/add_species.diamond.sh`](function_annotation_stats/00_run_sh/add_species.diamond.sh) | Command: `run_function_task_with_lock.sh genome_Mar.done genome_Mar.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/add_species.eggnog.sh`](function_annotation_stats/00_run_sh/add_species.eggnog.sh) | Command: `run_function_task_with_lock.sh genome_Mar.done genome_Mar.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/add_species.interproscan.sh`](function_annotation_stats/00_run_sh/add_species.interproscan.sh) | Command: `run_function_task_with_lock.sh genome_Mar.done genome_Mar.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/add_species.pfam.sh`](function_annotation_stats/00_run_sh/add_species.pfam.sh) | Command: `run_function_task_with_lock.sh genome_Mar.done genome_Mar.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/diamond.remaining.sh`](function_annotation_stats/00_run_sh/diamond.remaining.sh) | Command: `run_function_task_with_lock.sh genome_474.done genome_474.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/diamond.sh`](function_annotation_stats/00_run_sh/diamond.sh) |  |
| [`function_annotation_stats/00_run_sh/eggnog.remaining.sh`](function_annotation_stats/00_run_sh/eggnog.remaining.sh) | Command: `run_function_task_with_lock.sh genome_474.done genome_474.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/eggnog.sh`](function_annotation_stats/00_run_sh/eggnog.sh) |  |
| [`function_annotation_stats/00_run_sh/interproscan.remaining.sh`](function_annotation_stats/00_run_sh/interproscan.remaining.sh) | Command: `run_function_task_with_lock.sh genome_474.done genome_474.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/interproscan.sh`](function_annotation_stats/00_run_sh/interproscan.sh) |  |
| [`function_annotation_stats/00_run_sh/pfam.remaining.sh`](function_annotation_stats/00_run_sh/pfam.remaining.sh) | Command: `run_function_task_with_lock.sh genome_474.done genome_474.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/pfam.sh`](function_annotation_stats/00_run_sh/pfam.sh) |  |
| [`function_annotation_stats/02_cmds/m46_2_new_diamond.cmd`](function_annotation_stats/02_cmds/m46_2_new_diamond.cmd) |  |
| [`function_annotation_stats/02_cmds/m46_2_new_eggnog.cmd`](function_annotation_stats/02_cmds/m46_2_new_eggnog.cmd) |  |
| [`function_annotation_stats/02_cmds/m46_2_new_interproscan.cmd`](function_annotation_stats/02_cmds/m46_2_new_interproscan.cmd) |  |
| [`function_annotation_stats/02_cmds/m46_2_new_pfam.cmd`](function_annotation_stats/02_cmds/m46_2_new_pfam.cmd) |  |
| [`function_annotation_stats/add_three_species_function_annotation.sh`](function_annotation_stats/add_three_species_function_annotation.sh) | Command: `WORK_DIR="${WORK_DIR:-function_annotation_stats}"` |
| [`function_annotation_stats/extract_function_annotation_genelists.py`](function_annotation_stats/extract_function_annotation_genelists.py) | Build gene-list files and final coverage tables from functional annotation outputs. |
| [`function_annotation_stats/summarize_function_annotation_stats.py`](function_annotation_stats/summarize_function_annotation_stats.py) | Summarize functional annotation gene-list coverage for multiple genomes. |
| [`my_run/add_annotation_from_dat2.py`](my_run/add_annotation_from_dat2.py) |  |
| [`my_run/add_annotation_from_dat.py`](my_run/add_annotation_from_dat.py) |  |
| [`my_run/h1.change_flag_end.py`](my_run/h1.change_flag_end.py) | Fix stop-codon symbols in protein sequences |
| [`my_run/h1.diamond_mkdb.py`](my_run/h1.diamond_mkdb.py) | Build DIAMOND database and run alignment |
| [`my_run/h2.run_diamond.py`](my_run/h2.run_diamond.py) | Protein alignment with DIAMOND |
| [`my_run/h3.run_interproscan.py`](my_run/h3.run_interproscan.py) | GO annotation with InterProScan |
| [`my_run/h4.run_pfam_scan.py`](my_run/h4.run_pfam_scan.py) | Search the Pfam database with pfam_scan |
| [`my_run/h5.get_domin_interproscan.py`](my_run/h5.get_domin_interproscan.py) | Extract functional domain IDs per gene |
| [`my_run/h6.get_number_pfam.py`](my_run/h6.get_number_pfam.py) | Count genes annotated by Pfam |
| [`my_run/h7.get_genenum_database.py`](my_run/h7.get_genenum_database.py) | Count annotated genes per database |
| [`my_run/p1.filter_best_query.py`](my_run/p1.filter_best_query.py) | Keep the best hit per query from alignment results |
| [`my_run/p2.covert_uniport_dat_fasta.py`](my_run/p2.covert_uniport_dat_fasta.py) | Convert UniProt .dat to FASTA |
| [`my_run/s1.run_diamond_nr.py`](my_run/s1.run_diamond_nr.py) | Functional annotation of proteins with DIAMOND |
| [`my_run/s2.run_blastp.py`](my_run/s2.run_blastp.py) | Run blastp |
| [`my_run/s3.run_interproscan.py`](my_run/s3.run_interproscan.py) | Annotation with InterProScan |
| [`my_run/selenium_NR_easy.py`](my_run/selenium_NR_easy.py) | Fetch NCBI NR gene descriptions (Selenium) |
| [`my_run/selenium_uniprot_easy.py`](my_run/selenium_uniprot_easy.py) | Fetch UniProt gene descriptions (Selenium) |
