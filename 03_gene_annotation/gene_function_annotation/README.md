# gene_function_annotation

基因功能注释：DIAMOND（NR/SwissProt/TrEMBL）、eggNOG、InterProScan、Pfam、KOG

| 脚本 | 说明 |
|---|---|
| [`function_annotation_stats/00_run_sh/00_run_command_queue_no_parallel.sh`](function_annotation_stats/00_run_sh/00_run_command_queue_no_parallel.sh) | 运行：command queue no parallel（据文件名） |
| [`function_annotation_stats/00_run_sh/00_run_function_task_with_lock.sh`](function_annotation_stats/00_run_sh/00_run_function_task_with_lock.sh) | 运行：function task with lock（据文件名） |
| [`function_annotation_stats/00_run_sh/03_add_three_species_function_annotation.sh`](function_annotation_stats/00_run_sh/03_add_three_species_function_annotation.sh) | 命令：`WORK_DIR="${WORK_DIR:-function_annotation_stats}"` |
| [`function_annotation_stats/00_run_sh/04_add_M46_2_new_function_annotation.sh`](function_annotation_stats/00_run_sh/04_add_M46_2_new_function_annotation.sh) | 命令：`WORK_DIR=function_annotation_stats` |
| [`function_annotation_stats/00_run_sh/05_finalize_M46_2_new_function_annotation.sh`](function_annotation_stats/00_run_sh/05_finalize_M46_2_new_function_annotation.sh) | 收尾汇总：M46 2 new function annotation（据文件名） |
| [`function_annotation_stats/00_run_sh/add_species.diamond.sh`](function_annotation_stats/00_run_sh/add_species.diamond.sh) | 命令：`run_function_task_with_lock.sh genome_Mar.done genome_Mar.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/add_species.eggnog.sh`](function_annotation_stats/00_run_sh/add_species.eggnog.sh) | 命令：`run_function_task_with_lock.sh genome_Mar.done genome_Mar.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/add_species.interproscan.sh`](function_annotation_stats/00_run_sh/add_species.interproscan.sh) | 命令：`run_function_task_with_lock.sh genome_Mar.done genome_Mar.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/add_species.pfam.sh`](function_annotation_stats/00_run_sh/add_species.pfam.sh) | 命令：`run_function_task_with_lock.sh genome_Mar.done genome_Mar.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/diamond.remaining.sh`](function_annotation_stats/00_run_sh/diamond.remaining.sh) | 命令：`run_function_task_with_lock.sh genome_474.done genome_474.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/diamond.sh`](function_annotation_stats/00_run_sh/diamond.sh) |  |
| [`function_annotation_stats/00_run_sh/eggnog.remaining.sh`](function_annotation_stats/00_run_sh/eggnog.remaining.sh) | 命令：`run_function_task_with_lock.sh genome_474.done genome_474.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/eggnog.sh`](function_annotation_stats/00_run_sh/eggnog.sh) |  |
| [`function_annotation_stats/00_run_sh/interproscan.remaining.sh`](function_annotation_stats/00_run_sh/interproscan.remaining.sh) | 命令：`run_function_task_with_lock.sh genome_474.done genome_474.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/interproscan.sh`](function_annotation_stats/00_run_sh/interproscan.sh) |  |
| [`function_annotation_stats/00_run_sh/pfam.remaining.sh`](function_annotation_stats/00_run_sh/pfam.remaining.sh) | 命令：`run_function_task_with_lock.sh genome_474.done genome_474.lock bash -l…` |
| [`function_annotation_stats/00_run_sh/pfam.sh`](function_annotation_stats/00_run_sh/pfam.sh) |  |
| [`function_annotation_stats/02_cmds/m46_2_new_diamond.cmd`](function_annotation_stats/02_cmds/m46_2_new_diamond.cmd) |  |
| [`function_annotation_stats/02_cmds/m46_2_new_eggnog.cmd`](function_annotation_stats/02_cmds/m46_2_new_eggnog.cmd) |  |
| [`function_annotation_stats/02_cmds/m46_2_new_interproscan.cmd`](function_annotation_stats/02_cmds/m46_2_new_interproscan.cmd) |  |
| [`function_annotation_stats/02_cmds/m46_2_new_pfam.cmd`](function_annotation_stats/02_cmds/m46_2_new_pfam.cmd) |  |
| [`function_annotation_stats/add_three_species_function_annotation.sh`](function_annotation_stats/add_three_species_function_annotation.sh) | 命令：`WORK_DIR="${WORK_DIR:-function_annotation_stats}"` |
| [`function_annotation_stats/extract_function_annotation_genelists.py`](function_annotation_stats/extract_function_annotation_genelists.py) | Build gene-list files and final coverage tables from functional annotation outputs. |
| [`function_annotation_stats/summarize_function_annotation_stats.py`](function_annotation_stats/summarize_function_annotation_stats.py) | Summarize functional annotation gene-list coverage for multiple genomes. |
| [`my_run/add_annotation_from_dat2.py`](my_run/add_annotation_from_dat2.py) |  |
| [`my_run/add_annotation_from_dat.py`](my_run/add_annotation_from_dat.py) |  |
| [`my_run/h1.change_flag_end.py`](my_run/h1.change_flag_end.py) | 修改终止符号 |
| [`my_run/h1.diamond_mkdb.py`](my_run/h1.diamond_mkdb.py) | 使用diamond进行比对，数据库建库 |
| [`my_run/h2.run_diamond.py`](my_run/h2.run_diamond.py) | 使用diamond进行序列比对 |
| [`my_run/h3.run_interproscan.py`](my_run/h3.run_interproscan.py) | 使用interproscan进行GO注释 |
| [`my_run/h4.run_pfam_scan.py`](my_run/h4.run_pfam_scan.py) | 使用pfam_scan搜索pfam数据库 |
| [`my_run/h5.get_domin_interproscan.py`](my_run/h5.get_domin_interproscan.py) | 获取基因的功能结构域ID |
| [`my_run/h6.get_number_pfam.py`](my_run/h6.get_number_pfam.py) | 获取pfam注视到的基因数量 |
| [`my_run/h7.get_genenum_database.py`](my_run/h7.get_genenum_database.py) | 获取各个数据库中注释到基因数目 |
| [`my_run/p1.filter_best_query.py`](my_run/p1.filter_best_query.py) | 比对结果中筛选每个query的最佳subject |
| [`my_run/p2.covert_uniport_dat_fasta.py`](my_run/p2.covert_uniport_dat_fasta.py) | 将uniport dat格式转为fasta格式 |
| [`my_run/s1.run_diamond_nr.py`](my_run/s1.run_diamond_nr.py) | 使用diamond 对蛋白功能进行功能注释 |
| [`my_run/s2.run_blastp.py`](my_run/s2.run_blastp.py) | 运行：blastp（据文件名） |
| [`my_run/s3.run_interproscan.py`](my_run/s3.run_interproscan.py) | 使用interproscan进行注释 |
| [`my_run/selenium_NR_easy.py`](my_run/selenium_NR_easy.py) | 根据创建网址，爬取NCBI NR基因描述信息 |
| [`my_run/selenium_uniprot_easy.py`](my_run/selenium_uniprot_easy.py) | 根据创建网址，爬取uniprot基因描述信息 |
