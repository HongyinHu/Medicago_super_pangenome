# assembly_quality

组装质量评估（BUSCO、Merqury、reads 回比）

| 脚本 | 说明 |
|---|---|
| [`my_run/s1.run_busco_assembly.py`](my_run/s1.run_busco_assembly.py) | 运行busco评估基因组组装和蛋白注释 |
| [`my_run/s2.run_bwa_assembly.py`](my_run/s2.run_bwa_assembly.py) | 评估基因组组装的完整性 使用bwa-mem将二代基因组reads比对组装序列上 |
| [`my_run/s3.run_merqury.py`](my_run/s3.run_merqury.py) | 使用merqury评估基因组组织质量 |
| [`output/00_protein_busco_summary/summarize_protein_busco.py`](output/00_protein_busco_summary/summarize_protein_busco.py) | Collect one selected protein-mode BUSCO result per top-level genome directory. |
| [`output/X_genome_M46_2/scripts/assembly_stats.py`](output/X_genome_M46_2/scripts/assembly_stats.py) |  |
| [`output/X_genome_M46_2/scripts/finalize_report.py`](output/X_genome_M46_2/scripts/finalize_report.py) | 收尾汇总：report（据文件名） |
| [`output/X_genome_M46_2/scripts/run_busco.sh`](output/X_genome_M46_2/scripts/run_busco.sh) | 运行：busco（据文件名） |
| [`output/X_genome_M46_2/scripts/run_merqury.sh`](output/X_genome_M46_2/scripts/run_merqury.sh) | 运行：merqury（据文件名） |
| [`output/X_genome_M46_2/scripts/run_readmap.sh`](output/X_genome_M46_2/scripts/run_readmap.sh) | 运行：readmap（据文件名） |
