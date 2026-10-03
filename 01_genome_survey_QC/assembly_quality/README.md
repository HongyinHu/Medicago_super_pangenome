# assembly_quality

Assembly quality assessment (BUSCO, Merqury, short-read mapping)

| Script | Description |
|---|---|
| [`my_run/s1.run_busco_assembly.py`](my_run/s1.run_busco_assembly.py) | Run BUSCO on assemblies and protein annotations |
| [`my_run/s2.run_bwa_assembly.py`](my_run/s2.run_bwa_assembly.py) | Assess assembly completeness by mapping short reads with BWA-MEM |
| [`my_run/s3.run_merqury.py`](my_run/s3.run_merqury.py) | Assess assembly quality with Merqury |
| [`output/00_protein_busco_summary/summarize_protein_busco.py`](output/00_protein_busco_summary/summarize_protein_busco.py) | Collect one selected protein-mode BUSCO result per top-level genome directory. |
| [`output/X_genome_M46_2/scripts/assembly_stats.py`](output/X_genome_M46_2/scripts/assembly_stats.py) |  |
| [`output/X_genome_M46_2/scripts/finalize_report.py`](output/X_genome_M46_2/scripts/finalize_report.py) | Finalize report |
| [`output/X_genome_M46_2/scripts/run_busco.sh`](output/X_genome_M46_2/scripts/run_busco.sh) | Run busco |
| [`output/X_genome_M46_2/scripts/run_merqury.sh`](output/X_genome_M46_2/scripts/run_merqury.sh) | Run merqury |
| [`output/X_genome_M46_2/scripts/run_readmap.sh`](output/X_genome_M46_2/scripts/run_readmap.sh) | Run readmap |
