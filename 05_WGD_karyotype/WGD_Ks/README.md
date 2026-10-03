# WGD_Ks

Ks distributions of paralogs with WGDI to identify WGD events

| Script | Description |
|---|---|
| [`my_run/h1.select_chr.py`](my_run/h1.select_chr.py) | Remove non-essential chromosome fragments from WGDI input files |
| [`my_run/h2.rename_gff_geneid.py`](my_run/h2.rename_gff_geneid.py) | Rename gene IDs in GFF files |
| [`my_run/h3.filter_block_inf.py`](my_run/h3.filter_block_inf.py) | Filter Mpr_Mpr_block_information.csv by block size |
| [`my_run/s1.get_gff_lens.py`](my_run/s1.get_gff_lens.py) | Convert GFF to the WGDI format |
| [`my_run/s2.run_blastp.py`](my_run/s2.run_blastp.py) | BLASTP of protein sequences |
| [`my_run/s3.run_wgdi_same_sample.py`](my_run/s3.run_wgdi_same_sample.py) | Write WGDI configuration files and run each WGDI step |
| [`my_run/s4.run_wgdi_ksfigure.py`](my_run/s4.run_wgdi_ksfigure.py) | Write WGDI configuration files and run each WGDI step |
| [`my_run/s5.run_mcscan_cds_pep.py`](my_run/s5.run_mcscan_cds_pep.py) | Plot synteny with MCScan (JCVI) |
| [`my_run/s6.run_mcscan_ortholog.py`](my_run/s6.run_mcscan_ortholog.py) | Compute pairwise synteny between species |
| [`my_run/s7.run_mcscan_plot.py`](my_run/s7.run_mcscan_plot.py) | Compute pairwise synteny between species |
| [`output/ks_figure/total_ksfigure.conf`](output/ks_figure/total_ksfigure.conf) | Configuration file |
| [`output/ks_figure/total_ksfigure.M46_2.conf`](output/ks_figure/total_ksfigure.M46_2.conf) | Configuration file |
| [`output/N_genome_410_0/total.conf`](output/N_genome_410_0/total.conf) | Configuration file |
