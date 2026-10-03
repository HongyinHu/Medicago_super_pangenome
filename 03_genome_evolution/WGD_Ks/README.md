# WGD_Ks

WGDI 同源基因 Ks 分布与 WGD 识别

| 脚本 | 说明 |
|---|---|
| [`my_run/h1.select_chr.py`](my_run/h1.select_chr.py) | 为WGDI输入文件过滤部分非必要染色体片段 |
| [`my_run/h2.rename_gff_geneid.py`](my_run/h2.rename_gff_geneid.py) | 重新命名gff文件中的geneid名称 |
| [`my_run/h3.filter_block_inf.py`](my_run/h3.filter_block_inf.py) | 过滤Mpr_Mpr_block_information.csv大小 |
| [`my_run/s1.get_gff_lens.py`](my_run/s1.get_gff_lens.py) | 获取wgdi要求的gff格式 |
| [`my_run/s2.run_blastp.py`](my_run/s2.run_blastp.py) | 获取蛋白序列的blastp比对结果 |
| [`my_run/s3.run_wgdi_same_sample.py`](my_run/s3.run_wgdi_same_sample.py) | 书写wgdi配置文件和运行各个命令 |
| [`my_run/s4.run_wgdi_ksfigure.py`](my_run/s4.run_wgdi_ksfigure.py) | 书写wgdi配置文件和运行各个命令 |
| [`my_run/s5.run_mcscan_cds_pep.py`](my_run/s5.run_mcscan_cds_pep.py) | 使用mcscan绘制共线性图 |
| [`my_run/s6.run_mcscan_ortholog.py`](my_run/s6.run_mcscan_ortholog.py) | 获取两两物种间共线性 |
| [`my_run/s7.run_mcscan_plot.py`](my_run/s7.run_mcscan_plot.py) | 获取两两物种间共线性 |
| [`output/ks_figure/total_ksfigure.conf`](output/ks_figure/total_ksfigure.conf) |  |
| [`output/ks_figure/total_ksfigure.M46_2.conf`](output/ks_figure/total_ksfigure.M46_2.conf) |  |
| [`output/N_genome_410_0/total.conf`](output/N_genome_410_0/total.conf) |  |
