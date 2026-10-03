# resequencing_population

重测序 SNP/InDel 鉴定（GATK）与群体结构（PCA、ADMIXTURE、进化树、π/Fst）

| 脚本 | 说明 |
|---|---|
| [`2.call_SNP_new/00_reference_index/run_00_reference_index.sh`](2.call_SNP_new/00_reference_index/run_00_reference_index.sh) | 运行：00 reference index（据文件名） |
| [`2.call_SNP_new/01_sample_list/run_01_make_sample_list.sh`](2.call_SNP_new/01_sample_list/run_01_make_sample_list.sh) | 运行：01 make sample list（据文件名） |
| [`2.call_SNP_new/07_filter_snp/make_yue_biallelic_snp_for_gwas.sh`](2.call_SNP_new/07_filter_snp/make_yue_biallelic_snp_for_gwas.sh) | 生成：yue biallelic snp for gwas（据文件名） |
| [`2.call_SNP_new/08_filter_indel/run_yue_filter_indel_by_interval.sh`](2.call_SNP_new/08_filter_indel/run_yue_filter_indel_by_interval.sh) | 运行：yue filter indel by interval（据文件名） |
| [`2.call_SNP_new/99_run_all_stepwise/run_all_stepwise.sh`](2.call_SNP_new/99_run_all_stepwise/run_all_stepwise.sh) | 运行：all stepwise（据文件名） |
| [`test_population/01.population_analysis/01.run_population_analysis.sh`](test_population/01.population_analysis/01.run_population_analysis.sh) | 运行：population analysis（据文件名） |
| [`test_population/01.population_analysis/07_pi_fst_K4/99_scripts/plot_pi_fst_network.R`](test_population/01.population_analysis/07_pi_fst_K4/99_scripts/plot_pi_fst_network.R) | 绘图：pi fst network（据文件名） |
| [`test_population/01.population_analysis/07_tree_K4_no_outgroup/99_scripts/plot_ml_nj_trees.R`](test_population/01.population_analysis/07_tree_K4_no_outgroup/99_scripts/plot_ml_nj_trees.R) | 绘图：ml nj trees（据文件名） |
| [`test_population/01.population_analysis/08_tree_K4_with_outgroup/99_scripts/plot_ml_nj_trees_with_outgroup.R`](test_population/01.population_analysis/08_tree_K4_with_outgroup/99_scripts/plot_ml_nj_trees_with_outgroup.R) | 绘图：ml nj trees with outgroup（据文件名） |
| [`test_population/01.population_analysis/09_C3_C2_special_check.extract_admixture_K4.R`](test_population/01.population_analysis/09_C3_C2_special_check.extract_admixture_K4.R) |  |
| [`test_population/01.population_analysis/09_C3_C2_special_check.pca_knn_474.R`](test_population/01.population_analysis/09_C3_C2_special_check.pca_knn_474.R) |  |
| [`test_population/01.population_analysis/09_C3_C2_special_check.R`](test_population/01.population_analysis/09_C3_C2_special_check.R) |  |
| [`test_population/01.population_analysis/10_C3_runaway_confirm.R`](test_population/01.population_analysis/10_C3_runaway_confirm.R) |  |
| [`test_population/01.population_analysis/10_C3_runaway_confirm.refine.R`](test_population/01.population_analysis/10_C3_runaway_confirm.refine.R) |  |
| [`test_population/01.population_analysis/11_revised_C3_to_C2/99_scripts/fix_revised_pca_legend.R`](test_population/01.population_analysis/11_revised_C3_to_C2/99_scripts/fix_revised_pca_legend.R) | 修正：revised pca legend（据文件名） |
| [`test_population/01.population_analysis/11_revised_C3_to_C2/99_scripts/redraw_revised_C3_to_C2_plain.R`](test_population/01.population_analysis/11_revised_C3_to_C2/99_scripts/redraw_revised_C3_to_C2_plain.R) |  |
| [`test_population/01.population_analysis/11_revised_C3_to_C2/99_scripts/revised_groups_pca_trees.R`](test_population/01.population_analysis/11_revised_C3_to_C2/99_scripts/revised_groups_pca_trees.R) |  |
| [`test_population/01.population_analysis/11_revised_C3_to_C2/99_scripts/summarize_plot_pi_fst.R`](test_population/01.population_analysis/11_revised_C3_to_C2/99_scripts/summarize_plot_pi_fst.R) | 汇总统计：plot pi fst（据文件名） |
| [`test_population/01.population_analysis/99_scripts/plot_admixture.R`](test_population/01.population_analysis/99_scripts/plot_admixture.R) | 绘图：admixture（据文件名） |
| [`test_population/01.population_analysis/99_scripts/plot_pca.R`](test_population/01.population_analysis/99_scripts/plot_pca.R) | 绘图：pca（据文件名） |
| [`test_population/01.population_analysis/99_scripts/plot_pca_k4_outgroup.R`](test_population/01.population_analysis/99_scripts/plot_pca_k4_outgroup.R) | 绘图：pca k4 outgroup（据文件名） |
| [`test_population/01.population_analysis/99_scripts/plot_pca_k4_outgroup_paperstyle.R`](test_population/01.population_analysis/99_scripts/plot_pca_k4_outgroup_paperstyle.R) | 绘图：pca k4 outgroup paperstyle（据文件名） |
| [`test_population/01.population_analysis/99_scripts/vcf_to_phy.py`](test_population/01.population_analysis/99_scripts/vcf_to_phy.py) |  |
| [`test_population/01.population_analysis/run_revised_C3_to_C2_all.sh`](test_population/01.population_analysis/run_revised_C3_to_C2_all.sh) | 运行：revised C3 to C2 all（据文件名） |
| [`test_population/01.population_analysis/run_tree_K4_no_outgroup.sh`](test_population/01.population_analysis/run_tree_K4_no_outgroup.sh) | 运行：tree K4 no outgroup（据文件名） |
| [`test_population/01.population_analysis/run_tree_K4_with_outgroup.sh`](test_population/01.population_analysis/run_tree_K4_with_outgroup.sh) | 运行：tree K4 with outgroup（据文件名） |
