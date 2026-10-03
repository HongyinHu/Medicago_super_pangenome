# pod_coil

荚果螺旋：SHP 同源与共分离、SNP/InDel/SV GWAS 及模型敏感性

| 脚本 | 说明 |
|---|---|
| [`01_SHP_homology_and_coil_cosegregation/scripts/annotate_gene_snps.py`](01_SHP_homology_and_coil_cosegregation/scripts/annotate_gene_snps.py) | 注释：gene snps（据文件名） |
| [`01_SHP_homology_and_coil_cosegregation/scripts/association.py`](01_SHP_homology_and_coil_cosegregation/scripts/association.py) |  |
| [`01_SHP_homology_and_coil_cosegregation/scripts/build_rbh_tables.py`](01_SHP_homology_and_coil_cosegregation/scripts/build_rbh_tables.py) | 构建：rbh tables（据文件名） |
| [`01_SHP_homology_and_coil_cosegregation/scripts/count_shp_paralogs.sh`](01_SHP_homology_and_coil_cosegregation/scripts/count_shp_paralogs.sh) | 计数：shp paralogs（据文件名） |
| [`01_SHP_homology_and_coil_cosegregation/scripts/fasta_tools.py`](01_SHP_homology_and_coil_cosegregation/scripts/fasta_tools.py) |  |
| [`01_SHP_homology_and_coil_cosegregation/scripts/run_hq_stats.sh`](01_SHP_homology_and_coil_cosegregation/scripts/run_hq_stats.sh) | 运行：hq stats（据文件名） |
| [`01_SHP_homology_and_coil_cosegregation/scripts/run_population_only.sh`](01_SHP_homology_and_coil_cosegregation/scripts/run_population_only.sh) | 运行：population only（据文件名） |
| [`01_SHP_homology_and_coil_cosegregation/scripts/run_reference_tree_only.sh`](01_SHP_homology_and_coil_cosegregation/scripts/run_reference_tree_only.sh) | 运行：reference tree only（据文件名） |
| [`01_SHP_homology_and_coil_cosegregation/scripts/run_shp_analysis.sh`](01_SHP_homology_and_coil_cosegregation/scripts/run_shp_analysis.sh) | 运行：shp analysis（据文件名） |
| [`01_SHP_homology_and_coil_cosegregation/scripts/run_snp_annotation.sh`](01_SHP_homology_and_coil_cosegregation/scripts/run_snp_annotation.sh) | 运行：snp annotation（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/build_gemma_inputs.py`](02_SNP_INDEL_SV_GWAS_20260713/scripts/build_gemma_inputs.py) | 构建：gemma inputs（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/diagnose_gwas_calibration.py`](02_SNP_INDEL_SV_GWAS_20260713/scripts/diagnose_gwas_calibration.py) | 诊断：gwas calibration（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/prepare_strict_coil.py`](02_SNP_INDEL_SV_GWAS_20260713/scripts/prepare_strict_coil.py) | 准备输入：strict coil（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/reorder_kinship_for_sv.py`](02_SNP_INDEL_SV_GWAS_20260713/scripts/reorder_kinship_for_sv.py) | 重排：kinship for sv（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/run_gwas_kinship_only_sensitivity.sh`](02_SNP_INDEL_SV_GWAS_20260713/scripts/run_gwas_kinship_only_sensitivity.sh) | 运行：gwas kinship only sensitivity（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/run_snp_indel_gwas.sh`](02_SNP_INDEL_SV_GWAS_20260713/scripts/run_snp_indel_gwas.sh) | 运行：snp indel gwas（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/run_sv_followup.sh`](02_SNP_INDEL_SV_GWAS_20260713/scripts/run_sv_followup.sh) | 运行：sv followup（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/run_sv_kinship_sensitivity.sh`](02_SNP_INDEL_SV_GWAS_20260713/scripts/run_sv_kinship_sensitivity.sh) | 运行：sv kinship sensitivity（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/summarize_gwas.py`](02_SNP_INDEL_SV_GWAS_20260713/scripts/summarize_gwas.py) | 汇总统计：gwas（据文件名） |
| [`02_SNP_INDEL_SV_GWAS_20260713/scripts/validate_existing_sv_genotypes.py`](02_SNP_INDEL_SV_GWAS_20260713/scripts/validate_existing_sv_genotypes.py) | 验证：existing sv genotypes（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/analyze_sv_fisher_logistic.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/analyze_sv_fisher_logistic.py) | ALT-carrier Fisher screen and plotting for the supplementary SV analysis. |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/build_gwas_diagnostic_report.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/build_gwas_diagnostic_report.py) | 构建：gwas diagnostic report（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/compare_gwas_models.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/compare_gwas_models.py) | 比较：gwas models（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/compare_maf_sensitivity.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/compare_maf_sensitivity.py) | 比较：maf sensitivity（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/diagnose_model_calibration.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/diagnose_model_calibration.py) | Audit GWAS calibration without changing association results. |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/make_model_metadata.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/make_model_metadata.py) | 生成：model metadata（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/make_structure_diagnostics.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/make_structure_diagnostics.py) | Create PCA and SNP-GRM diagnostics for the 144-sample coil GWAS. |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/prepare_coil_ge1.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/prepare_coil_ge1.py) | 准备输入：coil ge1（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/prepare_modelA_maf005.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/prepare_modelA_maf005.sh) | 准备输入：modelA maf005（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_combined_finalize.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_combined_finalize.sh) | 运行：combined finalize（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_gmmat_gwas.R`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_gmmat_gwas.R) | 运行：gmmat gwas（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_gmmat_gwas_maf.R`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_gmmat_gwas_maf.R) | 运行：gmmat gwas maf（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_marker_gwas.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_marker_gwas.sh) | 运行：marker gwas（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_model_sensitivity_finalize.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_model_sensitivity_finalize.sh) | 运行：model sensitivity finalize（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_model_sensitivity_marker.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_model_sensitivity_marker.sh) | 运行：model sensitivity marker（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_model_sensitivity_sv.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_model_sensitivity_sv.sh) | 运行：model sensitivity sv（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_diagnostics.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_diagnostics.sh) | 运行：modelA maf diagnostics（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_finalize.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_finalize.sh) | 运行：modelA maf finalize（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_marker.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_marker.sh) | 运行：modelA maf marker（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_report.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_report.sh) | 运行：modelA maf report（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_sv.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_modelA_maf_sv.sh) | 运行：modelA maf sv（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_qh3_sv_genotype.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_qh3_sv_genotype.sh) | 运行：qh3 sv genotype（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_snp_indel_finalize.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_snp_indel_finalize.sh) | 运行：snp indel finalize（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_structure_144.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_structure_144.sh) | 运行：structure 144（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_sv_144.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_sv_144.sh) | 运行：sv 144（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_sv_fisher_logistic_supplement.sh`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/run_sv_fisher_logistic_supplement.sh) | 运行：sv fisher logistic supplement（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/summarize_gmmat.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/summarize_gmmat.py) | 汇总统计：gmmat（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/summarize_maf_cutoff.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/summarize_maf_cutoff.py) | 汇总统计：maf cutoff（据文件名） |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/summarize_model_sensitivity.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/summarize_model_sensitivity.py) | Compare Model A/B/C GWAS calibration and peak stability. |
| [`03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/write_model_sensitivity_report.py`](03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714/scripts/write_model_sensitivity_report.py) | Write an evidence-based recommendation from the completed A/B/C audit. |
| [`08_SHP_multiallelic_candidate_20260717/scripts/build_candidate_report.py`](08_SHP_multiallelic_candidate_20260717/scripts/build_candidate_report.py) | Join audited functional alleles to low-frequency GMMAT score outputs. |
| [`08_SHP_multiallelic_candidate_20260717/scripts/candidate_variant_tools.py`](08_SHP_multiallelic_candidate_20260717/scripts/candidate_variant_tools.py) | Small helpers for audited SHP/AG-like functional alleles. |
| [`08_SHP_multiallelic_candidate_20260717/scripts/plot_snp_unadjusted_vs_adjusted.R`](08_SHP_multiallelic_candidate_20260717/scripts/plot_snp_unadjusted_vs_adjusted.R) | 绘图：snp unadjusted vs adjusted（据文件名） |
| [`08_SHP_multiallelic_candidate_20260717/scripts/run.sh`](08_SHP_multiallelic_candidate_20260717/scripts/run.sh) | 运行（据文件名） |
| [`08_SHP_multiallelic_candidate_20260717/scripts/run_candidate_gmmat.R`](08_SHP_multiallelic_candidate_20260717/scripts/run_candidate_gmmat.R) | 运行：candidate gmmat（据文件名） |
| [`08_SHP_multiallelic_candidate_20260717/scripts/run_snp_unadjusted_vs_adjusted_plots.sh`](08_SHP_multiallelic_candidate_20260717/scripts/run_snp_unadjusted_vs_adjusted_plots.sh) | 运行：snp unadjusted vs adjusted plots（据文件名） |
| [`08_SHP_multiallelic_candidate_20260717/scripts/run_uncorrected_snp_assoc.sh`](08_SHP_multiallelic_candidate_20260717/scripts/run_uncorrected_snp_assoc.sh) | 运行：uncorrected snp assoc（据文件名） |
