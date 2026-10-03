# pod_spine

荚果刺：SV-表型共分离、短读长 SV、SNP/InDel/SV GWAS、候选基因 Chr23997 验证

| 脚本 | 说明 |
|---|---|
| [`03_candidate_gene_sv_cosegregation_Chr23997/scripts/analyze_Chr23997_sv_cosegregation.py`](03_candidate_gene_sv_cosegregation_Chr23997/scripts/analyze_Chr23997_sv_cosegregation.py) | 分析：Chr23997 sv cosegregation（据文件名） |
| [`03_candidate_gene_sv_cosegregation_Chr23997/scripts/make_igv_reports_Chr23997_sv.py`](03_candidate_gene_sv_cosegregation_Chr23997/scripts/make_igv_reports_Chr23997_sv.py) | 生成：igv reports Chr23997 sv（据文件名） |
| [`04_strict_event_flankQC_20260701/scripts/direction_free_cosegregation_from_flankQC.py`](04_strict_event_flankQC_20260701/scripts/direction_free_cosegregation_from_flankQC.py) | Direction-free cosegregation summary from per-sample flank-QC SV calls. |
| [`04_strict_event_flankQC_20260701/scripts/flank_qc_event_screen.py`](04_strict_event_flankQC_20260701/scripts/flank_qc_event_screen.py) |  |
| [`04_strict_event_flankQC_20260701/scripts/make_igv_reports_flankQC.py`](04_strict_event_flankQC_20260701/scripts/make_igv_reports_flankQC.py) | Build IGV HTML reports for strict flank-QC pod-spiny SV candidates. |
| [`04_strict_event_flankQC_20260701/scripts/run_flank_qc_test.sh`](04_strict_event_flankQC_20260701/scripts/run_flank_qc_test.sh) | 运行：flank qc test（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/check_chr23997_annotation.sh`](05_direction_free_sv_cosegregation_20260701/scripts/check_chr23997_annotation.sh) | 检查：chr23997 annotation（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/check_chr23997_annotation_hits.sh`](05_direction_free_sv_cosegregation_20260701/scripts/check_chr23997_annotation_hits.sh) | 检查：chr23997 annotation hits（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/find_chr23997_in_05.py`](05_direction_free_sv_cosegregation_20260701/scripts/find_chr23997_in_05.py) | 查找：chr23997 in 05（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/flank_qc_event_screen.py`](05_direction_free_sv_cosegregation_20260701/scripts/flank_qc_event_screen.py) |  |
| [`05_direction_free_sv_cosegregation_20260701/scripts/flank_qc_event_screen_dualref.py`](05_direction_free_sv_cosegregation_20260701/scripts/flank_qc_event_screen_dualref.py) |  |
| [`05_direction_free_sv_cosegregation_20260701/scripts/inspect_chr09113_function.sh`](05_direction_free_sv_cosegregation_20260701/scripts/inspect_chr09113_function.sh) | 检查：chr09113 function（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/inspect_chr09113_targeted.sh`](05_direction_free_sv_cosegregation_20260701/scripts/inspect_chr09113_targeted.sh) | 检查：chr09113 targeted（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/inspect_chr23997_caller_support.py`](05_direction_free_sv_cosegregation_20260701/scripts/inspect_chr23997_caller_support.py) | 检查：chr23997 caller support（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/inspect_chr23997_pansv.py`](05_direction_free_sv_cosegregation_20260701/scripts/inspect_chr23997_pansv.py) | 检查：chr23997 pansv（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/make_igv_reports_direction_free.py`](05_direction_free_sv_cosegregation_20260701/scripts/make_igv_reports_direction_free.py) | Build IGV HTML reports for strict flank-QC pod-spiny SV candidates. |
| [`05_direction_free_sv_cosegregation_20260701/scripts/rank_chr23997_tair_hits.sh`](05_direction_free_sv_cosegregation_20260701/scripts/rank_chr23997_tair_hits.sh) | 排序：chr23997 tair hits（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/run_make_igv_direction_free.sh`](05_direction_free_sv_cosegregation_20260701/scripts/run_make_igv_direction_free.sh) | 运行：make igv direction free（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/screen_direction_free_all_candidates.py`](05_direction_free_sv_cosegregation_20260701/scripts/screen_direction_free_all_candidates.py) | Direction-free SV/phenotype cosegregation screen for pod-spiny candidates. |
| [`05_direction_free_sv_cosegregation_20260701/scripts/summarize_05_ref_type.py`](05_direction_free_sv_cosegregation_20260701/scripts/summarize_05_ref_type.py) | 汇总统计：05 ref type（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/summarize_chr23997_records.py`](05_direction_free_sv_cosegregation_20260701/scripts/summarize_chr23997_records.py) | 汇总统计：chr23997 records（据文件名） |
| [`05_direction_free_sv_cosegregation_20260701/scripts/summarize_direction_free_05.py`](05_direction_free_sv_cosegregation_20260701/scripts/summarize_direction_free_05.py) | 汇总统计：direction free 05（据文件名） |
| [`07_trait_cosegregation_snakemake_20260701/config/config.yaml`](07_trait_cosegregation_snakemake_20260701/config/config.yaml) | Trait-guided SV cosegregation workflow configuration |
| [`07_trait_cosegregation_snakemake_20260701/config/project_medicago_config.yaml`](07_trait_cosegregation_snakemake_20260701/config/project_medicago_config.yaml) | Project-specific configuration for the Medicago pod spine analysis. |
| [`07_trait_cosegregation_snakemake_20260701/environment.yaml`](07_trait_cosegregation_snakemake_20260701/environment.yaml) |  |
| [`07_trait_cosegregation_snakemake_20260701/scripts/build_project_manifest_for_trait_workflow.sh`](07_trait_cosegregation_snakemake_20260701/scripts/build_project_manifest_for_trait_workflow.sh) | 构建：project manifest for trait workflow（据文件名） |
| [`07_trait_cosegregation_snakemake_20260701/scripts/draw_workflow_figure.py`](07_trait_cosegregation_snakemake_20260701/scripts/draw_workflow_figure.py) | 绘图：workflow figure（据文件名） |
| [`07_trait_cosegregation_snakemake_20260701/scripts/flank_coverage_qc.py`](07_trait_cosegregation_snakemake_20260701/scripts/flank_coverage_qc.py) |  |
| [`07_trait_cosegregation_snakemake_20260701/scripts/make_igv_manifest.py`](07_trait_cosegregation_snakemake_20260701/scripts/make_igv_manifest.py) | 生成：igv manifest（据文件名） |
| [`07_trait_cosegregation_snakemake_20260701/scripts/rank_candidates.py`](07_trait_cosegregation_snakemake_20260701/scripts/rank_candidates.py) | 排序：candidates（据文件名） |
| [`07_trait_cosegregation_snakemake_20260701/scripts/score_gene_interval_events.py`](07_trait_cosegregation_snakemake_20260701/scripts/score_gene_interval_events.py) | 打分：gene interval events（据文件名） |
| [`07_trait_cosegregation_snakemake_20260701/scripts/score_pav_cosegregation.py`](07_trait_cosegregation_snakemake_20260701/scripts/score_pav_cosegregation.py) | 打分：pav cosegregation（据文件名） |
| [`07_trait_cosegregation_snakemake_20260701/scripts/sv_utils.py`](07_trait_cosegregation_snakemake_20260701/scripts/sv_utils.py) |  |
| [`07_trait_cosegregation_snakemake_20260701/workflow/Snakefile`](07_trait_cosegregation_snakemake_20260701/workflow/Snakefile) |  |
| [`08_GWAS_pod_spine_SNP_SV_20260706/scripts/plot_snp_gwas_baseR.R`](08_GWAS_pod_spine_SNP_SV_20260706/scripts/plot_snp_gwas_baseR.R) | 绘图：snp gwas baseR（据文件名） |
| [`08_GWAS_pod_spine_SNP_SV_20260706/scripts/run_snp_gwas_plink2.sh`](08_GWAS_pod_spine_SNP_SV_20260706/scripts/run_snp_gwas_plink2.sh) | 运行：snp gwas plink2（据文件名） |
| [`08_GWAS_pod_spine_SNP_SV_20260706/scripts/run_snp_gwas_plink2_assoc_corrected.sh`](08_GWAS_pod_spine_SNP_SV_20260706/scripts/run_snp_gwas_plink2_assoc_corrected.sh) | 运行：snp gwas plink2 assoc corrected（据文件名） |
| [`08_GWAS_pod_spine_SNP_SV_20260706/scripts/run_snp_gwas_plink2_assoc_only.sh`](08_GWAS_pod_spine_SNP_SV_20260706/scripts/run_snp_gwas_plink2_assoc_only.sh) | 运行：snp gwas plink2 assoc only（据文件名） |
| [`08_GWAS_pod_spine_SNP_SV_20260706/scripts/summarize_plink2_gwas.py`](08_GWAS_pod_spine_SNP_SV_20260706/scripts/summarize_plink2_gwas.py) | 汇总统计：plink2 gwas（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/finalize_summary.sh`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/finalize_summary.sh) | 收尾汇总：summary（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/merge_bnd_support.py`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/merge_bnd_support.py) | 合并：bnd support（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/normalize_vcf.py`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/normalize_vcf.py) | 标准化：vcf（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_batch.sh`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_batch.sh) | 运行：batch（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_consensus_one.sh`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_consensus_one.sh) | 运行：consensus one（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_delly_one.sh`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_delly_one.sh) | 运行：delly one（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_manta_one.sh`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_manta_one.sh) | 运行：manta one（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_manta_smoove_batch.sh`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_manta_smoove_batch.sh) | 运行：manta smoove batch（据文件名） |
| [`09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_smoove_one.sh`](09_short_read_SV_Msa_20260706/medicago144_3caller_20260709/scripts/run_smoove_one.sh) | 运行：smoove one（据文件名） |
| [`09_short_read_SV_Msa_20260706/scripts/chr23997_depth_ratio.py`](09_short_read_SV_Msa_20260706/scripts/chr23997_depth_ratio.py) |  |
| [`09_short_read_SV_Msa_20260706/scripts/finalize_delly_DEL_genotypes.sh`](09_short_read_SV_Msa_20260706/scripts/finalize_delly_DEL_genotypes.sh) | 收尾汇总：delly DEL genotypes（据文件名） |
| [`09_short_read_SV_Msa_20260706/scripts/run_chr23997_depth_ratio_valid.py`](09_short_read_SV_Msa_20260706/scripts/run_chr23997_depth_ratio_valid.py) | 运行：chr23997 depth ratio valid（据文件名） |
| [`09_short_read_SV_Msa_20260706/scripts/run_delly_DEL_genotype_one.sh`](09_short_read_SV_Msa_20260706/scripts/run_delly_DEL_genotype_one.sh) | 运行：delly DEL genotype one（据文件名） |
| [`09_short_read_SV_Msa_20260706/scripts/run_delly_DEL_one.sh`](09_short_read_SV_Msa_20260706/scripts/run_delly_DEL_one.sh) | 运行：delly DEL one（据文件名） |
| [`09_short_read_SV_Msa_20260706/scripts/run_manta_one.sh`](09_short_read_SV_Msa_20260706/scripts/run_manta_one.sh) | 运行：manta one（据文件名） |
| [`09_short_read_SV_Msa_20260706/scripts/run_smoove_one.sh`](09_short_read_SV_Msa_20260706/scripts/run_smoove_one.sh) | 运行：smoove one（据文件名） |
| [`09_short_read_SV_Msa_20260706/scripts/summarize_chr23997_delly_DEL.py`](09_short_read_SV_Msa_20260706/scripts/summarize_chr23997_delly_DEL.py) | 汇总统计：chr23997 delly DEL（据文件名） |
| [`11_Chr23997_intron2_TE_check_20260707/scripts/check_Chr23997_intron2_TE.py`](11_Chr23997_intron2_TE_check_20260707/scripts/check_Chr23997_intron2_TE.py) | 检查：Chr23997 intron2 TE（据文件名） |
| [`12_Chr23997_Msa_R108_intron_insertion_TE_20260707/scripts/annotate_Msa_R108_homologous_insertions_TE.py`](12_Chr23997_Msa_R108_intron_insertion_TE_20260707/scripts/annotate_Msa_R108_homologous_insertions_TE.py) | 注释：Msa R108 homologous insertions TE（据文件名） |
| [`12_Chr23997_Msa_R108_intron_insertion_TE_20260707/scripts/extract_Msa_R108_intron_insertions_TE.py`](12_Chr23997_Msa_R108_intron_insertion_TE_20260707/scripts/extract_Msa_R108_intron_insertions_TE.py) | 提取：Msa R108 intron insertions TE（据文件名） |
| [`12_Chr23997_Msa_R108_intron_insertion_TE_20260707/scripts/locate_R108_homologous_insertion.py`](12_Chr23997_Msa_R108_intron_insertion_TE_20260707/scripts/locate_R108_homologous_insertion.py) | 定位：R108 homologous insertion（据文件名） |
| [`13_Chr23997_spiny_intron_window_comparison_20260708/scripts/annotate_spiny_Chr23997_windows_TE.py`](13_Chr23997_spiny_intron_window_comparison_20260708/scripts/annotate_spiny_Chr23997_windows_TE.py) | 注释：spiny Chr23997 windows TE（据文件名） |
| [`13_Chr23997_spiny_intron_window_comparison_20260708/scripts/compare_spiny_Chr23997_windows.py`](13_Chr23997_spiny_intron_window_comparison_20260708/scripts/compare_spiny_Chr23997_windows.py) | 比较：spiny Chr23997 windows（据文件名） |
| [`15_Chr23997_Msa_reference_similarity_RBH_20260708/scripts/plot_chr23997_msa_similarity_tracks.R`](15_Chr23997_Msa_reference_similarity_RBH_20260708/scripts/plot_chr23997_msa_similarity_tracks.R) | 绘图：chr23997 msa similarity tracks（据文件名） |
| [`15_Chr23997_Msa_reference_similarity_RBH_20260708/scripts/prepare_chr23997_similarity_tracks.py`](15_Chr23997_Msa_reference_similarity_RBH_20260708/scripts/prepare_chr23997_similarity_tracks.py) | 准备输入：chr23997 similarity tracks（据文件名） |
| [`16_Chr23997_Msa_reference_similarity_user_grouped_20260708/scripts/plot_chr23997_msa_similarity_tracks.R`](16_Chr23997_Msa_reference_similarity_user_grouped_20260708/scripts/plot_chr23997_msa_similarity_tracks.R) | 绘图：chr23997 msa similarity tracks（据文件名） |
| [`16_Chr23997_Msa_reference_similarity_user_grouped_20260708/scripts/prepare_chr23997_similarity_tracks.py`](16_Chr23997_Msa_reference_similarity_user_grouped_20260708/scripts/prepare_chr23997_similarity_tracks.py) | 准备输入：chr23997 similarity tracks（据文件名） |
| [`18_Chr23997_genomewide_cosegregation_rank_20260708/scripts/rank_sv_trait_cosegregation.py`](18_Chr23997_genomewide_cosegregation_rank_20260708/scripts/rank_sv_trait_cosegregation.py) | 排序：sv trait cosegregation（据文件名） |
| [`19_Chr23997_callset_truth_calibration_20260708/scripts/audit_exact_target_caller_records.py`](19_Chr23997_callset_truth_calibration_20260708/scripts/audit_exact_target_caller_records.py) | 核查：exact target caller records（据文件名） |
| [`19_Chr23997_callset_truth_calibration_20260708/scripts/build_chr23997_truth_calibration.py`](19_Chr23997_callset_truth_calibration_20260708/scripts/build_chr23997_truth_calibration.py) | 构建：chr23997 truth calibration（据文件名） |
| [`20_Chr23997_candidate_figure_20260708/scripts/draw_Chr23997_candidate_figure.R`](20_Chr23997_candidate_figure_20260708/scripts/draw_Chr23997_candidate_figure.R) | Base-R publication-style figure for the Chr23997 target intron deletion. |
| [`20_Chr23997_SNP_INDEL_SV_peak_20260708/scripts/plot_chr23997_region_gwas.py`](20_Chr23997_SNP_INDEL_SV_peak_20260708/scripts/plot_chr23997_region_gwas.py) | 绘图：chr23997 region gwas（据文件名） |
| [`20_Chr23997_SNP_INDEL_SV_peak_20260708/scripts/run_chr23997_peak_assoc.sh`](20_Chr23997_SNP_INDEL_SV_peak_20260708/scripts/run_chr23997_peak_assoc.sh) | 运行：chr23997 peak assoc（据文件名） |
| [`20_Chr23997_SNP_INDEL_SV_peak_20260708/scripts/summarize_gwas_region.py`](20_Chr23997_SNP_INDEL_SV_peak_20260708/scripts/summarize_gwas_region.py) | 汇总统计：gwas region（据文件名） |
| [`21_genomewide_SV_trait_cosegregation_scan_20260708/scripts/run_genomewide_sv_trait_scan.R`](21_genomewide_SV_trait_cosegregation_scan_20260708/scripts/run_genomewide_sv_trait_scan.R) | 运行：genomewide sv trait scan（据文件名） |
| [`22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/plot_raw_readmapping_sv_cosegregation.R`](22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/plot_raw_readmapping_sv_cosegregation.R) | 绘图：raw readmapping sv cosegregation（据文件名） |
| [`22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/plot_raw_readmapping_sv_cosegregation_genomewide.R`](22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/plot_raw_readmapping_sv_cosegregation_genomewide.R) | 绘图：raw readmapping sv cosegregation genomewide（据文件名） |
| [`22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/plot_raw_readmapping_sv_cosegregation_genomewide_jittered.R`](22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/plot_raw_readmapping_sv_cosegregation_genomewide_jittered.R) | 绘图：raw readmapping sv cosegregation genomewide jittered（据文件名） |
| [`22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/scan_raw_readmapping_sv_cosegregation.py`](22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/scan_raw_readmapping_sv_cosegregation.py) | 扫描：raw readmapping sv cosegregation（据文件名） |
| [`22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/scan_raw_readmapping_sv_cosegregation_genomewide.py`](22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/scan_raw_readmapping_sv_cosegregation_genomewide.py) | 扫描：raw readmapping sv cosegregation genomewide（据文件名） |
| [`22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/scan_raw_readmapping_sv_cosegregation_genomewide_fast.py`](22_readmapping_rawSV_trait_cosegregation_scan_20260708/scripts/scan_raw_readmapping_sv_cosegregation_genomewide_fast.py) | 扫描：raw readmapping sv cosegregation genomewide fast（据文件名） |
| [`24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/make_codon12_alignment.py`](24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/make_codon12_alignment.py) | 生成：codon12 alignment（据文件名） |
| [`24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/pad_fasta_alignment.py`](24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/pad_fasta_alignment.py) | 补齐：fasta alignment（据文件名） |
| [`24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/plot_newick_baseR.R`](24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/plot_newick_baseR.R) | 绘图：newick baseR（据文件名） |
| [`24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/prepare_chr23997_tree_inputs_with_outgroup.py`](24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/prepare_chr23997_tree_inputs_with_outgroup.py) | 准备输入：chr23997 tree inputs with outgroup（据文件名） |
| [`24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/run_iqtree_rooted.sh`](24_Chr23997_strict_RBH_gene_tree_with_outgroup_20260709/scripts/run_iqtree_rooted.sh) | 运行：iqtree rooted（据文件名） |
| [`25_SVGAP_supported_PAV_trait_scan_20260709/scripts/plot_svgap_supported_PAV_trait_cosegregation.R`](25_SVGAP_supported_PAV_trait_scan_20260709/scripts/plot_svgap_supported_PAV_trait_cosegregation.R) | 绘图：svgap supported PAV trait cosegregation（据文件名） |
| [`25_SVGAP_supported_PAV_trait_scan_20260709/scripts/plot_svgap_supported_PAV_trait_cosegregation_binned_bars.R`](25_SVGAP_supported_PAV_trait_scan_20260709/scripts/plot_svgap_supported_PAV_trait_cosegregation_binned_bars.R) | 绘图：svgap supported PAV trait cosegregation binned bars（据文件名） |
| [`25_SVGAP_supported_PAV_trait_scan_20260709/scripts/plot_svgap_supported_PAV_trait_cosegregation_hybrid_raster_scatter.R`](25_SVGAP_supported_PAV_trait_scan_20260709/scripts/plot_svgap_supported_PAV_trait_cosegregation_hybrid_raster_scatter.R) | 绘图：svgap supported PAV trait cosegregation hybrid raster scatter（据文件名） |
| [`25_SVGAP_supported_PAV_trait_scan_20260709/scripts/scan_svgap_supported_pav_trait_cosegregation.py`](25_SVGAP_supported_PAV_trait_scan_20260709/scripts/scan_svgap_supported_pav_trait_cosegregation.py) | 扫描：svgap supported pav trait cosegregation（据文件名） |
| [`26_Msa_SPL_family_phylogeny_20260710/scripts/run_msa_spl_family_phylogeny.sh`](26_Msa_SPL_family_phylogeny_20260710/scripts/run_msa_spl_family_phylogeny.sh) | 运行：msa spl family phylogeny（据文件名） |
| [`26_Msa_SPL_family_phylogeny_20260710/scripts/run_msa_spl_family_phylogeny_mt_split.sh`](26_Msa_SPL_family_phylogeny_20260710/scripts/run_msa_spl_family_phylogeny_mt_split.sh) | 运行：msa spl family phylogeny mt split（据文件名） |
| [`27_R108_SPL_family_phylogeny_20260710/scripts/run_r108_spl_family_phylogeny_mt.sh`](27_R108_SPL_family_phylogeny_20260710/scripts/run_r108_spl_family_phylogeny_mt.sh) | 运行：r108 spl family phylogeny mt（据文件名） |
| [`28_GWAS_spiny/01_input_audit/scripts/audit_sv_qc_and_chr23997.py`](28_GWAS_spiny/01_input_audit/scripts/audit_sv_qc_and_chr23997.py) | 核查：sv qc and chr23997（据文件名） |
| [`28_GWAS_spiny/01_input_audit/scripts/inspect_inputs.py`](28_GWAS_spiny/01_input_audit/scripts/inspect_inputs.py) | 检查：inputs（据文件名） |
| [`28_GWAS_spiny/01_input_audit/scripts/inspect_tools.sh`](28_GWAS_spiny/01_input_audit/scripts/inspect_tools.sh) | 检查：tools（据文件名） |
| [`28_GWAS_spiny/01_input_audit/scripts/prepare_samples.py`](28_GWAS_spiny/01_input_audit/scripts/prepare_samples.py) | 准备输入：samples（据文件名） |
| [`28_GWAS_spiny/02_snp_indel_qc/scripts/run_prepare_snp_indel.sh`](28_GWAS_spiny/02_snp_indel_qc/scripts/run_prepare_snp_indel.sh) | 运行：prepare snp indel（据文件名） |
| [`28_GWAS_spiny/02_snp_indel_qc/scripts/set_fam_phenotype.py`](28_GWAS_spiny/02_snp_indel_qc/scripts/set_fam_phenotype.py) | 设置：fam phenotype（据文件名） |
| [`28_GWAS_spiny/03_population_structure/scripts/build_structure_inputs.py`](28_GWAS_spiny/03_population_structure/scripts/build_structure_inputs.py) | 构建：structure inputs（据文件名） |
| [`28_GWAS_spiny/03_population_structure/scripts/run_03_then_04.sh`](28_GWAS_spiny/03_population_structure/scripts/run_03_then_04.sh) | 运行：03 then 04（据文件名） |
| [`28_GWAS_spiny/03_population_structure/scripts/run_population_structure.sh`](28_GWAS_spiny/03_population_structure/scripts/run_population_structure.sh) | 运行：population structure（据文件名） |
| [`28_GWAS_spiny/04_snp_indel_gwas/scripts/run_snp_indel_gwas.sh`](28_GWAS_spiny/04_snp_indel_gwas/scripts/run_snp_indel_gwas.sh) | 运行：snp indel gwas（据文件名） |
| [`28_GWAS_spiny/04_snp_indel_gwas/scripts/summarize_gemma.py`](28_GWAS_spiny/04_snp_indel_gwas/scripts/summarize_gemma.py) | 汇总统计：gemma（据文件名） |
| [`28_GWAS_spiny/05_sv_genotyping/scripts/filter_smoove_sites.py`](28_GWAS_spiny/05_sv_genotyping/scripts/filter_smoove_sites.py) | Filter Smoove/SVtools cohort sites without performing any SV clustering. |
| [`28_GWAS_spiny/05_sv_genotyping/scripts/remove_orphan_bnd.py`](28_GWAS_spiny/05_sv_genotyping/scripts/remove_orphan_bnd.py) | Remove BND records whose MATEID partner was filtered out, preserving order. |
| [`28_GWAS_spiny/05_sv_genotyping/scripts/run_chr23997_graphtyper.sh`](28_GWAS_spiny/05_sv_genotyping/scripts/run_chr23997_graphtyper.sh) | 运行：chr23997 graphtyper（据文件名） |
| [`28_GWAS_spiny/05_sv_genotyping/scripts/run_sv_cohort_qc.sh`](28_GWAS_spiny/05_sv_genotyping/scripts/run_sv_cohort_qc.sh) | 运行：sv cohort qc（据文件名） |
| [`28_GWAS_spiny/05_sv_genotyping/scripts/run_sv_genotype_task.sh`](28_GWAS_spiny/05_sv_genotyping/scripts/run_sv_genotype_task.sh) | 运行：sv genotype task（据文件名） |
| [`28_GWAS_spiny/05_sv_genotyping/scripts/run_sv_site_discovery.sh`](28_GWAS_spiny/05_sv_genotyping/scripts/run_sv_site_discovery.sh) | 运行：sv site discovery（据文件名） |
| [`28_GWAS_spiny/05_sv_genotyping/scripts/validate_bnd_pairs.py`](28_GWAS_spiny/05_sv_genotyping/scripts/validate_bnd_pairs.py) | Validate MATEID completeness and adjacency in a Smoove sites VCF. |
| [`28_GWAS_spiny/05_sv_genotyping/scripts/validate_existing_sv_genotypes.sh`](28_GWAS_spiny/05_sv_genotyping/scripts/validate_existing_sv_genotypes.sh) | 验证：existing sv genotypes（据文件名） |
| [`28_GWAS_spiny/06_sv_gwas/scripts/run_sv_gwas.sh`](28_GWAS_spiny/06_sv_gwas/scripts/run_sv_gwas.sh) | 运行：sv gwas（据文件名） |
| [`28_GWAS_spiny/06_sv_gwas/scripts/summarize_sv_gemma.py`](28_GWAS_spiny/06_sv_gwas/scripts/summarize_sv_gemma.py) | 汇总统计：sv gemma（据文件名） |
| [`28_GWAS_spiny/06_sv_gwas/scripts/vcf_to_bimbam.py`](28_GWAS_spiny/06_sv_gwas/scripts/vcf_to_bimbam.py) |  |
| [`28_GWAS_spiny/07_joint_gwas/scripts/finalize_stage07.sh`](28_GWAS_spiny/07_joint_gwas/scripts/finalize_stage07.sh) | 收尾汇总：stage07（据文件名） |
| [`28_GWAS_spiny/07_joint_gwas/scripts/integrate_gwas_results.py`](28_GWAS_spiny/07_joint_gwas/scripts/integrate_gwas_results.py) | 整合：gwas results（据文件名） |
| [`28_GWAS_spiny/07_joint_gwas/scripts/plot_joint_gwas.R`](28_GWAS_spiny/07_joint_gwas/scripts/plot_joint_gwas.R) | 绘图：joint gwas（据文件名） |
| [`28_GWAS_spiny/07_joint_gwas/scripts/run_gwas_sensitivity.sh`](28_GWAS_spiny/07_joint_gwas/scripts/run_gwas_sensitivity.sh) | 运行：gwas sensitivity（据文件名） |
| [`28_GWAS_spiny/07_joint_gwas/scripts/run_stage07_integration.sh`](28_GWAS_spiny/07_joint_gwas/scripts/run_stage07_integration.sh) | 运行：stage07 integration（据文件名） |
| [`28_GWAS_spiny/07_joint_gwas/scripts/summarize_sensitivity.py`](28_GWAS_spiny/07_joint_gwas/scripts/summarize_sensitivity.py) | 汇总统计：sensitivity（据文件名） |
| [`28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/extract_chr23997_bam_evidence.py`](28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/extract_chr23997_bam_evidence.py) | Extract independent short-read evidence for the validated Chr23997 deletion. |
| [`28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/finalize_target_analysis.sh`](28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/finalize_target_analysis.sh) | 收尾汇总：target analysis（据文件名） |
| [`28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/prepare_target_manifest.py`](28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/prepare_target_manifest.py) | Build the exact 143-sample manifest in the established GWAS order. |
| [`28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/run_target_genotype_task.sh`](28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/run_target_genotype_task.sh) | 运行：target genotype task（据文件名） |
| [`28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/run_target_pilot.sh`](28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/run_target_pilot.sh) | 运行：target pilot（据文件名） |
| [`28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/summarize_target_genotypes.py`](28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/scripts/summarize_target_genotypes.py) | Combine exact-site callers and independent BAM evidence without imputing failures. |
| [`28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/build_candidate_inputs.py`](28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/build_candidate_inputs.py) | Build typed Chr23997 candidate-SV inputs without imputing ambiguous calls. |
| [`28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/conditional_analysis.py`](28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/conditional_analysis.py) | Run local tag-SNP conditional logistic models on the common complete-case set. |
| [`28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/local_gwas_and_ld.py`](28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/local_gwas_and_ld.py) | Collect Chr23997 +/- window association results and un-imputed target-marker LD. |
| [`28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/plot_local_gwas.R`](28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/plot_local_gwas.R) | 绘图：local gwas（据文件名） |
| [`28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/run_candidate_association.py`](28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/run_candidate_association.py) | Run direct candidate-SV association analyses without imputing missing genotypes. |
| [`28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/run_candidate_gmmat.R`](28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/run_candidate_gmmat.R) | 运行：candidate gmmat（据文件名） |
| [`28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/write_validation_report.py`](28_GWAS_spiny/09_Chr23997_candidate_GWAS_validation_20260715/scripts/write_validation_report.py) | Write an auditable candidate-only GWAS validation report for Chr23997. |
| [`29_RcMYB106_promoter_SV_R108_M22_20260714/scripts/analyze_rcmyb106_promoter.py`](29_RcMYB106_promoter_SV_R108_M22_20260714/scripts/analyze_rcmyb106_promoter.py) | 分析：rcmyb106 promoter（据文件名） |
| [`29_RcMYB106_promoter_SV_R108_M22_20260714/scripts/compare_family_promoters.py`](29_RcMYB106_promoter_SV_R108_M22_20260714/scripts/compare_family_promoters.py) | 比较：family promoters（据文件名） |
| [`29_RcMYB106_promoter_SV_R108_M22_20260714/scripts/run_pipeline.sh`](29_RcMYB106_promoter_SV_R108_M22_20260714/scripts/run_pipeline.sh) | 运行：pipeline（据文件名） |
| [`29_RcMYB106_promoter_SV_R108_M22_20260714/scripts/screen_reciprocal_candidates.py`](29_RcMYB106_promoter_SV_R108_M22_20260714/scripts/screen_reciprocal_candidates.py) | 筛选：reciprocal candidates（据文件名） |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/chr23997_breakpoint_genotyper.py`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/chr23997_breakpoint_genotyper.py) | Genotype the Chr23997 intron-2 PAV directly from BAM split-read evidence. |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/chr23997_candidate_diagnostics.py`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/chr23997_candidate_diagnostics.py) | Population-aware diagnostics for the Chr23997 intron-2 core-loss state. |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/chr23997_haplotype_genotyper.py`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/chr23997_haplotype_genotyper.py) | Genotype the Chr23997 intron-2 deletion by local two-haplotype alignment. |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/prepare_chr23997_haplotypes.py`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/prepare_chr23997_haplotypes.py) | Create reference and exact-DEL local haplotypes for Chr23997. |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/run_chr23997_breakpoint_all.sh`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/run_chr23997_breakpoint_all.sh) | 运行：chr23997 breakpoint all（据文件名） |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/run_chr23997_breakpoint_pilot.sh`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/run_chr23997_breakpoint_pilot.sh) | 运行：chr23997 breakpoint pilot（据文件名） |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/run_chr23997_haplotype_task.sh`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/run_chr23997_haplotype_task.sh) | 运行：chr23997 haplotype task（据文件名） |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/run_chr23997_pilot.sh`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/run_chr23997_pilot.sh) | 运行：chr23997 pilot（据文件名） |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/summarize_chr23997_breakpoint_calls.py`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/summarize_chr23997_breakpoint_calls.py) | Merge Chr23997 core-loss genotypes and report the unadjusted contingency test. |
| [`30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/summarize_chr23997_haplotype_calls.py`](30_Chr23997_targeted_graphgenotype_GWAS_20260716/scripts/summarize_chr23997_haplotype_calls.py) | Merge per-sample Chr23997 local haplotype calls with phenotype metadata. |
| [`31_Chr23997_variable_core_haplotype_20260717/scripts/build_variable_core_haplotype.py`](31_Chr23997_variable_core_haplotype_20260717/scripts/build_variable_core_haplotype.py) | Summarize a variable-length Chr23997 intron-2 core-loss haplotype. |
| [`31_Chr23997_variable_core_haplotype_20260717/scripts/summarize_wgs_variable_core.py`](31_Chr23997_variable_core_haplotype_20260717/scripts/summarize_wgs_variable_core.py) | Summarize direct read-based variable core-loss calls across depth cutoffs. |
| [`32_PAV_trait_cosegregation_missing_aware_20260802/scripts/plot_ED7b_missing_aware_scan.R`](32_PAV_trait_cosegregation_missing_aware_20260802/scripts/plot_ED7b_missing_aware_scan.R) | Extended Data Fig. 7b (revised): genome-wide phenotype-segregation scan of high-confiden… |
| [`32_PAV_trait_cosegregation_missing_aware_20260802/scripts/reanalyze_missing_aware_pav.py`](32_PAV_trait_cosegregation_missing_aware_20260802/scripts/reanalyze_missing_aware_pav.py) | Missing-aware PAV/phenotype association reanalysis. |
| [`scripts/add_event_coordinates_to_priority.py`](scripts/add_event_coordinates_to_priority.py) |  |
| [`scripts/annotate_priority_with_genes.py`](scripts/annotate_priority_with_genes.py) | 注释：priority with genes（据文件名） |
| [`scripts/annotate_priority_with_te.py`](scripts/annotate_priority_with_te.py) | 注释：priority with te（据文件名） |
| [`scripts/make_igv_reports_for_29SV.py`](scripts/make_igv_reports_for_29SV.py) | Build igv-reports HTML pages for pod-spiny SV validation candidates. |
| [`scripts/make_igv_validation_package.py`](scripts/make_igv_validation_package.py) | 生成：igv validation package（据文件名） |
| [`scripts/make_manual_sv_review_table.py`](scripts/make_manual_sv_review_table.py) | 生成：manual sv review table（据文件名） |
| [`scripts/prioritize_spiny_candidates.py`](scripts/prioritize_spiny_candidates.py) | 优先级排序：spiny candidates（据文件名） |
| [`scripts/screen_direction_free_all_candidates.py`](scripts/screen_direction_free_all_candidates.py) | Direction-free SV/phenotype cosegregation screen for pod-spiny candidates. |
| [`scripts/screen_spiny_cosegregation.py`](scripts/screen_spiny_cosegregation.py) | 筛选：spiny cosegregation（据文件名） |
| [`scripts/summarize_igv_read_support.py`](scripts/summarize_igv_read_support.py) | 汇总统计：igv read support（据文件名） |
| [`scripts/tier_spiny_candidates.py`](scripts/tier_spiny_candidates.py) | 分级：spiny candidates（据文件名） |
