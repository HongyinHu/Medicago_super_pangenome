# panSV

SV 鉴定与泛 SV 构建：双参考 HiFi reads 多软件鉴定（Snakemake + Jasmine + Sniffles2 回填基因型）、SVGAP 组装法、整合与 QC、SV 特征/TE/表达分析

| 脚本 | 说明 |
|---|---|
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/config/config.yaml`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/config/config.yaml) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/build_chunk_genotype_tasks.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/build_chunk_genotype_tasks.py) | 构建：chunk genotype tasks（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/check_inputs.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/check_inputs.py) | 检查：inputs（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/concat_chunk_vcfs.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/concat_chunk_vcfs.py) | 拼接：chunk vcfs（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match.py) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match_publication.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match_publication.py) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match_strict.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match_strict.py) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/make_meta_pav_matrix.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/make_meta_pav_matrix.py) | 生成：meta pav matrix（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/make_pav_matrix.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/make_pav_matrix.py) | 生成：pav matrix（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/merge_genotyped_vcfs.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/merge_genotyped_vcfs.py) | 合并：genotyped vcfs（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/plot_msa_read_mapping_svtype_by_species.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/plot_msa_read_mapping_svtype_by_species.py) | 绘图：msa read mapping svtype by species（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/plot_msa_read_mapping_te_origin.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/plot_msa_read_mapping_te_origin.py) | 绘图：msa read mapping te origin（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_dualref_panSV.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_dualref_panSV.sh) | 运行：dualref panSV（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_dualref_panSV_chunked_genotype.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_dualref_panSV_chunked_genotype.sh) | 运行：dualref panSV chunked genotype（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_genotype_chunk_array_task.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_genotype_chunk_array_task.sh) | 运行：genotype chunk array task（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_publication_meta_panSV.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_publication_meta_panSV.sh) | 运行：publication meta panSV（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_strict_meta_panSV.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_strict_meta_panSV.sh) | 运行：strict meta panSV（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/sort_vcf.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/sort_vcf.py) | 排序：vcf（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/split_discovery_vcf.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/split_discovery_vcf.py) | 拆分：discovery vcf（据文件名） |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/vcf_filter.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/vcf_filter.py) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/Snakefile`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/Snakefile) |  |
| [`01.read_based_dualref_hifi/scripts/run_readbased_task.sh`](01.read_based_dualref_hifi/scripts/run_readbased_task.sh) | 运行：readbased task（据文件名） |
| [`01.read_based_dualref_hifi/scripts/run_worker.sh`](01.read_based_dualref_hifi/scripts/run_worker.sh) | 运行：worker（据文件名） |
| [`02.assembly_svgap_dualref/scripts/plot_msa_ref_svgap_svtype_counts.py`](02.assembly_svgap_dualref/scripts/plot_msa_ref_svgap_svtype_counts.py) | 绘图：msa ref svgap svtype counts（据文件名） |
| [`02.assembly_svgap_dualref/scripts/prepare_svgap_inputs.sh`](02.assembly_svgap_dualref/scripts/prepare_svgap_inputs.sh) | 准备输入：svgap inputs（据文件名） |
| [`02.assembly_svgap_dualref/scripts/run_svgap_ref_stages.sh`](02.assembly_svgap_dualref/scripts/run_svgap_ref_stages.sh) | 运行：svgap ref stages（据文件名） |
| [`02.assembly_svgap_dualref/scripts/run_svgap_wga_task.sh`](02.assembly_svgap_dualref/scripts/run_svgap_wga_task.sh) | 运行：svgap wga task（据文件名） |
| [`03.integrated_panSV/scripts/integrate_panSV.py`](03.integrated_panSV/scripts/integrate_panSV.py) | 整合：panSV（据文件名） |
| [`03.integrated_panSV/scripts/run_integrate_panSV.sh`](03.integrated_panSV/scripts/run_integrate_panSV.sh) | 运行：integrate panSV（据文件名） |
| [`04.integrated_panSV_Msa_single_ref/scripts/integrate_panSV_Msa_single_ref.py`](04.integrated_panSV_Msa_single_ref/scripts/integrate_panSV_Msa_single_ref.py) | 整合：panSV Msa single ref（据文件名） |
| [`04.integrated_panSV_Msa_single_ref/scripts/plot_msa_svtype_by_species.py`](04.integrated_panSV_Msa_single_ref/scripts/plot_msa_svtype_by_species.py) | 绘图：msa svtype by species（据文件名） |
| [`04.integrated_panSV_Msa_single_ref/scripts/run_integrate_panSV_Msa_single_ref.sh`](04.integrated_panSV_Msa_single_ref/scripts/run_integrate_panSV_Msa_single_ref.sh) | 运行：integrate panSV Msa single ref（据文件名） |
| [`05.integrated_panSV_R108_single_ref/scripts/integrate_panSV_R108_single_ref.py`](05.integrated_panSV_R108_single_ref/scripts/integrate_panSV_R108_single_ref.py) | 整合：panSV R108 single ref（据文件名） |
| [`05.integrated_panSV_R108_single_ref/scripts/run_integrate_panSV_R108_single_ref.sh`](05.integrated_panSV_R108_single_ref/scripts/run_integrate_panSV_R108_single_ref.sh) | 运行：integrate panSV R108 single ref（据文件名） |
| [`06.integrated_panSV_Msa_paper_style_unified/scripts/normalize_integrated_pansv_msa_genome474a_noA17.py`](06.integrated_panSV_Msa_paper_style_unified/scripts/normalize_integrated_pansv_msa_genome474a_noA17.py) | Normalize Msa single-reference integrated panSV tables. |
| [`06.integrated_panSV_Msa_paper_style_unified/scripts/validate_paper_style_unified_msa.py`](06.integrated_panSV_Msa_paper_style_unified/scripts/validate_paper_style_unified_msa.py) | Validate the genome_Msa paper-style unified panSV outputs. |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/build_per_species_merged_sv_msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/build_per_species_merged_sv_msa.py) | 构建：per species merged sv msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/build_species_level_pansv_msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/build_species_level_pansv_msa.py) | 构建：species level pansv msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/flank_coverage_filter_vcf.py`](07.read_mapping_flankQC_panSV_Msa/scripts/flank_coverage_filter_vcf.py) |  |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/integrate_svgap_read_primary_msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/integrate_svgap_read_primary_msa.py) | 整合：svgap read primary msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/make_exclude_genome_Msa_plot_inputs.py`](07.read_mapping_flankQC_panSV_Msa/scripts/make_exclude_genome_Msa_plot_inputs.py) | 生成：exclude genome Msa plot inputs（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/make_pav_from_jasmine_vcf.py`](07.read_mapping_flankQC_panSV_Msa/scripts/make_pav_from_jasmine_vcf.py) | 生成：pav from jasmine vcf（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/make_supplementary_table13_sv_distribution.py`](07.read_mapping_flankQC_panSV_Msa/scripts/make_supplementary_table13_sv_distribution.py) | 生成：supplementary table13 sv distribution（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie.py`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie.py) | 绘图：Msa svtype sorted bar pansv pie（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie.R) | 绘图：Msa svtype sorted bar pansv pie（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie_split_v4.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie_split_v4.R) | 绘图：Msa svtype sorted bar pansv pie split v4（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exact_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exact_exclude_Msa.R) | 绘图：pan core SV accumulation exact exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exact_exclude_Msa_split_v2.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exact_exclude_Msa_split_v2.R) | 绘图：pan core SV accumulation exact exclude Msa split v2（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exclude_Msa.R) | 绘图：pan core SV accumulation exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_panSV_DEL_INS_length_distribution_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_panSV_DEL_INS_length_distribution_exclude_Msa.R) | 绘图：panSV DEL INS length distribution exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size.R) | 绘图：PAV detection frequency by size（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size_DEL_INS.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size_DEL_INS.R) | Extended Data Fig. 6b (revised): PAV detection frequency by size, DEL/INS only. |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size_DEL_INS_figure.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size_DEL_INS_figure.R) | Extended Data Fig. 6b (revised) figure, drawn from the count table written by |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_gene_context_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_gene_context_exclude_Msa.R) | 绘图：PAV gene context exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_gene_structure_per_species_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_gene_structure_per_species_exclude_Msa.R) | 绘图：PAV gene structure per species exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_shared_specific_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_shared_specific_exclude_Msa.R) | 绘图：PAV shared specific exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_association_4class_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_association_4class_exclude_Msa.R) | 绘图：PAV TE association 4class exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_association_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_association_exclude_Msa.R) | 绘图：PAV TE association exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_gene_distance_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_gene_distance_exclude_Msa.R) | 绘图：PAV TE gene distance exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/run_Msa_flankQC_panSV.sh`](07.read_mapping_flankQC_panSV_Msa/scripts/run_Msa_flankQC_panSV.sh) | 运行：Msa flankQC panSV（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/sort_vcf.py`](07.read_mapping_flankQC_panSV_Msa/scripts/sort_vcf.py) | 排序：vcf（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_pan_core_SV_accumulation_exact_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_pan_core_SV_accumulation_exact_exclude_Msa.py) | 统计：pan core SV accumulation exact exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_pan_core_SV_accumulation_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_pan_core_SV_accumulation_exclude_Msa.py) | 统计：pan core SV accumulation exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_panSV_DEL_INS_length_distribution_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_panSV_DEL_INS_length_distribution_exclude_Msa.py) | 统计：panSV DEL INS length distribution exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_gene_context_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_gene_context_exclude_Msa.py) | 统计：PAV gene context exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_gene_structure_per_species_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_gene_structure_per_species_exclude_Msa.py) | 统计：PAV gene structure per species exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_shared_specific_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_shared_specific_exclude_Msa.py) | 统计：PAV shared specific exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_TE_association_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_TE_association_exclude_Msa.py) | 统计：PAV TE association exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_TE_gene_distance_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_TE_gene_distance_exclude_Msa.py) | 统计：PAV TE gene distance exclude Msa（据文件名） |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/vcf_filter_msa_with_tra.py`](07.read_mapping_flankQC_panSV_Msa/scripts/vcf_filter_msa_with_tra.py) |  |
| [`08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_mixed_illumina_expression_sv_distance.R`](08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_mixed_illumina_expression_sv_distance.R) | 绘图：mixed illumina expression sv distance（据文件名） |
| [`08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_pav_gene_structure_expression_divergence.R`](08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_pav_gene_structure_expression_divergence.R) | 绘图：pav gene structure expression divergence（据文件名） |
| [`08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_pav_gene_structure_fc_percent.R`](08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_pav_gene_structure_fc_percent.R) | 绘图：pav gene structure fc percent（据文件名） |
| [`08_mixed_illumina_expression_SV_distance_20260705/scripts/run_mixed_illumina_expression_sv_distance.sh`](08_mixed_illumina_expression_SV_distance_20260705/scripts/run_mixed_illumina_expression_sv_distance.sh) | 运行：mixed illumina expression sv distance（据文件名） |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_expression_level_by_sv.R`](09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_expression_level_by_sv.R) | 绘图：strict expression level by sv（据文件名） |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_ortholog_expression.R`](09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_ortholog_expression.R) | 绘图：strict ortholog expression（据文件名） |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_pearson_distance_by_sv.R`](09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_pearson_distance_by_sv.R) | 绘图：strict pearson distance by sv（据文件名） |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/run_strict_ortholog_expression_sv.sh`](09_strict_ortholog_expression_SV_distance_20260705/scripts/run_strict_ortholog_expression_sv.sh) | 运行：strict ortholog expression sv（据文件名） |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/strict_ortholog_expression_sv.py`](09_strict_ortholog_expression_SV_distance_20260705/scripts/strict_ortholog_expression_sv.py) |  |
| [`10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705/scripts/plot_read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson.R`](10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705/scripts/plot_read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson.R) | 绘图：read mapping primary exon CDS UTR promoter2kb pearson（据文件名） |
| [`10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705/scripts/plot_read_mapping_primary_gene_structure_pearson.R`](10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705/scripts/plot_read_mapping_primary_gene_structure_pearson.R) | 绘图：read mapping primary gene structure pearson（据文件名） |
