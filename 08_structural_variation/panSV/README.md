# panSV

SV calling and pan-SV construction: dual-reference HiFi read-based calling (Snakemake + Jasmine + Sniffles2 genotyping), SVGAP assembly-based calling, integration and QC, SV feature/TE/expression analyses

| Script | Description |
|---|---|
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/config/config.yaml`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/config/config.yaml) | Configuration (YAML) |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/build_chunk_genotype_tasks.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/build_chunk_genotype_tasks.py) | Build chunk genotype tasks |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/check_inputs.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/check_inputs.py) | Check inputs |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/concat_chunk_vcfs.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/concat_chunk_vcfs.py) | Concat chunk vcfs |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match.py) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match_publication.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match_publication.py) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match_strict.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/cross_ref_match_strict.py) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/make_meta_pav_matrix.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/make_meta_pav_matrix.py) | Make meta pav matrix |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/make_pav_matrix.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/make_pav_matrix.py) | Make pav matrix |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/merge_genotyped_vcfs.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/merge_genotyped_vcfs.py) | Merge genotyped vcfs |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/plot_msa_read_mapping_svtype_by_species.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/plot_msa_read_mapping_svtype_by_species.py) | Plot msa read mapping svtype by species |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/plot_msa_read_mapping_te_origin.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/plot_msa_read_mapping_te_origin.py) | Plot msa read mapping te origin |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_dualref_panSV.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_dualref_panSV.sh) | Run dualref panSV |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_dualref_panSV_chunked_genotype.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_dualref_panSV_chunked_genotype.sh) | Run dualref panSV chunked genotype |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_genotype_chunk_array_task.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_genotype_chunk_array_task.sh) | Run genotype chunk array task |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_publication_meta_panSV.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_publication_meta_panSV.sh) | Run publication meta panSV |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_strict_meta_panSV.sh`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/run_strict_meta_panSV.sh) | Run strict meta panSV |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/sort_vcf.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/sort_vcf.py) | Sort vcf |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/split_discovery_vcf.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/split_discovery_vcf.py) | Split discovery vcf |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/vcf_filter.py`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/scripts/vcf_filter.py) |  |
| [`01.read_based_dualref_hifi/05_dualref_panSV_20260626/Snakefile`](01.read_based_dualref_hifi/05_dualref_panSV_20260626/Snakefile) | Snakemake workflow |
| [`01.read_based_dualref_hifi/scripts/run_readbased_task.sh`](01.read_based_dualref_hifi/scripts/run_readbased_task.sh) | Run readbased task |
| [`01.read_based_dualref_hifi/scripts/run_worker.sh`](01.read_based_dualref_hifi/scripts/run_worker.sh) | Run worker |
| [`02.assembly_svgap_dualref/scripts/plot_msa_ref_svgap_svtype_counts.py`](02.assembly_svgap_dualref/scripts/plot_msa_ref_svgap_svtype_counts.py) | Plot msa ref svgap svtype counts |
| [`02.assembly_svgap_dualref/scripts/prepare_svgap_inputs.sh`](02.assembly_svgap_dualref/scripts/prepare_svgap_inputs.sh) | Prepare svgap inputs |
| [`02.assembly_svgap_dualref/scripts/run_svgap_ref_stages.sh`](02.assembly_svgap_dualref/scripts/run_svgap_ref_stages.sh) | Run svgap ref stages |
| [`02.assembly_svgap_dualref/scripts/run_svgap_wga_task.sh`](02.assembly_svgap_dualref/scripts/run_svgap_wga_task.sh) | Run svgap wga task |
| [`03.integrated_panSV/scripts/integrate_panSV.py`](03.integrated_panSV/scripts/integrate_panSV.py) | Integrate panSV |
| [`03.integrated_panSV/scripts/run_integrate_panSV.sh`](03.integrated_panSV/scripts/run_integrate_panSV.sh) | Run integrate panSV |
| [`04.integrated_panSV_Msa_single_ref/scripts/integrate_panSV_Msa_single_ref.py`](04.integrated_panSV_Msa_single_ref/scripts/integrate_panSV_Msa_single_ref.py) | Integrate panSV Msa single ref |
| [`04.integrated_panSV_Msa_single_ref/scripts/plot_msa_svtype_by_species.py`](04.integrated_panSV_Msa_single_ref/scripts/plot_msa_svtype_by_species.py) | Plot msa svtype by species |
| [`04.integrated_panSV_Msa_single_ref/scripts/run_integrate_panSV_Msa_single_ref.sh`](04.integrated_panSV_Msa_single_ref/scripts/run_integrate_panSV_Msa_single_ref.sh) | Run integrate panSV Msa single ref |
| [`05.integrated_panSV_R108_single_ref/scripts/integrate_panSV_R108_single_ref.py`](05.integrated_panSV_R108_single_ref/scripts/integrate_panSV_R108_single_ref.py) | Integrate panSV R108 single ref |
| [`05.integrated_panSV_R108_single_ref/scripts/run_integrate_panSV_R108_single_ref.sh`](05.integrated_panSV_R108_single_ref/scripts/run_integrate_panSV_R108_single_ref.sh) | Run integrate panSV R108 single ref |
| [`06.integrated_panSV_Msa_paper_style_unified/scripts/normalize_integrated_pansv_msa_genome474a_noA17.py`](06.integrated_panSV_Msa_paper_style_unified/scripts/normalize_integrated_pansv_msa_genome474a_noA17.py) | Normalize Msa single-reference integrated panSV tables. |
| [`06.integrated_panSV_Msa_paper_style_unified/scripts/validate_paper_style_unified_msa.py`](06.integrated_panSV_Msa_paper_style_unified/scripts/validate_paper_style_unified_msa.py) | Validate the genome_Msa paper-style unified panSV outputs. |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/build_per_species_merged_sv_msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/build_per_species_merged_sv_msa.py) | Build per species merged sv msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/build_species_level_pansv_msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/build_species_level_pansv_msa.py) | Build species level pansv msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/flank_coverage_filter_vcf.py`](07.read_mapping_flankQC_panSV_Msa/scripts/flank_coverage_filter_vcf.py) |  |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/integrate_svgap_read_primary_msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/integrate_svgap_read_primary_msa.py) | Integrate svgap read primary msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/make_exclude_genome_Msa_plot_inputs.py`](07.read_mapping_flankQC_panSV_Msa/scripts/make_exclude_genome_Msa_plot_inputs.py) | Make exclude genome Msa plot inputs |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/make_pav_from_jasmine_vcf.py`](07.read_mapping_flankQC_panSV_Msa/scripts/make_pav_from_jasmine_vcf.py) | Make pav from jasmine vcf |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/make_supplementary_table13_sv_distribution.py`](07.read_mapping_flankQC_panSV_Msa/scripts/make_supplementary_table13_sv_distribution.py) | Make supplementary table13 sv distribution |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie.py`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie.py) | Plot Msa svtype sorted bar pansv pie |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie.R) | Plot Msa svtype sorted bar pansv pie |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie_split_v4.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_Msa_svtype_sorted_bar_pansv_pie_split_v4.R) | Plot Msa svtype sorted bar pansv pie split v4 |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exact_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exact_exclude_Msa.R) | Plot pan core SV accumulation exact exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exact_exclude_Msa_split_v2.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exact_exclude_Msa_split_v2.R) | Plot pan core SV accumulation exact exclude Msa split v2 |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_pan_core_SV_accumulation_exclude_Msa.R) | Plot pan core SV accumulation exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_panSV_DEL_INS_length_distribution_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_panSV_DEL_INS_length_distribution_exclude_Msa.R) | Plot panSV DEL INS length distribution exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size.R) | Plot PAV detection frequency by size |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size_DEL_INS.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size_DEL_INS.R) | Extended Data Fig. 6b (revised): PAV detection frequency by size, DEL/INS only. |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size_DEL_INS_figure.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_detection_frequency_by_size_DEL_INS_figure.R) | Extended Data Fig. 6b (revised) figure, drawn from the count table written by |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_gene_context_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_gene_context_exclude_Msa.R) | Plot PAV gene context exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_gene_structure_per_species_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_gene_structure_per_species_exclude_Msa.R) | Plot PAV gene structure per species exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_shared_specific_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_shared_specific_exclude_Msa.R) | Plot PAV shared specific exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_association_4class_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_association_4class_exclude_Msa.R) | Plot PAV TE association 4class exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_association_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_association_exclude_Msa.R) | Plot PAV TE association exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_gene_distance_exclude_Msa.R`](07.read_mapping_flankQC_panSV_Msa/scripts/plot_PAV_TE_gene_distance_exclude_Msa.R) | Plot PAV TE gene distance exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/run_Msa_flankQC_panSV.sh`](07.read_mapping_flankQC_panSV_Msa/scripts/run_Msa_flankQC_panSV.sh) | Run Msa flankQC panSV |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/sort_vcf.py`](07.read_mapping_flankQC_panSV_Msa/scripts/sort_vcf.py) | Sort vcf |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_pan_core_SV_accumulation_exact_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_pan_core_SV_accumulation_exact_exclude_Msa.py) | Stat pan core SV accumulation exact exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_pan_core_SV_accumulation_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_pan_core_SV_accumulation_exclude_Msa.py) | Stat pan core SV accumulation exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_panSV_DEL_INS_length_distribution_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_panSV_DEL_INS_length_distribution_exclude_Msa.py) | Stat panSV DEL INS length distribution exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_gene_context_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_gene_context_exclude_Msa.py) | Stat PAV gene context exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_gene_structure_per_species_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_gene_structure_per_species_exclude_Msa.py) | Stat PAV gene structure per species exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_shared_specific_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_shared_specific_exclude_Msa.py) | Stat PAV shared specific exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_TE_association_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_TE_association_exclude_Msa.py) | Stat PAV TE association exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_TE_gene_distance_exclude_Msa.py`](07.read_mapping_flankQC_panSV_Msa/scripts/stat_PAV_TE_gene_distance_exclude_Msa.py) | Stat PAV TE gene distance exclude Msa |
| [`07.read_mapping_flankQC_panSV_Msa/scripts/vcf_filter_msa_with_tra.py`](07.read_mapping_flankQC_panSV_Msa/scripts/vcf_filter_msa_with_tra.py) |  |
| [`08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_mixed_illumina_expression_sv_distance.R`](08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_mixed_illumina_expression_sv_distance.R) | Plot mixed illumina expression sv distance |
| [`08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_pav_gene_structure_expression_divergence.R`](08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_pav_gene_structure_expression_divergence.R) | Plot pav gene structure expression divergence |
| [`08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_pav_gene_structure_fc_percent.R`](08_mixed_illumina_expression_SV_distance_20260705/scripts/plot_pav_gene_structure_fc_percent.R) | Plot pav gene structure fc percent |
| [`08_mixed_illumina_expression_SV_distance_20260705/scripts/run_mixed_illumina_expression_sv_distance.sh`](08_mixed_illumina_expression_SV_distance_20260705/scripts/run_mixed_illumina_expression_sv_distance.sh) | Run mixed illumina expression sv distance |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_expression_level_by_sv.R`](09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_expression_level_by_sv.R) | Plot strict expression level by sv |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_ortholog_expression.R`](09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_ortholog_expression.R) | Plot strict ortholog expression |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_pearson_distance_by_sv.R`](09_strict_ortholog_expression_SV_distance_20260705/scripts/plot_strict_pearson_distance_by_sv.R) | Plot strict pearson distance by sv |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/run_strict_ortholog_expression_sv.sh`](09_strict_ortholog_expression_SV_distance_20260705/scripts/run_strict_ortholog_expression_sv.sh) | Run strict ortholog expression sv |
| [`09_strict_ortholog_expression_SV_distance_20260705/scripts/strict_ortholog_expression_sv.py`](09_strict_ortholog_expression_SV_distance_20260705/scripts/strict_ortholog_expression_sv.py) |  |
| [`10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705/scripts/plot_read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson.R`](10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705/scripts/plot_read_mapping_primary_exon_CDS_UTR_promoter2kb_pearson.R) | Plot read mapping primary exon CDS UTR promoter2kb pearson |
| [`10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705/scripts/plot_read_mapping_primary_gene_structure_pearson.R`](10_gene_body_updown2kb_svlen50_tpmgt1_expression_20260705/scripts/plot_read_mapping_primary_gene_structure_pearson.R) | Plot read mapping primary gene structure pearson |
