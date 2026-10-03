# centromere

着丝粒分析：CENH3 ChIP-seq 定义功能着丝粒、卫星重复、甲基化、跨物种着丝粒共线性与演化、着丝粒 TE/LTR

| 脚本 | 说明 |
|---|---|
| [`my_run/h1.run_fastp_qc.py`](my_run/h1.run_fastp_qc.py) | fastp 质控和接头去除, fastp 的 -q 是“合格碱基的 Phred 阈值 |
| [`my_run/h2.bowtie2_to_ref.py`](my_run/h2.bowtie2_to_ref.py) | 过滤后的Chip seq数据比对到参考基因组 |
| [`my_run/h3.make_methylation_density_from_CX.py`](my_run/h3.make_methylation_density_from_CX.py) | 窗口甲基化统计脚本 |
| [`output_all/00_genome/h1.test_hifi_bam_for_cpG.sh`](output_all/00_genome/h1.test_hifi_bam_for_cpG.sh) | 单元测试：hifi bam for cpG（据文件名） |
| [`output_all/00_run0.sh`](output_all/00_run0.sh) | 命令：`SP=genome_Msa` |
| [`output_all/02_chipseq/genome_474/06_chipseq_peak/h1.filter_macsPeak.sh`](output_all/02_chipseq/genome_474/06_chipseq_peak/h1.filter_macsPeak.sh) | 过滤：macsPeak（据文件名） |
| [`output_all/03_repeat/genome_474.TRASH/h1.refind_monomer.sh`](output_all/03_repeat/genome_474.TRASH/h1.refind_monomer.sh) | 命令：`TRASH_BED=genome_474.TRASH_arrays.bed` |
| [`output_all/03_repeat/genome_474.TRASH/h2.refind_monomer2.sh`](output_all/03_repeat/genome_474.TRASH/h2.refind_monomer2.sh) | 先定义Chr6 CENH3 核心区 |
| [`output_all/03_repeat/genome_474.TRASH/h3.identity_chr_gap.sh`](output_all/03_repeat/genome_474.TRASH/h3.identity_chr_gap.sh) | 检查Chr6主峰区域是否组装完整 |
| [`output_all/03_repeat/genome_474.TRASH/h4.trf.sh`](output_all/03_repeat/genome_474.TRASH/h4.trf.sh) | 命令：`DAT=Chr6_93_100.5Mb.fa.2.7.7.80.10.50.2000.dat` |
| [`output_all/03_repeat/genome_R108.TRASH/h1.array_consensus2fa.py`](output_all/03_repeat/genome_R108.TRASH/h1.array_consensus2fa.py) |  |
| [`output_all/03_repeat/genome_R108.TRASH/h2.get_majority_consensus.py`](output_all/03_repeat/genome_R108.TRASH/h2.get_majority_consensus.py) | 获取：majority consensus（据文件名） |
| [`output_all/05_centromere_define/genome_474.satellite_validation/make_consensus_from_alignment.py`](output_all/05_centromere_define/genome_474.satellite_validation/make_consensus_from_alignment.py) | 生成：consensus from alignment（据文件名） |
| [`output_all/05_centromere_define/genome_474.TRF/trf_raw/run.sh`](output_all/05_centromere_define/genome_474.TRF/trf_raw/run.sh) | 运行（据文件名） |
| [`output_all/05_centromere_define/genome_Msa.satellite_validation/CenSat583_584/h1.blast2genome.sh`](output_all/05_centromere_define/genome_Msa.satellite_validation/CenSat583_584/h1.blast2genome.sh) | 命令：`SP=genome_Msa` |
| [`output_all/05_centromere_define/genome_Msa.satellite_validation/CenSat583_584_215m/h1.blast2genome.sh`](output_all/05_centromere_define/genome_Msa.satellite_validation/CenSat583_584_215m/h1.blast2genome.sh) | 命令：`SP=genome_Msa` |
| [`output_all/05_centromere_define/write_group.sh`](output_all/05_centromere_define/write_group.sh) | 写出：group（据文件名） |
| [`output_all/06_comparative/genome_474/h7.pyGenomeTracks_plot.sh`](output_all/06_comparative/genome_474/h7.pyGenomeTracks_plot.sh) | 命令：`SP=genome_474` |
| [`output_all/06_comparative/genome_474/tracks.ini`](output_all/06_comparative/genome_474/tracks.ini) |  |
| [`output_all/06_comparative/genome_A17/h5.cool2raw_cool.sh`](output_all/06_comparative/genome_A17/h5.cool2raw_cool.sh) | 命令：`cp genome_A17.40000.cool genome_A17.40000.raw.cool` |
| [`output_all/06_comparative/genome_Mpo/h1.sh`](output_all/06_comparative/genome_Mpo/h1.sh) | 命令：`CHR=Chr1` |
| [`output_all/06_comparative/genome_Mpo/tracks_2.ini`](output_all/06_comparative/genome_Mpo/tracks_2.ini) |  |
| [`output_all/06_comparative/genome_R108/define_pericentromere_local/h1.prepare_pericentromere_only_tracks.sh`](output_all/06_comparative/genome_R108/define_pericentromere_local/h1.prepare_pericentromere_only_tracks.sh) | 准备输入：pericentromere only tracks（据文件名） |
| [`output_all/06_comparative/genome_R108/define_pericentromere_local/h6.plot_each_chr.sh`](output_all/06_comparative/genome_R108/define_pericentromere_local/h6.plot_each_chr.sh) | 绘图：each chr（据文件名） |
| [`output_all/06_comparative/genome_R108/define_pericentromere_local/tracks.ini`](output_all/06_comparative/genome_R108/define_pericentromere_local/tracks.ini) |  |
| [`output_all/06_comparative/genome_R108/h1.centromere_region2track.sh`](output_all/06_comparative/genome_R108/h1.centromere_region2track.sh) | 命令：`python - <<'PY'` |
| [`output_all/06_comparative/genome_R108/h2.exam_matrix_value.sh`](output_all/06_comparative/genome_R108/h2.exam_matrix_value.sh) | 命令：`python - <<'PY'` |
| [`output_all/06_comparative/genome_R108/h3.exam_raw_balanced.sh`](output_all/06_comparative/genome_R108/h3.exam_raw_balanced.sh) | 命令：`python - <<'PY'` |
| [`output_all/06_comparative/genome_R108/h4.exam_NaN_bin.sh`](output_all/06_comparative/genome_R108/h4.exam_NaN_bin.sh) | 命令：`python - <<'PY'` |
| [`output_all/06_comparative/genome_R108/h5.cool2raw_cool.sh`](output_all/06_comparative/genome_R108/h5.cool2raw_cool.sh) | 命令：`cp genome_R108.40000.cool genome_R108.40000.raw.cool` |
| [`output_all/07_methylation/genome_474.hifi_methylation/weighted_methylation_window.py`](output_all/07_methylation/genome_474.hifi_methylation/weighted_methylation_window.py) |  |
| [`output_all/a1.sh`](output_all/a1.sh) | 只画主峰附近的trf repeat |
| [`output_all/a2.sh`](output_all/a2.sh) | 只画主峰附近的trf repeat |
| [`output_all/a3.sh`](output_all/a3.sh) | 只画主峰附近的trf repeat |
| [`output_all/a4.sh`](output_all/a4.sh) | 只画主峰附近的trf repeat |
| [`output_all/a5.sh`](output_all/a5.sh) | 只画主峰附近的trf repeat |
| [`output_all/a6.sh`](output_all/a6.sh) | 只画主峰附近的trf repeat |
| [`output_all/a7.sh`](output_all/a7.sh) | 只画主峰附近的trf repeat |
| [`output_all/a8.sh`](output_all/a8.sh) | 只画主峰附近的trf repeat |
| [`output_all/genome_474.summary_to_bed.py`](output_all/genome_474.summary_to_bed.py) |  |
| [`output_all/h1.run0.sh`](output_all/h1.run0.sh) | 命令：`SP=genome_Msa` |
| [`output_all/h2.run0.sh`](output_all/h2.run0.sh) | 命令：`SP=genome_Msa` |
| [`output_all/h3.run0.sh`](output_all/h3.run0.sh) | 命令：`SP=genome_474` |
| [`output_all/h4.run0.sh`](output_all/h4.run0.sh) | 让背景更干净，可以把低于某个 log2 值的小波动设为 0。比如低于 0.5 的信号全部隐藏 |
| [`output_all/h5.run0.sh`](output_all/h5.run0.sh) | 添加trf repeat 轨道 |
| [`output_all/h6.run0.sh`](output_all/h6.run0.sh) | 只画主峰附近的trf repeat |
| [`output_all/hhy.sh`](output_all/hhy.sh) | 命令：`sh ./a1.sh` |
| [`output_all/run0.sh`](output_all/run0.sh) | 命令：`SP=genome_Msa` |
| [`output_all/run1.sh`](output_all/run1.sh) | 环境CENH3_env |
| [`output_all/run2.sh`](output_all/run2.sh) | 命令：`SP=genome_474` |
| [`output_all/run3.sh`](output_all/run3.sh) | 命令：`SP=genome_474` |
| [`output_all/run4.sh`](output_all/run4.sh) | 命令：`SP=genome_474` |
| [`output_all/run5.sh`](output_all/run5.sh) | 命令：`SP=genome_474` |
| [`output_all/run6.sh`](output_all/run6.sh) | 命令：`SP=genome_474` |
| [`output_all/run7.sh`](output_all/run7.sh) | 命令：`SP=genome_474` |
| [`output_all/run8.py`](output_all/run8.py) |  |
| [`output_all/run9.py`](output_all/run9.py) |  |
| [`output_all/run10.sh`](output_all/run10.sh) | 命令：`SP=genome_474` |
| [`output_all/run11.sh`](output_all/run11.sh) | 指定候选单体生成一致性序列 |
| [`output_all/run12.sh`](output_all/run12.sh) | 提取主要单体的bed文件 |
| [`output_all/run13.sh`](output_all/run13.sh) | make_R108_density_tracks.sh |
| [`output_all/run14.sh`](output_all/run14.sh) | 计算甲基化水平 |
| [`output_all/run15.sh`](output_all/run15.sh) | 做成窗口化的Methylation density |
| [`output_all/run16.sh`](output_all/run16.sh) | 使用HIC数据分析A/B区室和TAD |
| [`output_all/run17.sh`](output_all/run17.sh) | 使用HIC数据分析A/B区室和TAD |
| [`output_all/run18.sh`](output_all/run18.sh) | 使用HIC数据分析A/B区室和TAD |
| [`output_all/run19.sh`](output_all/run19.sh) | 从 TRASH GFF 转成 array BED |
| [`output_all/run20.sh`](output_all/run20.sh) | 生成CENHH3功能着丝粒坐标 |
| [`output_all/run21.sh`](output_all/run21.sh) | 直接从 bigWig 里提取高信号区间，再合并成 CENH3 functional region |
| [`output_key/00_inventory/scripts/inventory_centromere_validation_inputs.py`](output_key/00_inventory/scripts/inventory_centromere_validation_inputs.py) | Inventory inputs for multi-species Mpo centromere-fate validation. |
| [`output_key/01_define_reference_CEN1_CEN8/scripts/project_reference_cens_to_mpo.py`](output_key/01_define_reference_CEN1_CEN8/scripts/project_reference_cens_to_mpo.py) |  |
| [`output_key/02_project_A17_R108_CENs_to_Mpo/scripts/project_cens_by_sequence_blocks.py`](output_key/02_project_A17_R108_CENs_to_Mpo/scripts/project_cens_by_sequence_blocks.py) |  |
| [`output_key/02_project_A17_R108_CENs_to_Mpo/scripts/project_reference_cens_to_mpo.py`](output_key/02_project_A17_R108_CENs_to_Mpo/scripts/project_reference_cens_to_mpo.py) |  |
| [`output_key/03_project_other_x8_CENs_to_Mpo/project_A17_cens_to_mpo_via_R108.py`](output_key/03_project_other_x8_CENs_to_Mpo/project_A17_cens_to_mpo_via_R108.py) |  |
| [`output_key/03_project_other_x8_CENs_to_Mpo/project_A17_to_R108_by_sequence_blocks.py`](output_key/03_project_other_x8_CENs_to_Mpo/project_A17_to_R108_by_sequence_blocks.py) |  |
| [`output_key/03_project_other_x8_CENs_to_Mpo/project_Msa_Mpo_reciprocal_by_sequence_blocks.py`](output_key/03_project_other_x8_CENs_to_Mpo/project_Msa_Mpo_reciprocal_by_sequence_blocks.py) |  |
| [`output_key/03_project_other_x8_CENs_to_Mpo/project_Msa_Mpo_reciprocal_curated_cens.py`](output_key/03_project_other_x8_CENs_to_Mpo/project_Msa_Mpo_reciprocal_curated_cens.py) | Project curated Msa CEN domains to Mpo and Mpo active CENs back to Msa. |
| [`output_key/03_project_other_x8_CENs_to_Mpo/project_x8_active_cens_to_R108.py`](output_key/03_project_other_x8_CENs_to_Mpo/project_x8_active_cens_to_R108.py) | Project active CENH3 cores from x=8 species to R108 ancestral CEN1-CEN8. |
| [`output_key/03_project_other_x8_CENs_to_Mpo/project_x8_active_cens_to_R108_candidates.py`](output_key/03_project_other_x8_CENs_to_Mpo/project_x8_active_cens_to_R108_candidates.py) | Long-format x=8 active CENH3 core projections to R108 CEN candidates. |
| [`output_key/04_reciprocal_project_Mpo_CENs_to_x8/project_Msa_Mpo_reciprocal_by_sequence_blocks.py`](output_key/04_reciprocal_project_Mpo_CENs_to_x8/project_Msa_Mpo_reciprocal_by_sequence_blocks.py) |  |
| [`output_key/04_reciprocal_project_Mpo_CENs_to_x8/reciprocal_mpo_cens_to_r108.py`](output_key/04_reciprocal_project_Mpo_CENs_to_x8/reciprocal_mpo_cens_to_r108.py) |  |
| [`output_key/05_fusion_chr_CEN_fate/cluster_projection_blocks.py`](output_key/05_fusion_chr_CEN_fate/cluster_projection_blocks.py) |  |
| [`output_key/05_fusion_chr_CEN_fate/summarize_centromere_fate.py`](output_key/05_fusion_chr_CEN_fate/summarize_centromere_fate.py) | 汇总统计：centromere fate（据文件名） |
| [`output_key/05_fusion_chr_CEN_fate/summarize_mpo_chr4_fusion_context.py`](output_key/05_fusion_chr_CEN_fate/summarize_mpo_chr4_fusion_context.py) | Summarize the Mpo Chr4 CEN4/CEN5 centromere-fate context. |
| [`output_key/05_fusion_chr_CEN_fate/summarize_mpo_chr4_local_synteny_context.py`](output_key/05_fusion_chr_CEN_fate/summarize_mpo_chr4_local_synteny_context.py) | Summarize local R108->Mpo sequence-block context around Mpo Chr4 CEN5 remnant. |
| [`output_key/06_repeat_gene_mappability_controls/measure_mpo_chip_input_coverage_control.py`](output_key/06_repeat_gene_mappability_controls/measure_mpo_chip_input_coverage_control.py) | Measure CENH3/Input coverage controls for the Mpo CEN5 candidate interval. |
| [`output_key/06_repeat_gene_mappability_controls/measure_mpo_interval_features.py`](output_key/06_repeat_gene_mappability_controls/measure_mpo_interval_features.py) |  |
| [`output_key/06_repeat_gene_mappability_controls/measure_mpo_orthogonal_chromatin_features.py`](output_key/06_repeat_gene_mappability_controls/measure_mpo_orthogonal_chromatin_features.py) | Measure orthogonal chromatin features around the Mpo CEN5 candidate. |
| [`output_key/07_publication_figures/plot_mpo_chr4_cen5_candidate.py`](output_key/07_publication_figures/plot_mpo_chr4_cen5_candidate.py) | 绘图：mpo chr4 cen5 candidate（据文件名） |
| [`output_key/08_multispecies_validation/scripts/plot_multispecies_CEN5_validation_matrix.py`](output_key/08_multispecies_validation/scripts/plot_multispecies_CEN5_validation_matrix.py) | Plot a compact multispecies CEN5 validation matrix. |
| [`output_key/08_multispecies_validation/scripts/summarize_multispecies_CEN5_validation.py`](output_key/08_multispecies_validation/scripts/summarize_multispecies_CEN5_validation.py) | Summarize multispecies evidence for the Mpo ancestral CEN5 fate story. |
| [`output_key/09_karyotype_model/scripts/plot_mpo_tad_based_karyotype_model.py`](output_key/09_karyotype_model/scripts/plot_mpo_tad_based_karyotype_model.py) | 绘图：mpo tad based karyotype model（据文件名） |
| [`output_key/09_karyotype_model/scripts/rename_rendered_karyotype_pngs.py`](output_key/09_karyotype_model/scripts/rename_rendered_karyotype_pngs.py) | 重命名：rendered karyotype pngs（据文件名） |
| [`output_key/09_karyotype_model/scripts/render_karyotype_pdfs.py`](output_key/09_karyotype_model/scripts/render_karyotype_pdfs.py) | Render local karyotype/chromosome-painting PDFs to PNG previews. |
| [`output_key/10_jcvi_mpo_synteny/summarize_jcvi_mpo_synteny.py`](output_key/10_jcvi_mpo_synteny/summarize_jcvi_mpo_synteny.py) | 汇总统计：jcvi mpo synteny（据文件名） |
| [`output_key/12_recall_Mpo_CENH3/scripts/call_mpo_cenh3_domains.py`](output_key/12_recall_Mpo_CENH3/scripts/call_mpo_cenh3_domains.py) |  |
| [`output_key/12_recall_Mpo_CENH3/scripts/plot_mpo_cenh3_recall_signal.py`](output_key/12_recall_Mpo_CENH3/scripts/plot_mpo_cenh3_recall_signal.py) | 绘图：mpo cenh3 recall signal（据文件名） |
| [`output_key/15_cen5_relic_raw_rerun/scripts/analyze_cen5_raw_rerun.py`](output_key/15_cen5_relic_raw_rerun/scripts/analyze_cen5_raw_rerun.py) | 分析：cen5 raw rerun（据文件名） |
| [`output_key/15_cen5_relic_raw_rerun/scripts/check_python_plot_envs.sh`](output_key/15_cen5_relic_raw_rerun/scripts/check_python_plot_envs.sh) | 检查：python plot envs（据文件名） |
| [`output_key/15_cen5_relic_raw_rerun/scripts/extract_raw_breakpoint_blocks.py`](output_key/15_cen5_relic_raw_rerun/scripts/extract_raw_breakpoint_blocks.py) | 提取：raw breakpoint blocks（据文件名） |
| [`output_key/15_cen5_relic_raw_rerun/scripts/plot_cen5_relic_raw_rerun.py`](output_key/15_cen5_relic_raw_rerun/scripts/plot_cen5_relic_raw_rerun.py) | 绘图：cen5 relic raw rerun（据文件名） |
| [`output_key/15_cen5_relic_raw_rerun/scripts/run_cen5_relic_raw_rerun.sh`](output_key/15_cen5_relic_raw_rerun/scripts/run_cen5_relic_raw_rerun.sh) | 运行：cen5 relic raw rerun（据文件名） |
| [`output_key/15_recall_Msa_474_CENH3/scripts/call_cenh3_domains_generic.py`](output_key/15_recall_Msa_474_CENH3/scripts/call_cenh3_domains_generic.py) |  |
| [`output_key/15_recall_Msa_474_CENH3/scripts/curate_cenh3_consensus.py`](output_key/15_recall_Msa_474_CENH3/scripts/curate_cenh3_consensus.py) |  |
| [`output_key/15_recall_Msa_474_CENH3/scripts/plot_cenh3_recall_signal_generic.py`](output_key/15_recall_Msa_474_CENH3/scripts/plot_cenh3_recall_signal_generic.py) | 绘图：cenh3 recall signal generic（据文件名） |
| [`output_key/15_recall_Msa_474_CENH3/scripts/run_recall_Msa_474_CENH3.sh`](output_key/15_recall_Msa_474_CENH3/scripts/run_recall_Msa_474_CENH3.sh) | 运行：recall Msa 474 CENH3（据文件名） |
| [`output_key/15_upload_tmp/analyze_cen5_raw_rerun.py`](output_key/15_upload_tmp/analyze_cen5_raw_rerun.py) | 分析：cen5 raw rerun（据文件名） |
| [`output_key/15_upload_tmp/run_cen5_relic_raw_rerun.sh`](output_key/15_upload_tmp/run_cen5_relic_raw_rerun.sh) | 运行：cen5 relic raw rerun（据文件名） |
| [`output_key/16_cen5_strong_validation/scripts/analyze_strong_validation.py`](output_key/16_cen5_strong_validation/scripts/analyze_strong_validation.py) | 分析：strong validation（据文件名） |
| [`output_key/16_cen5_strong_validation/scripts/check_inputs_tools.sh`](output_key/16_cen5_strong_validation/scripts/check_inputs_tools.sh) | 检查：inputs tools（据文件名） |
| [`output_key/16_cen5_strong_validation/scripts/hourly_monitor_strong_validation.sh`](output_key/16_cen5_strong_validation/scripts/hourly_monitor_strong_validation.sh) | 命令：`RUN=16_cen5_strong_validation` |
| [`output_key/16_cen5_strong_validation/scripts/make_strong_validation_figure.py`](output_key/16_cen5_strong_validation/scripts/make_strong_validation_figure.py) | 生成：strong validation figure（据文件名） |
| [`output_key/16_cen5_strong_validation/scripts/run_strong_validation.sh`](output_key/16_cen5_strong_validation/scripts/run_strong_validation.sh) | 运行：strong validation（据文件名） |
| [`output_key/16_CENH3_multitrack_evidence/pygenometracks_style/genome_474/tracks/Chr1.local/genome_474.Chr1.local.tracks.ini`](output_key/16_CENH3_multitrack_evidence/pygenometracks_style/genome_474/tracks/Chr1.local/genome_474.Chr1.local.tracks.ini) |  |
| [`output_key/16_CENH3_multitrack_evidence/pygenometracks_style/genome_Msa/tracks/Chr5.local/genome_Msa.Chr5.local.tracks.ini`](output_key/16_CENH3_multitrack_evidence/pygenometracks_style/genome_Msa/tracks/Chr5.local/genome_Msa.Chr5.local.tracks.ini) |  |
| [`output_key/16_CENH3_multitrack_evidence/scripts/build_cenh3_multitrack_evidence.py`](output_key/16_CENH3_multitrack_evidence/scripts/build_cenh3_multitrack_evidence.py) | 构建：cenh3 multitrack evidence（据文件名） |
| [`output_key/16_CENH3_multitrack_evidence/scripts/plot_brapa_style_cenh3_tracks.py`](output_key/16_CENH3_multitrack_evidence/scripts/plot_brapa_style_cenh3_tracks.py) | 绘图：brapa style cenh3 tracks（据文件名） |
| [`output_key/16_CENH3_multitrack_evidence/scripts/plot_publication_style_cenh3_tracks.py`](output_key/16_CENH3_multitrack_evidence/scripts/plot_publication_style_cenh3_tracks.py) | 绘图：publication style cenh3 tracks（据文件名） |
| [`output_key/16_CENH3_multitrack_evidence/scripts/plot_pygenometracks_cenh3_tracks.py`](output_key/16_CENH3_multitrack_evidence/scripts/plot_pygenometracks_cenh3_tracks.py) | 绘图：pygenometracks cenh3 tracks（据文件名） |
| [`output_key/17_main_figure_centromere_evolution.py`](output_key/17_main_figure_centromere_evolution.py) |  |
| [`output_key/17_main_figure_centromere_evolution/scripts/make_main_figure_centromere_evolution.py`](output_key/17_main_figure_centromere_evolution/scripts/make_main_figure_centromere_evolution.py) | 生成：main figure centromere evolution（据文件名） |
| [`output_key/17_main_figure_centromere_evolution/scripts/make_mpo_rta_cen5_model.py`](output_key/17_main_figure_centromere_evolution/scripts/make_mpo_rta_cen5_model.py) | 生成：mpo rta cen5 model（据文件名） |
| [`output_key/18_transition_repeat_intersect/scripts/analyze_transition_repeat_intersect.py`](output_key/18_transition_repeat_intersect/scripts/analyze_transition_repeat_intersect.py) | 分析：transition repeat intersect（据文件名） |
| [`output_key/18_transition_repeat_intersect/scripts/run_transition_repeat_intersect.sh`](output_key/18_transition_repeat_intersect/scripts/run_transition_repeat_intersect.sh) | 运行：transition repeat intersect（据文件名） |
| [`output_key/19_ltr_age_cen5_relics/scripts/analyze_ltr_age_cen5_relics.py`](output_key/19_ltr_age_cen5_relics/scripts/analyze_ltr_age_cen5_relics.py) | 分析：ltr age cen5 relics（据文件名） |
| [`output_key/20_ed5_cen5_relic_crossgenome/scripts/run_ed5_cen5_crossgenome.sh`](output_key/20_ed5_cen5_relic_crossgenome/scripts/run_ed5_cen5_crossgenome.sh) | 运行：ed5 cen5 crossgenome（据文件名） |
| [`output_key/20_ed5_cen5_relic_crossgenome/scripts/summarize_ed5_cen5_crossgenome.py`](output_key/20_ed5_cen5_relic_crossgenome/scripts/summarize_ed5_cen5_crossgenome.py) | 汇总统计：ed5 cen5 crossgenome（据文件名） |
| [`output_key/cluster_projection_blocks.py`](output_key/cluster_projection_blocks.py) |  |
| [`output_key/config.yaml`](output_key/config.yaml) |  |
| [`output_key/measure_mpo_interval_features.py`](output_key/measure_mpo_interval_features.py) |  |
| [`output_key/plot_mpo_chr4_cen5_candidate.py`](output_key/plot_mpo_chr4_cen5_candidate.py) | 绘图：mpo chr4 cen5 candidate（据文件名） |
| [`output_key/project_A17_cens_to_mpo_via_R108.py`](output_key/project_A17_cens_to_mpo_via_R108.py) |  |
| [`output_key/project_A17_to_R108_by_sequence_blocks.py`](output_key/project_A17_to_R108_by_sequence_blocks.py) |  |
| [`output_key/project_cens_by_sequence_blocks.py`](output_key/project_cens_by_sequence_blocks.py) |  |
| [`output_key/project_Msa_Mpo_reciprocal_by_sequence_blocks.py`](output_key/project_Msa_Mpo_reciprocal_by_sequence_blocks.py) |  |
| [`output_key/project_reference_cens_to_mpo.py`](output_key/project_reference_cens_to_mpo.py) |  |
| [`output_key/reciprocal_mpo_cens_to_r108.py`](output_key/reciprocal_mpo_cens_to_r108.py) |  |
| [`output_key/run_jcvi_A17_R108_only.sh`](output_key/run_jcvi_A17_R108_only.sh) | 运行：jcvi A17 R108 only（据文件名） |
| [`output_key/run_jcvi_all_today_curve.sh`](output_key/run_jcvi_all_today_curve.sh) | 运行：jcvi all today curve（据文件名） |
| [`output_key/run_jcvi_mpo_direct_pairs.sh`](output_key/run_jcvi_mpo_direct_pairs.sh) | 运行：jcvi mpo direct pairs（据文件名） |
| [`output_key/run_jcvi_mpo_synteny_karyotype.sh`](output_key/run_jcvi_mpo_synteny_karyotype.sh) | 运行：jcvi mpo synteny karyotype（据文件名） |
| [`output_key/run_recall_mpo_cenh3.sh`](output_key/run_recall_mpo_cenh3.sh) | 运行：recall mpo cenh3（据文件名） |
| [`output_key/summarize_centromere_fate.py`](output_key/summarize_centromere_fate.py) | 汇总统计：centromere fate（据文件名） |
| [`output_synthen/01_JCVI/01.run.sh`](output_synthen/01_JCVI/01.run.sh) | 批量修复多个物种的gff文件，使得生成的bed id和cds id 对应上 |
| [`output_synthen/01_JCVI/02.run.sh`](output_synthen/01_JCVI/02.run.sh) | 运行（据文件名） |
| [`output_synthen/01_JCVI/A17_R108_synteny_lz10_run/h2.plot_Sa_Sm_chr_synteny_CENH3.py`](output_synthen/01_JCVI/A17_R108_synteny_lz10_run/h2.plot_Sa_Sm_chr_synteny_CENH3.py) | 绘图：Sa Sm chr synteny CENH3（据文件名） |
| [`output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/h1.anchors_to_gene_links.py`](output_synthen/01_JCVI/R108_Mpo_synteny_lz10_run/h1.anchors_to_gene_links.py) |  |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/axt_to_sequence_blocks.py`](output_synthen/02_cent_synteny/genome_A17_genome_R108/axt_to_sequence_blocks.py) |  |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/chain_to_sequence_blocks.py`](output_synthen/02_cent_synteny/genome_A17_genome_R108/chain_to_sequence_blocks.py) |  |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/diagnose_centromere_blocks.py`](output_synthen/02_cent_synteny/genome_A17_genome_R108/diagnose_centromere_blocks.py) | 诊断：centromere blocks（据文件名） |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/normalize_paftools_maf.py`](output_synthen/02_cent_synteny/genome_A17_genome_R108/normalize_paftools_maf.py) | 标准化：paftools maf（据文件名） |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/paf_to_sequence_blocks.py`](output_synthen/02_cent_synteny/genome_A17_genome_R108/paf_to_sequence_blocks.py) |  |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/plot_A17_R108_three_alignment_methods.sh`](output_synthen/02_cent_synteny/genome_A17_genome_R108/plot_A17_R108_three_alignment_methods.sh) | 绘图：A17 R108 three alignment methods（据文件名） |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/plot_sequence_synteny_CENH3.py`](output_synthen/02_cent_synteny/genome_A17_genome_R108/plot_sequence_synteny_CENH3.py) | 绘图：sequence synteny CENH3（据文件名） |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/run_A17_R108_chainnet_pipeline.sh`](output_synthen/02_cent_synteny/genome_A17_genome_R108/run_A17_R108_chainnet_pipeline.sh) | 运行：A17 R108 chainnet pipeline（据文件名） |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/run_A17_R108_sequence_synteny.sh`](output_synthen/02_cent_synteny/genome_A17_genome_R108/run_A17_R108_sequence_synteny.sh) | 运行：A17 R108 sequence synteny（据文件名） |
| [`output_synthen/02_cent_synteny/genome_A17_genome_R108/summarize_sequence_blocks.py`](output_synthen/02_cent_synteny/genome_A17_genome_R108/summarize_sequence_blocks.py) | 汇总统计：sequence blocks（据文件名） |
| [`output_synthen/03_cent_TE/genome_474/CENH3_functional_centromere/analyze_LTRRT_cent_noncent.py`](output_synthen/03_cent_TE/genome_474/CENH3_functional_centromere/analyze_LTRRT_cent_noncent.py) | 分析：LTRRT cent noncent（据文件名） |
| [`output_synthen/03_cent_TE/genome_474/pericentromere_6Mb/analyze_LTRRT_cent_noncent.py`](output_synthen/03_cent_TE/genome_474/pericentromere_6Mb/analyze_LTRRT_cent_noncent.py) | 分析：LTRRT cent noncent（据文件名） |
| [`output_synthen/03_cent_TE/genome_R108/analyze_LTRRT_cent_noncent.py`](output_synthen/03_cent_TE/genome_R108/analyze_LTRRT_cent_noncent.py) | 分析：LTRRT cent noncent（据文件名） |
| [`output_synthen/03_cent_TE/plot_combined_CENH3core_LTRRT_insertion_time.py`](output_synthen/03_cent_TE/plot_combined_CENH3core_LTRRT_insertion_time.py) | 绘图：combined CENH3core LTRRT insertion time（据文件名） |
| [`output_synthen/03_cent_TE/plot_combined_pericentromere6Mb_LTRRT_insertion_time.py`](output_synthen/03_cent_TE/plot_combined_pericentromere6Mb_LTRRT_insertion_time.py) | 绘图：combined pericentromere6Mb LTRRT insertion time（据文件名） |
| [`output_synthen/04_cent_TE_type/analyze_functional_centromere_TE_types.py`](output_synthen/04_cent_TE_type/analyze_functional_centromere_TE_types.py) | 分析：functional centromere TE types（据文件名） |
| [`output_synthen/04_cent_TE_type/TEsorter_all_species/analyze_centromere_TE_types_tesorter.py`](output_synthen/04_cent_TE_type/TEsorter_all_species/analyze_centromere_TE_types_tesorter.py) | 分析：centromere TE types tesorter（据文件名） |
| [`output_synthen/04_cent_TE_type/TEsorter_all_species/combine_centromere_TE_types_tesorter.py`](output_synthen/04_cent_TE_type/TEsorter_all_species/combine_centromere_TE_types_tesorter.py) |  |
| [`output_synthen/04_cent_TE_type/TEsorter_R108_test/analyze_centromere_TE_types_tesorter.py`](output_synthen/04_cent_TE_type/TEsorter_R108_test/analyze_centromere_TE_types_tesorter.py) | 分析：centromere TE types tesorter（据文件名） |
| [`output_synthen/05_functional_centromere_intact_LTR_similarity/build_functional_centromere_intact_LTR_similarity.py`](output_synthen/05_functional_centromere_intact_LTR_similarity/build_functional_centromere_intact_LTR_similarity.py) | 构建：functional centromere intact LTR similarity（据文件名） |
| [`output_synthen/05_functional_centromere_intact_LTR_similarity/plot_functional_centromere_intact_LTR_similarity.py`](output_synthen/05_functional_centromere_intact_LTR_similarity/plot_functional_centromere_intact_LTR_similarity.py) | 绘图：functional centromere intact LTR similarity（据文件名） |
| [`output_synthen/05_functional_centromere_intact_LTR_similarity/run_05_functional_centromere_intact_LTR_similarity.sh`](output_synthen/05_functional_centromere_intact_LTR_similarity/run_05_functional_centromere_intact_LTR_similarity.sh) | 运行：05 functional centromere intact LTR similarity（据文件名） |
| [`output_synthen/05_functional_centromere_intact_LTR_similarity/run_05_functional_centromere_intact_LTR_similarity_blastn.sh`](output_synthen/05_functional_centromere_intact_LTR_similarity/run_05_functional_centromere_intact_LTR_similarity_blastn.sh) | 运行：05 functional centromere intact LTR similarity blastn（据文件名） |
| [`output_synthen/05_functional_centromere_intact_LTR_similarity/summarize_05_blastn_similarity.py`](output_synthen/05_functional_centromere_intact_LTR_similarity/summarize_05_blastn_similarity.py) | 汇总统计：05 blastn similarity（据文件名） |
| [`output_synthen/06_functional_centromere_sequence_composition/all_species_TableS11_style/calculate_functional_centromere_composition_tesorter_parallel.py`](output_synthen/06_functional_centromere_sequence_composition/all_species_TableS11_style/calculate_functional_centromere_composition_tesorter_parallel.py) |  |
| [`output_synthen/06_functional_centromere_sequence_composition/calculate_functional_centromere_sequence_composition.py`](output_synthen/06_functional_centromere_sequence_composition/calculate_functional_centromere_sequence_composition.py) |  |
| [`output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/calculate_R108_functional_centromere_composition_tesorter_parallel.py`](output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/calculate_R108_functional_centromere_composition_tesorter_parallel.py) |  |
| [`output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/compare_R108_priority_modes.py`](output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/compare_R108_priority_modes.py) | 比较：R108 priority modes（据文件名） |
| [`output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/format_R108_TableS11_sequence_composition.py`](output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/format_R108_TableS11_sequence_composition.py) | 格式化：R108 TableS11 sequence composition（据文件名） |
| [`output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/summarize_R108_sequence_composition_test.py`](output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/summarize_R108_sequence_composition_test.py) | 汇总统计：R108 sequence composition test（据文件名） |
| [`output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/TE_first/calculate_R108_functional_centromere_composition_tesorter_parallel.py`](output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/TE_first/calculate_R108_functional_centromere_composition_tesorter_parallel.py) |  |
| [`output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/TE_first/run_06_R108_te_first_patch.sh`](output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/TE_first/run_06_R108_te_first_patch.sh) | 运行：06 R108 te first patch（据文件名） |
| [`output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/TE_first/summarize_R108_sequence_composition_test.py`](output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test/TE_first/summarize_R108_sequence_composition_test.py) | 汇总统计：R108 sequence composition test（据文件名） |
| [`output_synthen/07_intact_LTRRT_TableS8_summary/calculate_intact_LTRRT_TableS8_summary.py`](output_synthen/07_intact_LTRRT_TableS8_summary/calculate_intact_LTRRT_TableS8_summary.py) |  |
| [`output_synthen/07_intact_LTRRT_TableS8_summary/prepare_passlist_internal_sequences.py`](output_synthen/07_intact_LTRRT_TableS8_summary/prepare_passlist_internal_sequences.py) | 准备输入：passlist internal sequences（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/plot_selected_centromere_full_length_LTRRT_similarity.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/plot_selected_centromere_full_length_LTRRT_similarity.py) | 绘图：selected centromere full length LTRRT similarity（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/prepare_selected_centromere_full_length_LTRRT.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/prepare_selected_centromere_full_length_LTRRT.py) | 准备输入：selected centromere full length LTRRT（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/04_full_length_iqtree_MFP_order/plot_full_length_ltrrt_similarity_publication_style.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/04_full_length_iqtree_MFP_order/plot_full_length_ltrrt_similarity_publication_style.py) | Plot full-length intact LTR-RT similarity as a tree-ordered triangular heatmap. |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/04_full_length_iqtree_MFP_order/plot_similarity_by_iqtree_order.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/04_full_length_iqtree_MFP_order/plot_similarity_by_iqtree_order.py) | 绘图：similarity by iqtree order（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/compute_alignment_similarity_matrix.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/compute_alignment_similarity_matrix.py) |  |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/plot_clustered_similarity_matrix.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/plot_clustered_similarity_matrix.py) | 绘图：clustered similarity matrix（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/plot_similarity_by_rt_tree_order.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/plot_similarity_by_rt_tree_order.py) | 绘图：similarity by rt tree order（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/plot_tableS8_rt_full_length_similarity.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/plot_tableS8_rt_full_length_similarity.py) | 绘图：tableS8 rt full length similarity（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/prepare_selected_RT_domains.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/prepare_selected_RT_domains.py) | 准备输入：selected RT domains（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/prepare_tableS8_rt_domain_ltrrt.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/prepare_tableS8_rt_domain_ltrrt.py) | 准备输入：tableS8 rt domain ltrrt（据文件名） |
| [`output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/prepare_tableS8_rt_full_length_ltrrt.py`](output_synthen/08_selected_functional_centromere_LTRRT_similarity/two_similarity_methods/prepare_tableS8_rt_full_length_ltrrt.py) | 准备输入：tableS8 rt full length ltrrt（据文件名） |
