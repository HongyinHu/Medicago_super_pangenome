# gene_family_pangenome

OrthoFinder 基因家族与泛基因组分类、物种树（IQ-TREE）、分化时间（MCMCTree）、CAFE5、Ka/Ks

| 脚本 | 说明 |
|---|---|
| [`my_run/b1.get_low_copy_orthology.py`](my_run/b1.get_low_copy_orthology.py) | Orthogroup_Sequences中获取低拷贝同源基因集 |
| [`my_run/b2.run_mafft.pl`](my_run/b2.run_mafft.pl) | 运行：mafft（据文件名） |
| [`my_run/b3.run_trimal.py`](my_run/b3.run_trimal.py) | 1. 将gblock过滤后的多序列比对结果，进一步使用trimal进行过滤 |
| [`my_run/b4.pal2nal_pep2cds.py`](my_run/b4.pal2nal_pep2cds.py) | 将文件夹下的pep转为对应的cds序列 |
| [`my_run/b5.get_codon12.py`](my_run/b5.get_codon12.py) | 提取密码子的第1，第2位 |
| [`my_run/b6.get_concatence_fa.py`](my_run/b6.get_concatence_fa.py) | 合并cds序列 |
| [`my_run/b7.run_iqtree.py`](my_run/b7.run_iqtree.py) | 1. 执行iqtree的sh文件, 为supertree 使用 |
| [`my_run/b8.remain_species_name_from_genetree.py`](my_run/b8.remain_species_name_from_genetree.py) | 删除基因树中的基因，保留物种名称 |
| [`my_run/b9.filter_aln_fa.py`](my_run/b9.filter_aln_fa.py) | 过滤溯祖树比对序列 |
| [`my_run/cafe5/h1.get_rapid_change_fam.py`](my_run/cafe5/h1.get_rapid_change_fam.py) | 根据cafe5的输出结果，筛选某一个物种的快速扩张或快速收缩的基因家族 |
| [`my_run/cafe5/h2.from_fam_to_geneid.py`](my_run/cafe5/h2.from_fam_to_geneid.py) | 从获得基因组家族中查找具体的geneID |
| [`my_run/d1.split_seq.py`](my_run/d1.split_seq.py) | 对每个单拷贝基因家族进行处理，创建该家族单独文件夹，包含pep和cds序列 |
| [`my_run/d2.align.py`](my_run/d2.align.py) | 多序列比对并获取密码子序列 |
| [`my_run/d3.get_concatence.py`](my_run/d3.get_concatence.py) | 合并cds序列 |
| [`my_run/d4.extract_4Dsite.py`](my_run/d4.extract_4Dsite.py) | 提取四重简并为点从fasta文件中 |
| [`my_run/d5.orthofinder2cafe_input.py`](my_run/d5.orthofinder2cafe_input.py) | 将Orthogroups.GeneCount.tsv整理成cafe的输入格式 |
| [`my_run/h1.format_pep_prefixID_orthofinder_data.py`](my_run/h1.format_pep_prefixID_orthofinder_data.py) | 为系统发育树候选物种蛋白序列添加物种名前缀 |
| [`my_run/h1.remove_gene_tree.py`](my_run/h1.remove_gene_tree.py) | 删除基因树中的基因，保留物种名称 |
| [`my_run/h2.statistic_interproscan_anno.py`](my_run/h2.statistic_interproscan_anno.py) | 统计interproscan注释结果 |
| [`my_run/h3.format_Orthotsv_PanGP.py`](my_run/h3.format_Orthotsv_PanGP.py) | 修改格式Orthotsv适应PanGP软件 |
| [`my_run/h4.orthology_sort.py`](my_run/h4.orthology_sort.py) | 对同源基因家族进行分类 |
| [`my_run/h5.gen_cafe_input.py`](my_run/h5.gen_cafe_input.py) | 将Orthogroups.GeneCount.tsv文件修改为cafe输入文件格式 |
| [`my_run/h6.heatmap_pangenome_present_absent.py`](my_run/h6.heatmap_pangenome_present_absent.py) | 转换格式绘制家族存在缺失热图 |
| [`my_run/p1.plot_simluted_core_pan.py`](my_run/p1.plot_simluted_core_pan.py) | 绘图：simluted core pan（据文件名） |
| [`my_run/p2.plot_bar_frequency_family.py`](my_run/p2.plot_bar_frequency_family.py) | 绘制基因家族中物种数量分布条形图 |
| [`my_run/p3.plot_species_family.py`](my_run/p3.plot_species_family.py) | 绘制每个物种中的分类分布 |
| [`my_run/p4.plot_statistic_anno.py`](my_run/p4.plot_statistic_anno.py) | 绘制interpro注释结果 |
| [`my_run/p5.plot_kaks_res.py`](my_run/p5.plot_kaks_res.py) | 绘制不同泛基因组家族中kaks的比较箱线图 |
| [`my_run/s1.run_orthofinder.py`](my_run/s1.run_orthofinder.py) | 获取系统发育树使用的单拷贝基因家族 |
| [`my_run/s2.run_iqtree.py`](my_run/s2.run_iqtree.py) | 获取单拷贝基因家族序列后，进行系统发育建树 |
| [`my_run/s3.core_genes_class.py`](my_run/s3.core_genes_class.py) | 以下是一个Python脚本，可以从OrthoFinder结果中区分泛基因组中的四类基因（core, soft-core, shell和specific）： |
| [`my_run/s4.simluted_core_pan.py`](my_run/s4.simluted_core_pan.py) | 计算在随意组合情况下，pan和core的基因家族数量变化趋势 |
| [`my_run/s5.frequency_family_number.py`](my_run/s5.frequency_family_number.py) | 计算基因家族包含物种数量频率分布 |
| [`my_run/s6.run_interproscan_annotation.py`](my_run/s6.run_interproscan_annotation.py) | 对每个基因组蛋白序列运行interproscan进行interpro结构域注释 |
| [`my_run/s7.calcalator_ka_ks.py`](my_run/s7.calcalator_ka_ks.py) | 使用orthofinder聚类的同源基因家族簇，计算四类范基因组四类基因集的ka/ks |
| [`my_run/s8.kaks_revision_para.py`](my_run/s8.kaks_revision_para.py) | 产生两两组合的基因对，用于计算成对的Ka/KS |
| [`my_run/s9.run_ParaAT_kaks.py`](my_run/s9.run_ParaAT_kaks.py) | 使用ParaAT批量计算kaks |
| [`my_run/s10.get_res_kaks.py`](my_run/s10.get_res_kaks.py) | 提取ka/ks结果 |
| [`output/2.iqtree/genetrees/h1.remove_gene_tree.py`](output/2.iqtree/genetrees/h1.remove_gene_tree.py) | 删除基因树中的基因，保留物种名称 |
| [`output/2.iqtree/genetrees/s12.pep_gblocks_wrapper.py`](output/2.iqtree/genetrees/s12.pep_gblocks_wrapper.py) | assuming that all alignment files end with ".aln" |
| [`output/2.iqtree/genetrees/s13.gblock_trimal.py`](output/2.iqtree/genetrees/s13.gblock_trimal.py) | 1. 将gblock过滤后的多序列比对结果，进一步使用trimal进行过滤 |
| [`output/2.iqtree/genetrees/s14.run_supertree_iqtree.py`](output/2.iqtree/genetrees/s14.run_supertree_iqtree.py) | 1. 执行iqtree的sh文件, 为supertree 使用 |
| [`output/2.iqtree_cds/genetrees/gCF_sCF/h1.remove_gene_tree.py`](output/2.iqtree_cds/genetrees/gCF_sCF/h1.remove_gene_tree.py) | 删除基因树中的基因，保留物种名称 |
| [`output/2.iqtree_cds/genetrees/s13.gblock_trimal.py`](output/2.iqtree_cds/genetrees/s13.gblock_trimal.py) | 1. 将gblock过滤后的多序列比对结果，进一步使用trimal进行过滤 |
| [`output/2.iqtree_no_Gma/codon12_gCF_sCF/h1.remove_gene_tree.py`](output/2.iqtree_no_Gma/codon12_gCF_sCF/h1.remove_gene_tree.py) | 删除基因树中的基因，保留物种名称 |
| [`output/2.iqtree_no_Gma/concatence_gCF_sCF/CDS_gcf_scf/h1.remove_gene_tree.py`](output/2.iqtree_no_Gma/concatence_gCF_sCF/CDS_gcf_scf/h1.remove_gene_tree.py) | 删除基因树中的基因，保留物种名称 |
| [`output/2.iqtree_no_Gma/concatence_gCF_sCF/codon12_gcf_scf/h1.remove_gene_tree.py`](output/2.iqtree_no_Gma/concatence_gCF_sCF/codon12_gcf_scf/h1.remove_gene_tree.py) | 删除基因树中的基因，保留物种名称 |
| [`output/2.iqtree_no_Gma/get_concatence.py`](output/2.iqtree_no_Gma/get_concatence.py) | 合并cds序列 |
| [`output/2.iqtree_no_Gma/run_iqtree.py`](output/2.iqtree_no_Gma/run_iqtree.py) | 1. 执行iqtree的sh文件, 为supertree 使用 |
| [`output/2.iqtree_no_Gma/run_trimal.py`](output/2.iqtree_no_Gma/run_trimal.py) | 1. 将gblock过滤后的多序列比对结果，进一步使用trimal进行过滤 |
| [`output/2.iqtree_no_ZM4_outgroup/s14.run_supertree_iqtree.py`](output/2.iqtree_no_ZM4_outgroup/s14.run_supertree_iqtree.py) | 1. 执行iqtree的sh文件, 为supertree 使用 |
| [`output/3.interpro_anno/h1.get_pep_fa.py`](output/3.interpro_anno/h1.get_pep_fa.py) | 先获取泛基因组类别的基因列表，后根据基因列表获取序列 |
| [`output/3.interpro_anno/h2.get_pep_from_id.py`](output/3.interpro_anno/h2.get_pep_from_id.py) | 先获取泛基因组类别的基因列表，后根据基因列表获取序列 |
| [`output/5.PAML/cds_paml2/baseml/baseml.ctl`](output/5.PAML/cds_paml2/baseml/baseml.ctl) |  |
| [`output/5.PAML/cds_paml2/mcmctree2/mcmctree.ctl`](output/5.PAML/cds_paml2/mcmctree2/mcmctree.ctl) |  |
| [`output/5.PAML/cds_paml2/mcmctree/mcmctree.ctl`](output/5.PAML/cds_paml2/mcmctree/mcmctree.ctl) |  |
| [`output/5.PAML/cds_paml2/mcmctree/tmp0001.ctl`](output/5.PAML/cds_paml2/mcmctree/tmp0001.ctl) |  |
| [`output/5.PAML/cds_paml/baseml/baseml.ctl`](output/5.PAML/cds_paml/baseml/baseml.ctl) |  |
| [`output/5.PAML/cds_paml/data/4.4Dsites.pl`](output/5.PAML/cds_paml/data/4.4Dsites.pl) |  |
| [`output/5.PAML/cds_paml/data/fasta2phy.py`](output/5.PAML/cds_paml/data/fasta2phy.py) | fasta2phy.py |
| [`output/5.PAML/cds_paml/mcmctree2/mcmctree.ctl`](output/5.PAML/cds_paml/mcmctree2/mcmctree.ctl) |  |
| [`output/5.PAML/cds_paml/mcmctree/mcmctree.ctl`](output/5.PAML/cds_paml/mcmctree/mcmctree.ctl) |  |
| [`output/5.PAML/cds_paml/mcmctree/tmp0001.ctl`](output/5.PAML/cds_paml/mcmctree/tmp0001.ctl) |  |
| [`output/6.cafe/run_cafe.sh`](output/6.cafe/run_cafe.sh) | 运行：cafe（据文件名） |
| [`output/7.RNA_tpm/h1.get_geneTPM.py`](output/7.RNA_tpm/h1.get_geneTPM.py) | 1. 从stringtie ballgown 中产生的样品注释文件中提取转录本的表达量 |
| [`output/7.RNA_tpm/h2.get_panGP_tpm.py`](output/7.RNA_tpm/h2.get_panGP_tpm.py) | 获取每个物种在泛基因组中的各种类型的基因表达水平 |
| [`output/8.nucleotide_diversity/h1.get_genefam_gene.py`](output/8.nucleotide_diversity/h1.get_genefam_gene.py) | 获取对应基因家族的基因id |
| [`output/8.nucleotide_diversity/h2.cal_fam_pai.py`](output/8.nucleotide_diversity/h2.cal_fam_pai.py) | 计算泛基因组家族类型的核酸多样性 |
| [`output/8.nucleotide_diversity/nucleotide_diversity.R`](output/8.nucleotide_diversity/nucleotide_diversity.R) | 不同类型的泛基因组家族类型的核酸多样性计算 |
| [`output_M46_new/scripts/audit_sequence_characters.sh`](output_M46_new/scripts/audit_sequence_characters.sh) | 核查：sequence characters（据文件名） |
| [`output_M46_new/scripts/postprocess_pan_genome_B.R`](output_M46_new/scripts/postprocess_pan_genome_B.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_M46_new/scripts/postprocess_pan_genome_new.R`](output_M46_new/scripts/postprocess_pan_genome_new.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_M46_new/scripts/postprocess_table_only.R`](output_M46_new/scripts/postprocess_table_only.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_M46_new/scripts/postprocess_table_source.R`](output_M46_new/scripts/postprocess_table_source.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_M46_new/scripts/prepare_clean_input.sh`](output_M46_new/scripts/prepare_clean_input.sh) | 准备输入：clean input（据文件名） |
| [`output_M46_new/scripts/run_input_qc.sh`](output_M46_new/scripts/run_input_qc.sh) | 运行：input qc（据文件名） |
| [`output_M46_new/scripts/run_orthofinder.sh`](output_M46_new/scripts/run_orthofinder.sh) | 运行：orthofinder（据文件名） |
| [`output_pangenome_new/scripts/audit_sequence_characters.sh`](output_pangenome_new/scripts/audit_sequence_characters.sh) | 核查：sequence characters（据文件名） |
| [`output_pangenome_new/scripts/build_supplementary_workbook.mjs`](output_pangenome_new/scripts/build_supplementary_workbook.mjs) | 构建：supplementary workbook（据文件名） |
| [`output_pangenome_new/scripts/finalize_completed_status.sh`](output_pangenome_new/scripts/finalize_completed_status.sh) | 收尾汇总：completed status（据文件名） |
| [`output_pangenome_new/scripts/postprocess_pan_genome.R`](output_pangenome_new/scripts/postprocess_pan_genome.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_pangenome_new/scripts/prepare_clean_input.sh`](output_pangenome_new/scripts/prepare_clean_input.sh) | 准备输入：clean input（据文件名） |
| [`output_pangenome_new/scripts/qa_pan_genome_outputs.R`](output_pangenome_new/scripts/qa_pan_genome_outputs.R) |  |
| [`output_pangenome_new/scripts/run_input_qc.sh`](output_pangenome_new/scripts/run_input_qc.sh) | 运行：input qc（据文件名） |
| [`output_pangenome_new/scripts/run_orthofinder.sh`](output_pangenome_new/scripts/run_orthofinder.sh) | 运行：orthofinder（据文件名） |
