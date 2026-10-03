# gene_family_pangenome

OrthoFinder gene families and pangenome categories, species tree (IQ-TREE), divergence times (MCMCTree), CAFE5, Ka/Ks

| Script | Description |
|---|---|
| [`my_run/b1.get_low_copy_orthology.py`](my_run/b1.get_low_copy_orthology.py) | Extract low-copy orthogroups from Orthogroup_Sequences |
| [`my_run/b2.run_mafft.pl`](my_run/b2.run_mafft.pl) | Run mafft |
| [`my_run/b3.run_trimal.py`](my_run/b3.run_trimal.py) | Further filter Gblocks-trimmed multiple sequence alignments with trimAl |
| [`my_run/b4.pal2nal_pep2cds.py`](my_run/b4.pal2nal_pep2cds.py) | Convert protein alignments to codon (CDS) alignments |
| [`my_run/b5.get_codon12.py`](my_run/b5.get_codon12.py) | Extract codon positions 1 and 2 |
| [`my_run/b6.get_concatence_fa.py`](my_run/b6.get_concatence_fa.py) | Concatenate CDS alignments |
| [`my_run/b7.run_iqtree.py`](my_run/b7.run_iqtree.py) | Generate and run IQ-TREE commands for the supertree analysis |
| [`my_run/b8.remain_species_name_from_genetree.py`](my_run/b8.remain_species_name_from_genetree.py) | Strip gene IDs from gene trees, keeping species names |
| [`my_run/b9.filter_aln_fa.py`](my_run/b9.filter_aln_fa.py) | Filter alignments for coalescent (ASTRAL) trees |
| [`my_run/cafe5/h1.get_rapid_change_fam.py`](my_run/cafe5/h1.get_rapid_change_fam.py) | Select rapidly expanding/contracting families for a species from CAFE5 output |
| [`my_run/cafe5/h2.from_fam_to_geneid.py`](my_run/cafe5/h2.from_fam_to_geneid.py) | Retrieve gene IDs for selected gene families |
| [`my_run/d1.split_seq.py`](my_run/d1.split_seq.py) | Create a folder with protein and CDS sequences for each single-copy family |
| [`my_run/d2.align.py`](my_run/d2.align.py) | Multiple sequence alignment and codon alignment |
| [`my_run/d3.get_concatence.py`](my_run/d3.get_concatence.py) | Concatenate CDS alignments |
| [`my_run/d4.extract_4Dsite.py`](my_run/d4.extract_4Dsite.py) | Extract fourfold degenerate sites from FASTA |
| [`my_run/d5.orthofinder2cafe_input.py`](my_run/d5.orthofinder2cafe_input.py) | Convert Orthogroups.GeneCount.tsv to CAFE input |
| [`my_run/h1.format_pep_prefixID_orthofinder_data.py`](my_run/h1.format_pep_prefixID_orthofinder_data.py) | Add species prefixes to protein IDs for phylogenetic analysis |
| [`my_run/h1.remove_gene_tree.py`](my_run/h1.remove_gene_tree.py) | Strip gene IDs from gene trees, keeping species names |
| [`my_run/h2.statistic_interproscan_anno.py`](my_run/h2.statistic_interproscan_anno.py) | Summarize InterProScan annotation |
| [`my_run/h3.format_Orthotsv_PanGP.py`](my_run/h3.format_Orthotsv_PanGP.py) | Reformat Orthogroups.tsv for PanGP |
| [`my_run/h4.orthology_sort.py`](my_run/h4.orthology_sort.py) | Classify gene families |
| [`my_run/h5.gen_cafe_input.py`](my_run/h5.gen_cafe_input.py) | Convert Orthogroups.GeneCount.tsv to CAFE input |
| [`my_run/h6.heatmap_pangenome_present_absent.py`](my_run/h6.heatmap_pangenome_present_absent.py) | Prepare and plot gene-family presence/absence heatmap |
| [`my_run/p1.plot_simluted_core_pan.py`](my_run/p1.plot_simluted_core_pan.py) | Plot simluted core pan |
| [`my_run/p2.plot_bar_frequency_family.py`](my_run/p2.plot_bar_frequency_family.py) | Bar plot of species counts per gene family |
| [`my_run/p3.plot_species_family.py`](my_run/p3.plot_species_family.py) | Plot pangenome category composition per species |
| [`my_run/p4.plot_statistic_anno.py`](my_run/p4.plot_statistic_anno.py) | Plot InterPro annotation results |
| [`my_run/p5.plot_kaks_res.py`](my_run/p5.plot_kaks_res.py) | Boxplot of Ka/Ks across pangenome family categories |
| [`my_run/s1.run_orthofinder.py`](my_run/s1.run_orthofinder.py) | Run OrthoFinder to obtain single-copy families for phylogeny |
| [`my_run/s2.run_iqtree.py`](my_run/s2.run_iqtree.py) | Build phylogeny from single-copy gene families |
| [`my_run/s3.core_genes_class.py`](my_run/s3.core_genes_class.py) | Classify OrthoFinder gene families into core, soft-core, shell and specific categories |
| [`my_run/s4.simluted_core_pan.py`](my_run/s4.simluted_core_pan.py) | Simulate pan/core family counts over random genome combinations |
| [`my_run/s5.frequency_family_number.py`](my_run/s5.frequency_family_number.py) | Frequency distribution of species per gene family |
| [`my_run/s6.run_interproscan_annotation.py`](my_run/s6.run_interproscan_annotation.py) | Run InterProScan domain annotation on each proteome |
| [`my_run/s7.calcalator_ka_ks.py`](my_run/s7.calcalator_ka_ks.py) | Compute Ka/Ks for the four pangenome gene categories from OrthoFinder families |
| [`my_run/s8.kaks_revision_para.py`](my_run/s8.kaks_revision_para.py) | Generate pairwise gene pairs for Ka/Ks calculation |
| [`my_run/s9.run_ParaAT_kaks.py`](my_run/s9.run_ParaAT_kaks.py) | Batch Ka/Ks calculation with ParaAT |
| [`my_run/s10.get_res_kaks.py`](my_run/s10.get_res_kaks.py) | Extract Ka/Ks results |
| [`output/2.iqtree/genetrees/h1.remove_gene_tree.py`](output/2.iqtree/genetrees/h1.remove_gene_tree.py) | Strip gene IDs from gene trees, keeping species names |
| [`output/2.iqtree/genetrees/s12.pep_gblocks_wrapper.py`](output/2.iqtree/genetrees/s12.pep_gblocks_wrapper.py) | assuming that all alignment files end with ".aln" |
| [`output/2.iqtree/genetrees/s13.gblock_trimal.py`](output/2.iqtree/genetrees/s13.gblock_trimal.py) | Further filter Gblocks-trimmed multiple sequence alignments with trimAl |
| [`output/2.iqtree/genetrees/s14.run_supertree_iqtree.py`](output/2.iqtree/genetrees/s14.run_supertree_iqtree.py) | Generate and run IQ-TREE commands for the supertree analysis |
| [`output/2.iqtree_cds/genetrees/gCF_sCF/h1.remove_gene_tree.py`](output/2.iqtree_cds/genetrees/gCF_sCF/h1.remove_gene_tree.py) | Strip gene IDs from gene trees, keeping species names |
| [`output/2.iqtree_cds/genetrees/s13.gblock_trimal.py`](output/2.iqtree_cds/genetrees/s13.gblock_trimal.py) | Further filter Gblocks-trimmed multiple sequence alignments with trimAl |
| [`output/2.iqtree_no_Gma/codon12_gCF_sCF/h1.remove_gene_tree.py`](output/2.iqtree_no_Gma/codon12_gCF_sCF/h1.remove_gene_tree.py) | Strip gene IDs from gene trees, keeping species names |
| [`output/2.iqtree_no_Gma/concatence_gCF_sCF/CDS_gcf_scf/h1.remove_gene_tree.py`](output/2.iqtree_no_Gma/concatence_gCF_sCF/CDS_gcf_scf/h1.remove_gene_tree.py) | Strip gene IDs from gene trees, keeping species names |
| [`output/2.iqtree_no_Gma/concatence_gCF_sCF/codon12_gcf_scf/h1.remove_gene_tree.py`](output/2.iqtree_no_Gma/concatence_gCF_sCF/codon12_gcf_scf/h1.remove_gene_tree.py) | Strip gene IDs from gene trees, keeping species names |
| [`output/2.iqtree_no_Gma/get_concatence.py`](output/2.iqtree_no_Gma/get_concatence.py) | Concatenate CDS alignments |
| [`output/2.iqtree_no_Gma/run_iqtree.py`](output/2.iqtree_no_Gma/run_iqtree.py) | Generate and run IQ-TREE commands for the supertree analysis |
| [`output/2.iqtree_no_Gma/run_trimal.py`](output/2.iqtree_no_Gma/run_trimal.py) | Further filter Gblocks-trimmed multiple sequence alignments with trimAl |
| [`output/2.iqtree_no_ZM4_outgroup/s14.run_supertree_iqtree.py`](output/2.iqtree_no_ZM4_outgroup/s14.run_supertree_iqtree.py) | Generate and run IQ-TREE commands for the supertree analysis |
| [`output/3.interpro_anno/h1.get_pep_fa.py`](output/3.interpro_anno/h1.get_pep_fa.py) | Get gene lists per pangenome category and extract their sequences |
| [`output/3.interpro_anno/h2.get_pep_from_id.py`](output/3.interpro_anno/h2.get_pep_from_id.py) | Get gene lists per pangenome category and extract their sequences |
| [`output/5.PAML/cds_paml2/baseml/baseml.ctl`](output/5.PAML/cds_paml2/baseml/baseml.ctl) | PAML control file |
| [`output/5.PAML/cds_paml2/mcmctree2/mcmctree.ctl`](output/5.PAML/cds_paml2/mcmctree2/mcmctree.ctl) | PAML control file |
| [`output/5.PAML/cds_paml2/mcmctree/mcmctree.ctl`](output/5.PAML/cds_paml2/mcmctree/mcmctree.ctl) | PAML control file |
| [`output/5.PAML/cds_paml2/mcmctree/tmp0001.ctl`](output/5.PAML/cds_paml2/mcmctree/tmp0001.ctl) | PAML control file |
| [`output/5.PAML/cds_paml/baseml/baseml.ctl`](output/5.PAML/cds_paml/baseml/baseml.ctl) | PAML control file |
| [`output/5.PAML/cds_paml/data/4.4Dsites.pl`](output/5.PAML/cds_paml/data/4.4Dsites.pl) |  |
| [`output/5.PAML/cds_paml/data/fasta2phy.py`](output/5.PAML/cds_paml/data/fasta2phy.py) | fasta2phy.py |
| [`output/5.PAML/cds_paml/mcmctree2/mcmctree.ctl`](output/5.PAML/cds_paml/mcmctree2/mcmctree.ctl) | PAML control file |
| [`output/5.PAML/cds_paml/mcmctree/mcmctree.ctl`](output/5.PAML/cds_paml/mcmctree/mcmctree.ctl) | PAML control file |
| [`output/5.PAML/cds_paml/mcmctree/tmp0001.ctl`](output/5.PAML/cds_paml/mcmctree/tmp0001.ctl) | PAML control file |
| [`output/6.cafe/run_cafe.sh`](output/6.cafe/run_cafe.sh) | Run cafe |
| [`output/7.RNA_tpm/h1.get_geneTPM.py`](output/7.RNA_tpm/h1.get_geneTPM.py) | Extract transcript expression from StringTie Ballgown sample annotation files |
| [`output/7.RNA_tpm/h2.get_panGP_tpm.py`](output/7.RNA_tpm/h2.get_panGP_tpm.py) | Gene expression per pangenome category for each species |
| [`output/8.nucleotide_diversity/h1.get_genefam_gene.py`](output/8.nucleotide_diversity/h1.get_genefam_gene.py) | Get gene IDs for gene families |
| [`output/8.nucleotide_diversity/h2.cal_fam_pai.py`](output/8.nucleotide_diversity/h2.cal_fam_pai.py) | Compute nucleotide diversity per pangenome family category |
| [`output/8.nucleotide_diversity/nucleotide_diversity.R`](output/8.nucleotide_diversity/nucleotide_diversity.R) | Compute nucleotide diversity for each pangenome family category |
| [`output_M46_new/scripts/audit_sequence_characters.sh`](output_M46_new/scripts/audit_sequence_characters.sh) | Audit sequence characters |
| [`output_M46_new/scripts/postprocess_pan_genome_B.R`](output_M46_new/scripts/postprocess_pan_genome_B.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_M46_new/scripts/postprocess_pan_genome_new.R`](output_M46_new/scripts/postprocess_pan_genome_new.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_M46_new/scripts/postprocess_table_only.R`](output_M46_new/scripts/postprocess_table_only.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_M46_new/scripts/postprocess_table_source.R`](output_M46_new/scripts/postprocess_table_source.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_M46_new/scripts/prepare_clean_input.sh`](output_M46_new/scripts/prepare_clean_input.sh) | Prepare clean input |
| [`output_M46_new/scripts/run_input_qc.sh`](output_M46_new/scripts/run_input_qc.sh) | Run input qc |
| [`output_M46_new/scripts/run_orthofinder.sh`](output_M46_new/scripts/run_orthofinder.sh) | Run orthofinder |
| [`output_pangenome_new/scripts/audit_sequence_characters.sh`](output_pangenome_new/scripts/audit_sequence_characters.sh) | Audit sequence characters |
| [`output_pangenome_new/scripts/build_supplementary_workbook.mjs`](output_pangenome_new/scripts/build_supplementary_workbook.mjs) | Build supplementary workbook |
| [`output_pangenome_new/scripts/finalize_completed_status.sh`](output_pangenome_new/scripts/finalize_completed_status.sh) | Finalize completed status |
| [`output_pangenome_new/scripts/postprocess_pan_genome.R`](output_pangenome_new/scripts/postprocess_pan_genome.R) | Post-process an 18-genome OrthoFinder run and draw manuscript figures. |
| [`output_pangenome_new/scripts/prepare_clean_input.sh`](output_pangenome_new/scripts/prepare_clean_input.sh) | Prepare clean input |
| [`output_pangenome_new/scripts/qa_pan_genome_outputs.R`](output_pangenome_new/scripts/qa_pan_genome_outputs.R) |  |
| [`output_pangenome_new/scripts/run_input_qc.sh`](output_pangenome_new/scripts/run_input_qc.sh) | Run input qc |
| [`output_pangenome_new/scripts/run_orthofinder.sh`](output_pangenome_new/scripts/run_orthofinder.sh) | Run orthofinder |
