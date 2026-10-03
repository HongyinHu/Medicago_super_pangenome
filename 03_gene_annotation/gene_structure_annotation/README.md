# gene_structure_annotation

基因结构注释：BRAKER3、GeMoMa、PASA、EVidenceModeler 整合

| 脚本 | 说明 |
|---|---|
| [`00_inputs/my_species_RNA/genome_Msa1/illumina-seq/download_jing_ye_illumina_RNA.sh`](00_inputs/my_species_RNA/genome_Msa1/illumina-seq/download_jing_ye_illumina_RNA.sh) | 下载：jing ye illumina RNA（据文件名） |
| [`00_logs_recovery/workflow_patch_staging_20260710/build_ready_evm_inputs.sh`](00_logs_recovery/workflow_patch_staging_20260710/build_ready_evm_inputs.sh) | 构建：ready evm inputs（据文件名） |
| [`00_logs_recovery/workflow_patch_staging_20260710/prepare_evm_sample.sh`](00_logs_recovery/workflow_patch_staging_20260710/prepare_evm_sample.sh) | 准备输入：evm sample（据文件名） |
| [`00_logs_recovery/workflow_patch_staging_20260710/run_evm_sample.sh`](00_logs_recovery/workflow_patch_staging_20260710/run_evm_sample.sh) | 运行：evm sample（据文件名） |
| [`00_logs_recovery/workflow_patch_staging_20260710/run_gemoma_strict_v2.sh`](00_logs_recovery/workflow_patch_staging_20260710/run_gemoma_strict_v2.sh) | 运行：gemoma strict v2（据文件名） |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/cds_with_upstream_support.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/cds_with_upstream_support.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/combine_gff_records.pl`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/combine_gff_records.pl) | Alex Lomsadze, Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/count_cds_overlaps.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/count_cds_overlaps.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/flag_top_proteins.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/flag_top_proteins.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/gff_from_region_to_contig.pl`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/gff_from_region_to_contig.pl) | Alex Lomsadze |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/make_chains.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/make_chains.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/nucseq_for_selected_genes.pl`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/nucseq_for_selected_genes.pl) | Alex Lomsadze |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/print_high_confidence.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/print_high_confidence.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/print_longest_isoform.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/print_longest_isoform.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/proteins_from_gtf.pl`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/proteins_from_gtf.pl) | Alex Lomsadze |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/prothint2augustus.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/prothint2augustus.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/prothint.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/prothint.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/run_spliced_alignment.pl`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/run_spliced_alignment.pl) | Tomas Bruna, Alex Lomsadze |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/run_spliced_alignment_pbs.pl`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/run_spliced_alignment_pbs.pl) | Tomas Bruna, Alex Lomsadze |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/select_best_proteins.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/select_best_proteins.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/select_for_next_iteration.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/select_for_next_iteration.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/spaln_to_gff.py`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/spaln_to_gff.py) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/spalnBatch.sh`](02_abinitio_prediction/BRAKER3_protein/run/prothint_patch_skip_malformed_spaln_20260627_064103/bin/spalnBatch.sh) | Tomas Bruna |
| [`02_abinitio_prediction/BRAKER3_protein/run/run_braker3_protein_only_worker.sh`](02_abinitio_prediction/BRAKER3_protein/run/run_braker3_protein_only_worker.sh) | 运行：braker3 protein only worker（据文件名） |
| [`02_abinitio_prediction/BRAKER3_protein/run/run_braker3_protein_worker.sh`](02_abinitio_prediction/BRAKER3_protein/run/run_braker3_protein_worker.sh) | 运行：braker3 protein worker（据文件名） |
| [`02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_one_genome_Msa1_20260705.sh`](02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_one_genome_Msa1_20260705.sh) | 运行：braker3 rna protein one genome Msa1 20260705（据文件名） |
| [`02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_safe_worker.sh`](02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_safe_worker.sh) | 运行：braker3 rna protein safe worker（据文件名） |
| [`02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_worker.sh`](02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_worker.sh) | 运行：braker3 rna protein worker（据文件名） |
| [`02_abinitio_prediction/run/build_related_protein_for_braker3.sh`](02_abinitio_prediction/run/build_related_protein_for_braker3.sh) | 构建：related protein for braker3（据文件名） |
| [`02_abinitio_prediction/run/run_braker3_pasa_worker.sh`](02_abinitio_prediction/run/run_braker3_pasa_worker.sh) | 运行：braker3 pasa worker（据文件名） |
| [`03_homology_prediction/GeMoMa_strict_v2/run/run_gemoma_strict_v2.sh`](03_homology_prediction/GeMoMa_strict_v2/run/run_gemoma_strict_v2.sh) | 运行：gemoma strict v2（据文件名） |
| [`04_transcript_prediction/PASA/genome_395/genome_395.alignAssembly.config`](04_transcript_prediction/PASA/genome_395/genome_395.alignAssembly.config) | templated variables to be replaced exist as <__var_name__> |
| [`04_transcript_prediction/PASA/run/run_pasa_one_sample_20260704.sh`](04_transcript_prediction/PASA/run/run_pasa_one_sample_20260704.sh) | 运行：pasa one sample 20260704（据文件名） |
| [`04_transcript_prediction/PASA/run/run_pasa_worker.sh`](04_transcript_prediction/PASA/run/run_pasa_worker.sh) | 运行：pasa worker（据文件名） |
| [`04_transcript_prediction/PASA_inputs/pacbio_bam_processing/run/build_pacbio_bam_manifest.sh`](04_transcript_prediction/PASA_inputs/pacbio_bam_processing/run/build_pacbio_bam_manifest.sh) | 构建：pacbio bam manifest（据文件名） |
| [`04_transcript_prediction/PASA_inputs/pacbio_bam_processing/run/run_pacbio_isoseq_worker.sh`](04_transcript_prediction/PASA_inputs/pacbio_bam_processing/run/run_pacbio_isoseq_worker.sh) | 运行：pacbio isoseq worker（据文件名） |
| [`04_transcript_prediction/PASA_inputs/run/collapse_isoforms_by_sam_no_rep.py`](04_transcript_prediction/PASA_inputs/run/collapse_isoforms_by_sam_no_rep.py) |  |
| [`04_transcript_prediction/PASA_inputs/run/collect_all_pasa_transcripts.sh`](04_transcript_prediction/PASA_inputs/run/collect_all_pasa_transcripts.sh) | 收集结果：all pasa transcripts（据文件名） |
| [`04_transcript_prediction/PASA_inputs/run/collect_pasa_transcripts.sh`](04_transcript_prediction/PASA_inputs/run/collect_pasa_transcripts.sh) | 收集结果：pasa transcripts（据文件名） |
| [`04_transcript_prediction/PASA_inputs/run/prepare_existing_isoseq_fasta.sh`](04_transcript_prediction/PASA_inputs/run/prepare_existing_isoseq_fasta.sh) | 准备输入：existing isoseq fasta（据文件名） |
| [`04_transcript_prediction/PASA_inputs/run/run_collect_all_pasa_transcripts_locked_20260711.sh`](04_transcript_prediction/PASA_inputs/run/run_collect_all_pasa_transcripts_locked_20260711.sh) | 运行：collect all pasa transcripts locked 20260711（据文件名） |
| [`04_transcript_prediction/PASA_inputs/run/run_stringtie_for_pasa_worker.sh`](04_transcript_prediction/PASA_inputs/run/run_stringtie_for_pasa_worker.sh) | 运行：stringtie for pasa worker（据文件名） |
| [`04_transcript_prediction/run/build_rnaseq_pairs.py`](04_transcript_prediction/run/build_rnaseq_pairs.py) | 构建：rnaseq pairs（据文件名） |
| [`04_transcript_prediction/run/run_hisat2_worker.sh`](04_transcript_prediction/run/run_hisat2_worker.sh) | 运行：hisat2 worker（据文件名） |
| [`05_EVM_integration/run/build_ready_evm_inputs.sh`](05_EVM_integration/run/build_ready_evm_inputs.sh) | 构建：ready evm inputs（据文件名） |
| [`05_EVM_integration/run/filter_evm_models_v2.py`](05_EVM_integration/run/filter_evm_models_v2.py) | 过滤：evm models v2（据文件名） |
| [`05_EVM_integration/run/prepare_evm_sample.sh`](05_EVM_integration/run/prepare_evm_sample.sh) | 准备输入：evm sample（据文件名） |
| [`05_EVM_integration/run/run_evm_sample.sh`](05_EVM_integration/run/run_evm_sample.sh) | 运行：evm sample（据文件名） |
| [`05_EVM_integration/run/run_evm_worker.sh`](05_EVM_integration/run/run_evm_worker.sh) | 运行：evm worker（据文件名） |
| [`scripts/collect_current_evm_stats.py`](scripts/collect_current_evm_stats.py) | 收集结果：current evm stats（据文件名） |
| [`scripts/link_my_species_RNA.sh`](scripts/link_my_species_RNA.sh) | 建立软链接：my species RNA（据文件名） |
| [`scripts/link_project_T2T_RNA_extra.sh`](scripts/link_project_T2T_RNA_extra.sh) | 建立软链接：project T2T RNA extra（据文件名） |
| [`scripts/verify_current_evm_stats.sh`](scripts/verify_current_evm_stats.sh) | 命令：`base=N_1.coding_gene_anno` |
