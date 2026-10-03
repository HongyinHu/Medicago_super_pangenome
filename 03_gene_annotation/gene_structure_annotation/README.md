# gene_structure_annotation

Gene structure annotation: BRAKER3, GeMoMa, PASA and EVidenceModeler

| Script | Description |
|---|---|
| [`00_inputs/my_species_RNA/genome_Msa1/illumina-seq/download_jing_ye_illumina_RNA.sh`](00_inputs/my_species_RNA/genome_Msa1/illumina-seq/download_jing_ye_illumina_RNA.sh) | Download jing ye illumina RNA |
| [`00_logs_recovery/workflow_patch_staging_20260710/build_ready_evm_inputs.sh`](00_logs_recovery/workflow_patch_staging_20260710/build_ready_evm_inputs.sh) | Build ready evm inputs |
| [`00_logs_recovery/workflow_patch_staging_20260710/prepare_evm_sample.sh`](00_logs_recovery/workflow_patch_staging_20260710/prepare_evm_sample.sh) | Prepare evm sample |
| [`00_logs_recovery/workflow_patch_staging_20260710/run_evm_sample.sh`](00_logs_recovery/workflow_patch_staging_20260710/run_evm_sample.sh) | Run evm sample |
| [`00_logs_recovery/workflow_patch_staging_20260710/run_gemoma_strict_v2.sh`](00_logs_recovery/workflow_patch_staging_20260710/run_gemoma_strict_v2.sh) | Run gemoma strict v2 |
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
| [`02_abinitio_prediction/BRAKER3_protein/run/run_braker3_protein_only_worker.sh`](02_abinitio_prediction/BRAKER3_protein/run/run_braker3_protein_only_worker.sh) | Run braker3 protein only worker |
| [`02_abinitio_prediction/BRAKER3_protein/run/run_braker3_protein_worker.sh`](02_abinitio_prediction/BRAKER3_protein/run/run_braker3_protein_worker.sh) | Run braker3 protein worker |
| [`02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_one_genome_Msa1_20260705.sh`](02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_one_genome_Msa1_20260705.sh) | Run braker3 rna protein one genome Msa1 20260705 |
| [`02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_safe_worker.sh`](02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_safe_worker.sh) | Run braker3 rna protein safe worker |
| [`02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_worker.sh`](02_abinitio_prediction/BRAKER3_RNA_protein/run/run_braker3_rna_protein_worker.sh) | Run braker3 rna protein worker |
| [`02_abinitio_prediction/run/build_related_protein_for_braker3.sh`](02_abinitio_prediction/run/build_related_protein_for_braker3.sh) | Build related protein for braker3 |
| [`02_abinitio_prediction/run/run_braker3_pasa_worker.sh`](02_abinitio_prediction/run/run_braker3_pasa_worker.sh) | Run braker3 pasa worker |
| [`03_homology_prediction/GeMoMa_strict_v2/run/run_gemoma_strict_v2.sh`](03_homology_prediction/GeMoMa_strict_v2/run/run_gemoma_strict_v2.sh) | Run gemoma strict v2 |
| [`04_transcript_prediction/PASA/genome_395/genome_395.alignAssembly.config`](04_transcript_prediction/PASA/genome_395/genome_395.alignAssembly.config) | templated variables to be replaced exist as <__var_name__> |
| [`04_transcript_prediction/PASA/run/run_pasa_one_sample_20260704.sh`](04_transcript_prediction/PASA/run/run_pasa_one_sample_20260704.sh) | Run pasa one sample 20260704 |
| [`04_transcript_prediction/PASA/run/run_pasa_worker.sh`](04_transcript_prediction/PASA/run/run_pasa_worker.sh) | Run pasa worker |
| [`04_transcript_prediction/PASA_inputs/pacbio_bam_processing/run/build_pacbio_bam_manifest.sh`](04_transcript_prediction/PASA_inputs/pacbio_bam_processing/run/build_pacbio_bam_manifest.sh) | Build pacbio bam manifest |
| [`04_transcript_prediction/PASA_inputs/pacbio_bam_processing/run/run_pacbio_isoseq_worker.sh`](04_transcript_prediction/PASA_inputs/pacbio_bam_processing/run/run_pacbio_isoseq_worker.sh) | Run pacbio isoseq worker |
| [`04_transcript_prediction/PASA_inputs/run/collapse_isoforms_by_sam_no_rep.py`](04_transcript_prediction/PASA_inputs/run/collapse_isoforms_by_sam_no_rep.py) |  |
| [`04_transcript_prediction/PASA_inputs/run/collect_all_pasa_transcripts.sh`](04_transcript_prediction/PASA_inputs/run/collect_all_pasa_transcripts.sh) | Collect all pasa transcripts |
| [`04_transcript_prediction/PASA_inputs/run/collect_pasa_transcripts.sh`](04_transcript_prediction/PASA_inputs/run/collect_pasa_transcripts.sh) | Collect pasa transcripts |
| [`04_transcript_prediction/PASA_inputs/run/prepare_existing_isoseq_fasta.sh`](04_transcript_prediction/PASA_inputs/run/prepare_existing_isoseq_fasta.sh) | Prepare existing isoseq fasta |
| [`04_transcript_prediction/PASA_inputs/run/run_collect_all_pasa_transcripts_locked_20260711.sh`](04_transcript_prediction/PASA_inputs/run/run_collect_all_pasa_transcripts_locked_20260711.sh) | Run collect all pasa transcripts locked 20260711 |
| [`04_transcript_prediction/PASA_inputs/run/run_stringtie_for_pasa_worker.sh`](04_transcript_prediction/PASA_inputs/run/run_stringtie_for_pasa_worker.sh) | Run stringtie for pasa worker |
| [`04_transcript_prediction/run/build_rnaseq_pairs.py`](04_transcript_prediction/run/build_rnaseq_pairs.py) | Build rnaseq pairs |
| [`04_transcript_prediction/run/run_hisat2_worker.sh`](04_transcript_prediction/run/run_hisat2_worker.sh) | Run hisat2 worker |
| [`05_EVM_integration/run/build_ready_evm_inputs.sh`](05_EVM_integration/run/build_ready_evm_inputs.sh) | Build ready evm inputs |
| [`05_EVM_integration/run/filter_evm_models_v2.py`](05_EVM_integration/run/filter_evm_models_v2.py) | Filter evm models v2 |
| [`05_EVM_integration/run/prepare_evm_sample.sh`](05_EVM_integration/run/prepare_evm_sample.sh) | Prepare evm sample |
| [`05_EVM_integration/run/run_evm_sample.sh`](05_EVM_integration/run/run_evm_sample.sh) | Run evm sample |
| [`05_EVM_integration/run/run_evm_worker.sh`](05_EVM_integration/run/run_evm_worker.sh) | Run evm worker |
| [`scripts/collect_current_evm_stats.py`](scripts/collect_current_evm_stats.py) | Collect current evm stats |
| [`scripts/link_my_species_RNA.sh`](scripts/link_my_species_RNA.sh) | Link my species RNA |
| [`scripts/link_project_T2T_RNA_extra.sh`](scripts/link_project_T2T_RNA_extra.sh) | Link project T2T RNA extra |
| [`scripts/verify_current_evm_stats.sh`](scripts/verify_current_evm_stats.sh) | Command: `base=N_1.coding_gene_anno` |
