# TE_landscape_LTR_insertion

TE landscape, shared/specific TEs, intact-LTR clustering and insertion times

| Script | Description |
|---|---|
| [`my_run/a1.get_LTR_fa.py`](my_run/a1.get_LTR_fa.py) | Extract full-length LTR sequences |
| [`my_run/a2.cluster_file_build.py`](my_run/a2.cluster_file_build.py) | Generate pairwise file combinations |
| [`my_run/a3.TE_vmatch_cluster.py`](my_run/a3.TE_vmatch_cluster.py) | Batch vmatch runs for large-scale LTR clustering |
| [`my_run/a4.get_percentage_syntenic.py`](my_run/a4.get_percentage_syntenic.py) | Compute syntenic proportion for each pairwise comparison |
| [`my_run/a5.format_res.py`](my_run/a5.format_res.py) | Reformat results |
| [`my_run/distance_estimation_pairwise_nucleotide.mao`](my_run/distance_estimation_pairwise_nucleotide.mao) | MEGA analysis options |
| [`my_run/h1.TE_lanscape_bedtools.py`](my_run/h1.TE_lanscape_bedtools.py) | Distribution of TE classes within 10 kb of protein-coding genes |
| [`my_run/h2.get_distance_near_gene.py`](my_run/h2.get_distance_near_gene.py) | Distance between each TE class and flanking genes |
| [`my_run/h3.get_maf_block.py`](my_run/h3.get_maf_block.py) | Extract alignment blocks shared by at least two species and merge with bedtools |
| [`my_run/h4.covert_stand_bed.py`](my_run/h4.covert_stand_bed.py) | Convert minus-strand coordinates to plus-strand coordinates |
| [`my_run/h5.get_shared_TE.py`](my_run/h5.get_shared_TE.py) | Extract shared TE regions |
| [`my_run/h6.covert_stand_bed_indir.py`](my_run/h6.covert_stand_bed_indir.py) | Convert minus-strand coordinates to plus-strand coordinates |
| [`my_run/h7.bed_sort_merge.py`](my_run/h7.bed_sort_merge.py) | Sort and merge BED files |
| [`my_run/h8.bedtools_indersect.py`](my_run/h8.bedtools_indersect.py) | Sort and merge BED files |
| [`my_run/h9.get_regions_len_indir.py`](my_run/h9.get_regions_len_indir.py) | Sum interval lengths |
| [`my_run/h10.get_intersect_bed.py`](my_run/h10.get_intersect_bed.py) | Intersect BED files, reporting A features fully covered by B |
| [`my_run/M6CC_distance.mao`](my_run/M6CC_distance.mao) | MEGA analysis options |
| [`my_run/M6CC_muscle.mao`](my_run/M6CC_muscle.mao) | MEGA analysis options |
| [`my_run/muscle_align_nucleotide.mao`](my_run/muscle_align_nucleotide.mao) | MEGA analysis options |
| [`my_run/p1.get_bed_from_cactus_shared.py`](my_run/p1.get_bed_from_cactus_shared.py) | Extract BED of regions shared among species from the Cactus alignment |
| [`my_run/p2.get_bed_from_cactus_specific.py`](my_run/p2.get_bed_from_cactus_specific.py) | Extract BED of species-specific regions from the Cactus alignment |
| [`my_run/p3.covert_stand_bed.py`](my_run/p3.covert_stand_bed.py) | Convert minus-strand coordinates to plus-strand coordinates |
| [`my_run/p4.get_genome_bed.py`](my_run/p4.get_genome_bed.py) | Build genome BED to derive non-shared regions |
| [`my_run/p5.get_regions_len.py`](my_run/p5.get_regions_len.py) | Sum interval lengths |
| [`my_run/p6.classic_TE.py`](my_run/p6.classic_TE.py) | Split TE annotation by class into BED files |
| [`my_run/p7.get_rate_TE.py`](my_run/p7.get_rate_TE.py) | Proportion of each TE class |
| [`my_run/p8.get_intact_LTR.py`](my_run/p8.get_intact_LTR.py) | Extract intact full-length LTRs from genome.fa.mod.LTR.intact.gff3 |
| [`my_run/p9.get_shared_intact_TE_fa.py`](my_run/p9.get_shared_intact_TE_fa.py) | Extract LTR sequences overlapping shared regions |
| [`my_run/p10.cal_K_value.py`](my_run/p10.cal_K_value.py) | Compute K (divergence between 5' and 3' LTRs) with MEGA-CC |
| [`my_run/p11.get_k_Value_from_result.py`](my_run/p11.get_k_Value_from_result.py) | Extract K values (LTR divergence) from result files |
| [`my_run/p12.get_intact_LTR_type.py`](my_run/p12.get_intact_LTR_type.py) | Extract intact full-length LTRs from genome.fa.mod.LTR.intact.gff3 by type (Copia/Gypsy/unknown) |
| [`my_run/p13.get_shared_intact_TE_fa_type.py`](my_run/p13.get_shared_intact_TE_fa_type.py) | Extract LTR sequences overlapping shared regions |
| [`my_run/s1.TE_landscape_around_genes2.py`](my_run/s1.TE_landscape_around_genes2.py) | Distribution of TE classes within 10 kb of protein-coding genes |
| [`my_run/s1.TE_landscape_around_genes.py`](my_run/s1.TE_landscape_around_genes.py) | Distribution of TE classes within 10 kb of protein-coding genes |
| [`my_run/s2.pipline_intact_LTR.py`](my_run/s2.pipline_intact_LTR.py) |  |
| [`my_run/s3.pipline_intact_LTR_specific.py`](my_run/s3.pipline_intact_LTR_specific.py) |  |
| [`my_run/s4.run_M46_2_LTR_time.sh`](my_run/s4.run_M46_2_LTR_time.sh) | Recompute intact-LTR K values / insertion times for the new M46 assembly (M46_2), |
| [`output_LTR_time_M46_2/Mfi_46.intact_LTR/p12.get_intact_LTR_type.M46_2.py`](output_LTR_time_M46_2/Mfi_46.intact_LTR/p12.get_intact_LTR_type.M46_2.py) | Extract intact full-length LTRs from genome.fa.mod.LTR.intact.gff3 by type (Copia/Gypsy/unknown) |
