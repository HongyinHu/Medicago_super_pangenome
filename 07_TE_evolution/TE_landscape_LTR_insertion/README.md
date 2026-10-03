# TE_landscape_LTR_insertion

TE 景观、共享/特异 TE、完整 LTR 聚类与插入时间

| 脚本 | 说明 |
|---|---|
| [`my_run/a1.get_LTR_fa.py`](my_run/a1.get_LTR_fa.py) | 获取fl_LTR序列 |
| [`my_run/a2.cluster_file_build.py`](my_run/a2.cluster_file_build.py) | 获取文件两两组合 |
| [`my_run/a3.TE_vmatch_cluster.py`](my_run/a3.TE_vmatch_cluster.py) | 批量运行vmatch用于大规模聚类计算 |
| [`my_run/a4.get_percentage_syntenic.py`](my_run/a4.get_percentage_syntenic.py) | 获取每个两两成组结果的共线性比例 |
| [`my_run/a5.format_res.py`](my_run/a5.format_res.py) | 将结果进行特定格式格式化 |
| [`my_run/distance_estimation_pairwise_nucleotide.mao`](my_run/distance_estimation_pairwise_nucleotide.mao) |  |
| [`my_run/h1.TE_lanscape_bedtools.py`](my_run/h1.TE_lanscape_bedtools.py) | 蛋白编码基因前后10kp附近,各个类型TE的分布情况 |
| [`my_run/h2.get_distance_near_gene.py`](my_run/h2.get_distance_near_gene.py) | 获取各类型TE和gene的上下游距离 |
| [`my_run/h3.get_maf_block.py`](my_run/h3.get_maf_block.py) | 获取多个基因组序列比对blocks, 每个block至少共享在两个物种中，最后使用bedtools进行合并处理 |
| [`my_run/h4.covert_stand_bed.py`](my_run/h4.covert_stand_bed.py) | 转换链信息，将负链的位置坐标转到相应的正链位置坐标上 |
| [`my_run/h5.get_shared_TE.py`](my_run/h5.get_shared_TE.py) | 获取shared TE 区域 |
| [`my_run/h6.covert_stand_bed_indir.py`](my_run/h6.covert_stand_bed_indir.py) | 转换链信息，将负链的位置坐标转到相应的正链位置坐标上 |
| [`my_run/h7.bed_sort_merge.py`](my_run/h7.bed_sort_merge.py) | bed文件进行sort和merge操作 |
| [`my_run/h8.bedtools_indersect.py`](my_run/h8.bedtools_indersect.py) | bed文件进行sort和merge操作 |
| [`my_run/h9.get_regions_len_indir.py`](my_run/h9.get_regions_len_indir.py) | 获取区间长度累加之和 |
| [`my_run/h10.get_intersect_bed.py`](my_run/h10.get_intersect_bed.py) | 获取beb文件交集，如果A文件完全被B文件覆盖，则输出A文件坐标 |
| [`my_run/M6CC_distance.mao`](my_run/M6CC_distance.mao) |  |
| [`my_run/M6CC_muscle.mao`](my_run/M6CC_muscle.mao) |  |
| [`my_run/muscle_align_nucleotide.mao`](my_run/muscle_align_nucleotide.mao) |  |
| [`my_run/p1.get_bed_from_cactus_shared.py`](my_run/p1.get_bed_from_cactus_shared.py) | 从cactus多序列比对结果提取物种共享片段的bed文件 |
| [`my_run/p2.get_bed_from_cactus_specific.py`](my_run/p2.get_bed_from_cactus_specific.py) | 从cactus多序列比对结果提取物种物种特异性片段的bed文件 |
| [`my_run/p3.covert_stand_bed.py`](my_run/p3.covert_stand_bed.py) | 转换链信息，将负链的位置坐标转到相应的正链位置坐标上 |
| [`my_run/p4.get_genome_bed.py`](my_run/p4.get_genome_bed.py) | 获取基因组bed文件，为差异就是非共享区域 |
| [`my_run/p5.get_regions_len.py`](my_run/p5.get_regions_len.py) | 获取区间长度累加之和 |
| [`my_run/p6.classic_TE.py`](my_run/p6.classic_TE.py) | 对TE注释文件中的类别进行分类，获得相应的bed文件 |
| [`my_run/p7.get_rate_TE.py`](my_run/p7.get_rate_TE.py) | 获取各类型TE的比例 |
| [`my_run/p8.get_intact_LTR.py`](my_run/p8.get_intact_LTR.py) | 获取全长完整的LTR从genome.fa.mod.LTR.intact.gff3文件中 |
| [`my_run/p9.get_shared_intact_TE_fa.py`](my_run/p9.get_shared_intact_TE_fa.py) | 获取在共享区域的LTR的序列，并集有一段在即可 |
| [`my_run/p10.cal_K_value.py`](my_run/p10.cal_K_value.py) | 使用MEGAcc计算K值(sequence divergence rate between its 5` and 3` LTRs) |
| [`my_run/p11.get_k_Value_from_result.py`](my_run/p11.get_k_Value_from_result.py) | 从结果文件提取K值 |
| [`my_run/p12.get_intact_LTR_type.py`](my_run/p12.get_intact_LTR_type.py) | 获取全长完整的LTR从genome.fa.mod.LTR.intact.gff3文件中,分三种类型 |
| [`my_run/p13.get_shared_intact_TE_fa_type.py`](my_run/p13.get_shared_intact_TE_fa_type.py) | 获取在共享区域的LTR的序列，并集有一段在即可 |
| [`my_run/s1.TE_landscape_around_genes2.py`](my_run/s1.TE_landscape_around_genes2.py) | 蛋白编码基因前后10kp附近，各个类型TE的分布情况 |
| [`my_run/s1.TE_landscape_around_genes.py`](my_run/s1.TE_landscape_around_genes.py) | 蛋白编码基因前后10kp附近，各个类型TE的分布情况 |
| [`my_run/s2.pipline_intact_LTR.py`](my_run/s2.pipline_intact_LTR.py) |  |
| [`my_run/s3.pipline_intact_LTR_specific.py`](my_run/s3.pipline_intact_LTR_specific.py) |  |
| [`my_run/s4.run_M46_2_LTR_time.sh`](my_run/s4.run_M46_2_LTR_time.sh) | Recompute intact-LTR K values / insertion times for the new M46 assembly (M46_2), |
| [`output_LTR_time_M46_2/Mfi_46.intact_LTR/p12.get_intact_LTR_type.M46_2.py`](output_LTR_time_M46_2/Mfi_46.intact_LTR/p12.get_intact_LTR_type.M46_2.py) | 获取全长完整的LTR从genome.fa.mod.LTR.intact.gff3文件中,分三种类型 |
