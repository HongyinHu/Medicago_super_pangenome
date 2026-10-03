# Medicago super-pangenome

苜蓿属（*Medicago*）超级泛基因组研究的分析代码，包括基因组评估与注释、泛基因组与系统发育、基因组演化（WGD、核型、着丝粒、TE）、结构变异与泛 SV，以及群体重测序与荚果刺、荚果螺旋性状分析。

## 目录

```text
01_genome_assembly_annotation/
    genome_survey/
    assembly_quality/
    TE_annotation/
    gene_structure_annotation/
    gene_function_annotation/
02_pangenome_phylogeny/
    gene_family_pangenome/
03_genome_evolution/
    WGD_Ks/
    karyotype_reconstruction/
    centromere/
    TE_landscape_LTR_insertion/
    TE_subtype_soloLTR/
04_structural_variation/
    panSV/
05_population_traits/
    resequencing_population/
    pod_spine/
    pod_spine_Chr23997_SPL/
    pod_coil/
```

## 分析模块

| 类别 | 模块 | 内容 |
|---|---|---|
| 基因组评估与注释 | [genome_survey](01_genome_assembly_annotation/genome_survey) | k-mer 基因组 survey（jellyfish + GenomeScope2） |
|  | [assembly_quality](01_genome_assembly_annotation/assembly_quality) | 组装质量评估（BUSCO、Merqury、reads 回比） |
|  | [TE_annotation](01_genome_assembly_annotation/TE_annotation) | 重复序列注释：单基因组 EDTA → panEDTA → 软屏蔽 → LAI → TEsorter |
|  | [gene_structure_annotation](01_genome_assembly_annotation/gene_structure_annotation) | 基因结构注释：BRAKER3、GeMoMa、PASA、EVidenceModeler 整合 |
|  | [gene_function_annotation](01_genome_assembly_annotation/gene_function_annotation) | 基因功能注释：DIAMOND（NR/SwissProt/TrEMBL）、eggNOG、InterProScan、Pfam、KOG |
| 泛基因组与系统发育 | [gene_family_pangenome](02_pangenome_phylogeny/gene_family_pangenome) | OrthoFinder 基因家族与泛基因组分类、物种树（IQ-TREE）、分化时间（MCMCTree）、CAFE5、Ka/Ks |
| 基因组演化：WGD、核型、着丝粒与 TE | [WGD_Ks](03_genome_evolution/WGD_Ks) | WGDI 同源基因 Ks 分布与 WGD 识别 |
|  | [karyotype_reconstruction](03_genome_evolution/karyotype_reconstruction) | WGDI 祖先核型重建、染色体断点与核型演化路径 |
|  | [centromere](03_genome_evolution/centromere) | 着丝粒分析：CENH3 ChIP-seq 定义功能着丝粒、卫星重复、甲基化、跨物种着丝粒共线性与演化、着丝粒 TE/LTR |
|  | [TE_landscape_LTR_insertion](03_genome_evolution/TE_landscape_LTR_insertion) | TE 景观、共享/特异 TE、完整 LTR 聚类与插入时间 |
|  | [TE_subtype_soloLTR](03_genome_evolution/TE_subtype_soloLTR) | TE 亚型统计、solo-LTR/完整 LTR 比例与插入时间作图 |
| 结构变异与泛 SV | [panSV](04_structural_variation/panSV) | SV 鉴定与泛 SV 构建：双参考 HiFi reads 多软件鉴定（Snakemake + Jasmine + Sniffles2 回填基因型）、SVGAP 组装法、整合与 QC、SV 特征/TE/表达分析 |
| 群体重测序与性状分析 | [resequencing_population](05_population_traits/resequencing_population) | 重测序 SNP/InDel 鉴定（GATK）与群体结构（PCA、ADMIXTURE、进化树、π/Fst） |
|  | [pod_spine](05_population_traits/pod_spine) | 荚果刺：SV-表型共分离、短读长 SV、SNP/InDel/SV GWAS、候选基因 Chr23997 验证 |
|  | [pod_spine_Chr23997_SPL](05_population_traits/pod_spine_Chr23997_SPL) | 荚果刺候选基因 Chr23997 的同源建树与微共线性（SPL 家族归属） |
|  | [pod_coil](05_population_traits/pod_coil) | 荚果螺旋：SHP 同源与共分离、SNP/InDel/SV GWAS 及模型敏感性 |

## 说明

- 每个模块目录下的 `README.md` 列出了该模块的脚本及简要说明。
- 早期模块（含 `my_run/` 目录）中的 Python 脚本大多只打印要执行的命令，用法为 `python my_run/s1.xxx.py <参数> > run.sh && bash run.sh`；脚本前缀 `s` 为主流程步骤，`h` 为辅助脚本，`p` 为统计与作图。
- 脚本中的路径已替换为占位符 `path/to/project`（项目目录）和 `path/to/home`（软件与 conda 环境所在目录），使用前请改为本地路径。
- 仓库只包含代码，不含数据和结果。
