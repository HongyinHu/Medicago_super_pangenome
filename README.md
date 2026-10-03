# Medicago super-pangenome

苜蓿属（*Medicago*）超级泛基因组研究的分析代码，包括基因组评估、重复序列与基因注释、泛基因组与系统发育、全基因组复制与核型演化、着丝粒、转座子演化以及结构变异与泛 SV。

## 目录

```text
01_genome_survey_QC/
    genome_survey/
    assembly_quality/
02_repeat_annotation/
    TE_annotation/
03_gene_annotation/
    gene_structure_annotation/
    gene_function_annotation/
04_pangenome_phylogeny/
    gene_family_pangenome/
05_WGD_karyotype/
    WGD_Ks/
    karyotype_reconstruction/
06_centromere/
    centromere/
07_TE_evolution/
    TE_landscape_LTR_insertion/
    TE_subtype_soloLTR/
08_structural_variation/
    panSV/
```

## 分析模块

| 类别 | 模块 | 内容 |
|---|---|---|
| 基因组 survey 与组装质量评估 | [genome_survey](01_genome_survey_QC/genome_survey) | k-mer 基因组 survey（jellyfish + GenomeScope2） |
|  | [assembly_quality](01_genome_survey_QC/assembly_quality) | 组装质量评估（BUSCO、Merqury、reads 回比） |
| 重复序列注释 | [TE_annotation](02_repeat_annotation/TE_annotation) | 单基因组 EDTA → panEDTA 泛 TE 库 → 软屏蔽 → LAI → TEsorter |
| 基因结构与功能注释 | [gene_structure_annotation](03_gene_annotation/gene_structure_annotation) | 基因结构注释：BRAKER3、GeMoMa、PASA、EVidenceModeler 整合 |
|  | [gene_function_annotation](03_gene_annotation/gene_function_annotation) | 基因功能注释：DIAMOND（NR/SwissProt/TrEMBL）、eggNOG、InterProScan、Pfam、KOG |
| 泛基因组与系统发育 | [gene_family_pangenome](04_pangenome_phylogeny/gene_family_pangenome) | OrthoFinder 基因家族与泛基因组分类、物种树（IQ-TREE）、分化时间（MCMCTree）、CAFE5、Ka/Ks |
| 全基因组复制与核型演化 | [WGD_Ks](05_WGD_karyotype/WGD_Ks) | WGDI 同源基因 Ks 分布与 WGD 识别 |
|  | [karyotype_reconstruction](05_WGD_karyotype/karyotype_reconstruction) | WGDI 祖先核型重建、染色体断点与核型演化路径 |
| 着丝粒分析 | [centromere](06_centromere/centromere) | CENH3 ChIP-seq 定义功能着丝粒、卫星重复、甲基化、跨物种着丝粒共线性与演化、着丝粒 TE/LTR |
| 转座子演化 | [TE_landscape_LTR_insertion](07_TE_evolution/TE_landscape_LTR_insertion) | TE 景观、共享/特异 TE、完整 LTR 聚类与插入时间 |
|  | [TE_subtype_soloLTR](07_TE_evolution/TE_subtype_soloLTR) | TE 亚型统计、solo-LTR/完整 LTR 比例与插入时间作图 |
| 结构变异与泛 SV | [panSV](08_structural_variation/panSV) | SV 鉴定与泛 SV 构建：双参考 HiFi reads 多软件鉴定（Snakemake + Jasmine + Sniffles2 回填基因型）、SVGAP 组装法、整合与 QC、SV 特征/TE/表达分析 |

## 说明

- 每个模块目录下的 `README.md` 列出了该模块的脚本及简要说明。
- 早期模块（含 `my_run/` 目录）中的 Python 脚本大多只打印要执行的命令，用法为 `python my_run/s1.xxx.py <参数> > run.sh && bash run.sh`；脚本前缀 `s` 为主流程步骤，`h` 为辅助脚本，`p` 为统计与作图。
- 脚本中的路径已替换为占位符 `path/to/project`（项目目录）和 `path/to/home`（软件与 conda 环境所在目录），使用前请改为本地路径。
- 仓库只包含代码，不含数据和结果。
