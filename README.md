# Medicago super-pangenome

Analysis code for the genus-wide *Medicago* super-pangenome study, covering genome assessment, repeat and gene annotation, pangenome and phylogeny, WGD and karyotype evolution, centromeres, TE evolution and structural variation.

## Repository structure

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

## Analysis modules

| Category | Module | Description |
|---|---|---|
| Genome survey & assembly QC | [genome_survey](01_genome_survey_QC/genome_survey) | k-mer genome survey (jellyfish + GenomeScope2) |
|  | [assembly_quality](01_genome_survey_QC/assembly_quality) | Assembly quality assessment (BUSCO, Merqury, short-read mapping) |
| Repeat annotation | [TE_annotation](02_repeat_annotation/TE_annotation) | Per-genome EDTA → panEDTA pan-TE library → soft-masking → LAI → TEsorter |
| Gene annotation | [gene_structure_annotation](03_gene_annotation/gene_structure_annotation) | Gene structure annotation: BRAKER3, GeMoMa, PASA and EVidenceModeler |
|  | [gene_function_annotation](03_gene_annotation/gene_function_annotation) | Functional annotation: DIAMOND (NR/Swiss-Prot/TrEMBL), eggNOG, InterProScan, Pfam, KOG |
| Pangenome & phylogeny | [gene_family_pangenome](04_pangenome_phylogeny/gene_family_pangenome) | OrthoFinder gene families and pangenome categories, species tree (IQ-TREE), divergence times (MCMCTree), CAFE5, Ka/Ks |
| WGD & karyotype evolution | [WGD_Ks](05_WGD_karyotype/WGD_Ks) | Ks distributions of paralogs with WGDI to identify WGD events |
|  | [karyotype_reconstruction](05_WGD_karyotype/karyotype_reconstruction) | Ancestral karyotype reconstruction with WGDI, chromosome breakpoints and karyotype evolution |
| Centromere | [centromere](06_centromere/centromere) | Functional centromeres from CENH3 ChIP-seq, satellite repeats, methylation, cross-species centromere synteny and evolution, centromeric TEs/LTRs |
| TE evolution | [TE_landscape_LTR_insertion](07_TE_evolution/TE_landscape_LTR_insertion) | TE landscape, shared/specific TEs, intact-LTR clustering and insertion times |
|  | [TE_subtype_soloLTR](07_TE_evolution/TE_subtype_soloLTR) | TE subtype statistics, solo-LTR/intact-LTR ratios and insertion-time plots |
| Structural variation | [panSV](08_structural_variation/panSV) | SV calling and pan-SV construction: dual-reference HiFi read-based calling (Snakemake + Jasmine + Sniffles2 genotyping), SVGAP assembly-based calling, integration and QC, SV feature/TE/expression analyses |

## Notes

- Each module directory contains a `README.md` listing its scripts with short descriptions.
- In modules with a `my_run/` directory, most Python scripts print the commands to run rather than executing them, e.g. `python my_run/s1.xxx.py <args> > run.sh && bash run.sh`. Script prefixes: `s` = main pipeline step, `h` = helper, `p` = statistics/plotting.
- Absolute paths have been replaced with placeholders: `path/to/project` (project directory), `path/to/home` (software and conda environments), `path/to/data` and `path/to/software`. Adjust them to your local paths before use.
- Some code comments are in Chinese.
- This repository contains code only; no data or results are included.
