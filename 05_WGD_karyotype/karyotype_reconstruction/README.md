# karyotype_reconstruction

Ancestral karyotype reconstruction with WGDI, chromosome breakpoints and karyotype evolution

| Script | Description |
|---|---|
| [`my_run/h1.select_chr.py`](my_run/h1.select_chr.py) | Remove non-essential chromosome fragments from WGDI input files |
| [`my_run/h2.rename_gff_geneid.py`](my_run/h2.rename_gff_geneid.py) | Rename gene IDs in GFF files |
| [`my_run/h3.filter_block_inf.py`](my_run/h3.filter_block_inf.py) | Filter Mpr_Mpr_block_information.csv by block size |
| [`my_run/s1.get_gff_lens.py`](my_run/s1.get_gff_lens.py) | Prepare WGDI-like inputs: |
| [`my_run/s2.run_blastp.py`](my_run/s2.run_blastp.py) | BLASTP of protein sequences |
| [`my_run/s3.run_wgdi_same_sample.py`](my_run/s3.run_wgdi_same_sample.py) | Write WGDI configuration files and run each WGDI step |
| [`my_run/s4.run_wgdi_ksfigure.py`](my_run/s4.run_wgdi_ksfigure.py) | Write WGDI configuration files and run each WGDI step |
| [`my_run/s5.run_mcscan_cds_pep.py`](my_run/s5.run_mcscan_cds_pep.py) | Plot synteny with MCScan (JCVI) |
| [`my_run/s6.run_mcscan_ortholog.py`](my_run/s6.run_mcscan_ortholog.py) | Compute pairwise synteny between species |
| [`my_run/s7.run_mcscan_plot.py`](my_run/s7.run_mcscan_plot.py) | Compute pairwise synteny between species |
| [`output/genome_410_genome_Msa/aak.dotplot2.conf`](output/genome_410_genome_Msa/aak.dotplot2.conf) | Configuration file |
| [`output/genome_410_genome_Msa/aak.dotplot.conf`](output/genome_410_genome_Msa/aak.dotplot.conf) | Configuration file |
| [`output/genome_410_genome_Msa/ak.conf`](output/genome_410_genome_Msa/ak.conf) | Configuration file |
| [`output/genome_410_genome_Msa/ances_wgdi/total.conf`](output/genome_410_genome_Msa/ances_wgdi/total.conf) | Configuration file |
| [`output/genome_410_genome_Msa/k.conf`](output/genome_410_genome_Msa/k.conf) | Configuration file |
| [`output/genome_410_genome_Msa/km.conf`](output/genome_410_genome_Msa/km.conf) | Configuration file |
| [`output/genome_410_genome_Msa/total.conf`](output/genome_410_genome_Msa/total.conf) | Configuration file |
| [`output/genome_nan_genome_Msa/ances_wgdi/total.conf`](output/genome_nan_genome_Msa/ances_wgdi/total.conf) | Configuration file |
| [`output/genome_qing_genome_Msa/ak.conf`](output/genome_qing_genome_Msa/ak.conf) | Configuration file |
| [`output/genome_qing_genome_Msa/icl.conf`](output/genome_qing_genome_Msa/icl.conf) | Configuration file |
| [`output_ED5a_redo_20260929/02_wgdi/Malbus/c.conf`](output_ED5a_redo_20260929/02_wgdi/Malbus/c.conf) | Configuration file |
| [`output_ED5a_redo_20260929/02_wgdi/Malbus/k.conf`](output_ED5a_redo_20260929/02_wgdi/Malbus/k.conf) | Configuration file |
| [`output_ED5a_redo_20260929/02_wgdi/Malbus/km.conf`](output_ED5a_redo_20260929/02_wgdi/Malbus/km.conf) | Configuration file |
| [`output_ED5a_redo_20260929/scripts/amk7_chk.py`](output_ED5a_redo_20260929/scripts/amk7_chk.py) | Does the chr7 block of Mrut/Marc still contain AMK7 genes 6343-6406 |
| [`output_ED5a_redo_20260929/scripts/amk7_end.py`](output_ED5a_redo_20260929/scripts/amk7_end.py) | Where does the distal end of AMK7 (genes 6200-6466) sit in each genome? |
| [`output_ED5a_redo_20260929/scripts/breakpoint_table.py`](output_ED5a_redo_20260929/scripts/breakpoint_table.py) | Supplementary table: inter-AMK breakpoints of M. polymorpha located on two |
| [`output_ED5a_redo_20260929/scripts/cen5_map.sh`](output_ED5a_redo_20260929/scripts/cen5_map.sh) | Map the T2T M. sativa subsp. caerulea annotation (CENH3 coordinates) onto the AMK |
| [`output_ED5a_redo_20260929/scripts/cen8_check.sh`](output_ED5a_redo_20260929/scripts/cen8_check.sh) | Is the AMK8 breakpoint of M. praecox (AMK8 3659/3661) at the ancestral CEN8? |
| [`output_ED5a_redo_20260929/scripts/cen8_region.py`](output_ED5a_redo_20260929/scripts/cen8_region.py) | Where are the AMK8 genes around CEN8 (AMK8 2350-2600) in M. praecox? |
| [`output_ED5a_redo_20260929/scripts/check_cen5.py`](output_ED5a_redo_20260929/scripts/check_cen5.py) | Do the AMK5 breakpoints of M. praecox Chr5 (and M. polymorpha) fall at CEN5? |
| [`output_ED5a_redo_20260929/scripts/dotplot_evidence.py`](output_ED5a_redo_20260929/scripts/dotplot_evidence.py) | Dot-plot evidence for a karyotype change (layout after the AZAK/TLI/Aal figure). |
| [`output_ED5a_redo_20260929/scripts/dotplot_zoom.py`](output_ED5a_redo_20260929/scripts/dotplot_zoom.py) | Zoomed dot plots for small inter-chromosomal translocations. |
| [`output_ED5a_redo_20260929/scripts/gap_diag.py`](output_ED5a_redo_20260929/scripts/gap_diag.py) | List every internal gap between km_result segments and explain each with the |
| [`output_ED5a_redo_20260929/scripts/junctions.py`](output_ED5a_redo_20260929/scripts/junctions.py) | Inter-AMK junctions (adjacent major segments from different AMK chromosomes) |
| [`output_ED5a_redo_20260929/scripts/karyo_struct.py`](output_ED5a_redo_20260929/scripts/karyo_struct.py) | Per-genome chromosome composition in AMK terms, from the filtered |
| [`output_ED5a_redo_20260929/scripts/karyotype_backbone_pathways_R.R`](output_ED5a_redo_20260929/scripts/karyotype_backbone_pathways_R.R) | R-only redraw of the AMK karyotype backbone and the evidence-ranked |
| [`output_ED5a_redo_20260929/scripts/ks_survey.py`](output_ED5a_redo_20260929/scripts/ks_survey.py) | Survey block-level and pair-level Ks of every species-vs-AMK WGDI run. |
| [`output_ED5a_redo_20260929/scripts/map_t2t_msa.sh`](output_ED5a_redo_20260929/scripts/map_t2t_msa.sh) | Map T2T M. sativa subsp. caerulea genes (CENH3 coordinate system) onto the AMK gene |
| [`output_ED5a_redo_20260929/scripts/missing_arm.py`](output_ED5a_redo_20260929/scripts/missing_arm.py) | Are AMK4 1-2478 (Marc) and AMK2 1-2735 (Mlup) absent from the assembled |
| [`output_ED5a_redo_20260929/scripts/one_to_one.py`](output_ED5a_redo_20260929/scripts/one_to_one.py) | Genome-wide one-to-one filter on WGDI block tables. |
| [`output_ED5a_redo_20260929/scripts/plot_ed5a.py`](output_ED5a_redo_20260929/scripts/plot_ed5a.py) | Redraw ED5a as a radial layout: AMK in the centre, 19 genomes around it. |
| [`output_ED5a_redo_20260929/scripts/plot_ed5a_r.R`](output_ED5a_redo_20260929/scripts/plot_ed5a_r.R) | ED5a radial redraw in base R. All drawing and visual QA stay in R. |
| [`output_ED5a_redo_20260929/scripts/plot_karyo_evolution.py`](output_ED5a_redo_20260929/scripts/plot_karyo_evolution.py) | Karyotype evolution of Medicago from the AMK, drawn on the species tree |
| [`output_ED5a_redo_20260929/scripts/plot_karyo_pathways.py`](output_ED5a_redo_20260929/scripts/plot_karyo_pathways.py) | Panel b: step-by-step inter-chromosomal rearrangement pathways from the AMK. |
| [`output_ED5a_redo_20260929/scripts/plot_karyo_scheme.py`](output_ED5a_redo_20260929/scripts/plot_karyo_scheme.py) | Karyotype-evolution scheme (Brassicaceae ACBK style, not a dated binary tree). |
| [`output_ED5a_redo_20260929/scripts/plot_ks_M46_compare.py`](output_ED5a_redo_20260929/scripts/plot_ks_M46_compare.py) | Intra-genome Ks distribution of M. fischeriana: old (M46_1) vs new (M46_2) assembly. |
| [`output_ED5a_redo_20260929/scripts/prepare_confs.py`](output_ED5a_redo_20260929/scripts/prepare_confs.py) | Stage per-species WGDI runs for the filtered ED5a karyotype mapping. |
| [`output_ED5a_redo_20260929/scripts/run_diamond_Mpol_R108.sh`](output_ED5a_redo_20260929/scripts/run_diamond_Mpol_R108.sh) | M. polymorpha (T2T, genome_nan annotation used in the karyotype runs) vs |
| [`output_ED5a_redo_20260929/scripts/run_diamond_Mpra_R108.sh`](output_ED5a_redo_20260929/scripts/run_diamond_Mpra_R108.sh) | M. praecox proteins vs M. truncatula R108 proteins (same R108 annotation as ED5a, |
| [`output_ED5a_redo_20260929/scripts/run_diamond_x8_R108.sh`](output_ED5a_redo_20260929/scripts/run_diamond_x8_R108.sh) | Proteins of the x = 8 genomes with inter-chromosomal changes vs M. truncatula R108 |
| [`output_ED5a_redo_20260929/scripts/run_species.sh`](output_ED5a_redo_20260929/scripts/run_species.sh) | Run WGDI -c -> one-to-one -> -km -> -k for one staged species directory. |
| [`output_ED5a_redo_20260929/scripts/run_wgd_M46_2.sh`](output_ED5a_redo_20260929/scripts/run_wgd_M46_2.sh) | Intra-genome Ks distribution of M. fischeriana with the NEW assembly (genome_M46_2), |
| [`output_ED5a_redo_20260929/scripts/small_insert.py`](output_ED5a_redo_20260929/scripts/small_insert.py) | Find short foreign-AMK segments embedded inside a chromosome (flanked on both |
| [`output_ED5a_redo_20260929/scripts/summarize.py`](output_ED5a_redo_20260929/scripts/summarize.py) | Per-species filtering / mapping summary for the ED5a redo. |
| [`output_finally/ED5_panel_a_source_20260727/extract_ed5_panel_a_from_wgdi.py`](output_finally/ED5_panel_a_source_20260727/extract_ed5_panel_a_from_wgdi.py) | Extract ED5 panel-a source data from the underlying WGDI analysis outputs. |
| [`output_finally/Mpo_diagnosis_20260618/plot_Mpo_diagnosis.R`](output_finally/Mpo_diagnosis_20260618/plot_Mpo_diagnosis.R) | Plot Mpo diagnosis |
| [`output_finally/prepare_genome_nan_wgdi_20260618.py`](output_finally/prepare_genome_nan_wgdi_20260618.py) | Prepare genome nan wgdi 20260618 |
| [`output_Kary/genome_395/total.conf`](output_Kary/genome_395/total.conf) | Configuration file |
| [`output_Kary/genome_M46_2/prepare_inputs.sh`](output_Kary/genome_M46_2/prepare_inputs.sh) | Prepare inputs |
| [`output_Kary/genome_M46_2/run_analysis.sh`](output_Kary/genome_M46_2/run_analysis.sh) | Run analysis |
| [`output_Kary/genome_M46_2/run_wgdi.sh`](output_Kary/genome_M46_2/run_wgdi.sh) | Run wgdi |
| [`output_Kary/genome_Mpo_rerun_20260618/total_karyotype_png.conf`](output_Kary/genome_Mpo_rerun_20260618/total_karyotype_png.conf) | Configuration file |
