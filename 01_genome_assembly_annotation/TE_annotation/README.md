# TE_annotation

重复序列注释：单基因组 EDTA → panEDTA → 软屏蔽 → LAI → TEsorter

| 脚本 | 说明 |
|---|---|
| [`02_EDTA_single/LAI_recalc_single_EDTA_20260625/run/collect_LAI_single_summary.py`](02_EDTA_single/LAI_recalc_single_EDTA_20260625/run/collect_LAI_single_summary.py) | 收集结果：LAI single summary（据文件名） |
| [`02_EDTA_single/LAI_recalc_single_EDTA_20260625/run/prepare_manifest.py`](02_EDTA_single/LAI_recalc_single_EDTA_20260625/run/prepare_manifest.py) | 准备输入：manifest（据文件名） |
| [`02_EDTA_single/LAI_recalc_single_EDTA_20260625/run/run_LAI_single_worker.sh`](02_EDTA_single/LAI_recalc_single_EDTA_20260625/run/run_LAI_single_worker.sh) | 运行：LAI single worker（据文件名） |
| [`02_EDTA_single/LAI_recalc_single_EDTA_IDfixed_20260625/run/collect_LAI_single_idfixed_summary.py`](02_EDTA_single/LAI_recalc_single_EDTA_IDfixed_20260625/run/collect_LAI_single_idfixed_summary.py) | 收集结果：LAI single idfixed summary（据文件名） |
| [`02_EDTA_single/LAI_recalc_single_EDTA_IDfixed_20260625/run/prepare_single_idfixed_inputs.py`](02_EDTA_single/LAI_recalc_single_EDTA_IDfixed_20260625/run/prepare_single_idfixed_inputs.py) | 准备输入：single idfixed inputs（据文件名） |
| [`02_EDTA_single/LAI_recalc_single_EDTA_IDfixed_20260625/run/run_LAI_single_idfixed_worker.sh`](02_EDTA_single/LAI_recalc_single_EDTA_IDfixed_20260625/run/run_LAI_single_idfixed_worker.sh) | 运行：LAI single idfixed worker（据文件名） |
| [`03_panEDTA/EDTA_runner/panEDTA.sh`](03_panEDTA/EDTA_runner/panEDTA.sh) | This is the serial version of panEDTA |
| [`03_panEDTA/run_direct_homology_monitor.sh`](03_panEDTA/run_direct_homology_monitor.sh) | 运行：direct homology monitor（据文件名） |
| [`03_panEDTA/run_direct_homology_worker.sh`](03_panEDTA/run_direct_homology_worker.sh) | 运行：direct homology worker（据文件名） |
| [`03_panEDTA/run_panEDTA.sh`](03_panEDTA/run_panEDTA.sh) | 运行：panEDTA（据文件名） |
| [`03_panEDTA/run_panEDTA_resume_after_homology.sh`](03_panEDTA/run_panEDTA_resume_after_homology.sh) | 运行：panEDTA resume after homology（据文件名） |
| [`03_panEDTA/run_panEDTA_structural_dynamic_worker.sh`](03_panEDTA/run_panEDTA_structural_dynamic_worker.sh) | 运行：panEDTA structural dynamic worker（据文件名） |
| [`03_panEDTA/run_panEDTA_structural_dynamic_worker_fast.sh`](03_panEDTA/run_panEDTA_structural_dynamic_worker_fast.sh) | 运行：panEDTA structural dynamic worker fast（据文件名） |
| [`03_panEDTA/run_panEDTA_structural_monitor.sh`](03_panEDTA/run_panEDTA_structural_monitor.sh) | 运行：panEDTA structural monitor（据文件名） |
| [`03_panEDTA/run_panEDTA_structural_reconcile.sh`](03_panEDTA/run_panEDTA_structural_reconcile.sh) | 运行：panEDTA structural reconcile（据文件名） |
| [`03_panEDTA/run_panEDTA_structural_worker.sh`](03_panEDTA/run_panEDTA_structural_worker.sh) | 运行：panEDTA structural worker（据文件名） |
| [`05_softmask_panEDTA/regenerate_softmask_lowercase.sh`](05_softmask_panEDTA/regenerate_softmask_lowercase.sh) | 重新生成：softmask lowercase（据文件名） |
| [`05_softmask_panEDTA/run_softmask_after_ready.sh`](05_softmask_panEDTA/run_softmask_after_ready.sh) | 运行：softmask after ready（据文件名） |
| [`07_LAI_panEDTA_IDfixed_fullpass/collect_LAI_fullpass_summary.sh`](07_LAI_panEDTA_IDfixed_fullpass/collect_LAI_fullpass_summary.sh) | 收集结果：LAI fullpass summary（据文件名） |
| [`07_LAI_panEDTA_IDfixed_fullpass/fast_rerun_20260630_474ab/run_one_fast.sh`](07_LAI_panEDTA_IDfixed_fullpass/fast_rerun_20260630_474ab/run_one_fast.sh) | 运行：one fast（据文件名） |
| [`07_LAI_panEDTA_IDfixed_fullpass/prepare_idfixed_inputs_fullpass.sh`](07_LAI_panEDTA_IDfixed_fullpass/prepare_idfixed_inputs_fullpass.sh) | 准备输入：idfixed inputs fullpass（据文件名） |
| [`07_LAI_panEDTA_IDfixed_fullpass/run_LAI_fullpass_worker.sh`](07_LAI_panEDTA_IDfixed_fullpass/run_LAI_fullpass_worker.sh) | 运行：LAI fullpass worker（据文件名） |
| [`08_TEsorter_EDTA_final/run/build_manifest.sh`](08_TEsorter_EDTA_final/run/build_manifest.sh) | 构建：manifest（据文件名） |
| [`08_TEsorter_EDTA_final/run/collect_TEsorter_summary.py`](08_TEsorter_EDTA_final/run/collect_TEsorter_summary.py) | 收集结果：TEsorter summary（据文件名） |
| [`08_TEsorter_EDTA_final/run/TEsorter_worker.sh`](08_TEsorter_EDTA_final/run/TEsorter_worker.sh) | 命令：`NODE_LABEL=${1:-$(hostname)}` |
| [`make_singleEDTA_parallel_cmds.sh`](make_singleEDTA_parallel_cmds.sh) | 生成：singleEDTA parallel cmds（据文件名） |
| [`run_all_singleEDTA_cmd.sh`](run_all_singleEDTA_cmd.sh) | 运行：all singleEDTA cmd（据文件名） |
