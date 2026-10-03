#!/usr/bin/env python3
from pathlib import Path
import re

base = Path('path/to/project/6.genome_quality_assess/output/X_genome_M46_2')
flagstat = base / 'assembly_evaluate/ReadMap/flagstat.txt'
text = flagstat.read_text(errors='replace')
def find(pattern, label):
    m = re.search(pattern, text, flags=re.MULTILINE)
    if not m:
        raise ValueError(f'Could not parse {label} from flagstat.txt')
    return m.groups()
total, = find(r'^(\d+) \+ \d+ in total', 'total reads')
mapped, mapped_pct = find(r'^(\d+) \+ \d+ mapped \(([0-9.]+)%', 'mapped reads')
proper, proper_pct = find(r'^(\d+) \+ \d+ properly paired \(([0-9.]+)%', 'properly paired reads')
total, mapped, proper = map(int, (total, mapped, proper))
report = base / 'assessment_report.md'
s = report.read_text()
s = s.replace('状态：BUSCO 与 Merqury 已完成；survey reads 的 BWA 比对仍在运行，待 flagstat 生成后补入最终读段比对指标。本报告中的其他指标来自已完成的结果文件。', '状态：组装统计、组装与注释 BUSCO、Merqury 和 survey reads 比对均已完成。')
s = s.replace('| Survey reads 比对 | 运行中 | 完成后补入总 reads、mapped 与 properly paired 比例 |', f'| Survey reads 比对 | {mapped_pct}% mapped；{proper_pct}% properly paired | {total:,} reads；mapped {mapped:,}；properly paired {proper:,} |')
s = s.replace('| Survey reads BWA | 运行中 | 99.71% mapped；99.16% properly paired | M46_2 结果待补齐 |', f'| Survey reads BWA | {mapped_pct}% mapped；{proper_pct}% properly paired | 99.71% mapped；99.16% properly paired | 同为 BWA 0.7.17 与 survey reads 全量比对 |')
s = s.replace('当前尚未结束。', '已完成；flagstat 统计基于成对 survey FASTQ。')
s = s.replace('当前等待 `flagstat.txt` 完成', '最终统计见 `flagstat.txt`')
report.write_text(s)
(base / 'assessment_status.txt').write_text(
    'Assessment complete.\n'
    f'BWA survey reads: {total} total; {mapped} mapped ({mapped_pct}%); {proper} properly paired ({proper_pct}%).\n'
    'BUSCO assembly: C96.9%, S88.2%, D8.7%, F0.7%, M2.4%.\n'
    'BUSCO annotation: C94.7%, S85.3%, D9.4%, F1.1%, M4.2%.\n'
    'Merqury k=21: QV40.2527; completeness63.779%.\n'
)
(base / 'comparison_summary.tsv').write_text(
    'metric\tM46_2\tM22_0\tnote\n'
    'assembly_BUSCO_C\t96.9% (S88.2,D8.7)\t99.4% (S95.4,D4.0)\tBUSCO 5.6.0 vs 5.4.3; same 2020 lineage\n'
    'annotation_BUSCO_C\t94.7% (S85.3,D9.4)\t94.7% (S91.0,D3.7)\tBUSCO versions differ; same lineage\n'
    f'survey_reads_mapped\t{mapped_pct}%\t99.71%\tBWA 0.7.17, survey reads\n'
    f'survey_reads_properly_paired\t{proper_pct}%\t99.16%\tBWA 0.7.17, survey reads\n'
    'Merqury_QV\t40.2527 (k=21)\t48.283 (k=19)\tk differs; not directly comparable\n'
    'Merqury_completeness\t63.779% (k=21)\t98.8316% (k=19)\tk and filtering differ; not directly comparable\n'
)
print(f'Finalized report: total={total}, mapped={mapped_pct}%, properly_paired={proper_pct}%')