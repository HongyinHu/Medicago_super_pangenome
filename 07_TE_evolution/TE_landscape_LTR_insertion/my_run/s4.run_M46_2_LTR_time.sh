#!/bin/bash
# Recompute intact-LTR K values / insertion times for the new M46 assembly (M46_2),
# using the same pipeline as 17.TE_landscape_2/output_LTR_time (p12 -> p13 -> MEGA K2P -> p11)
# note: 32-bit megacc prints "libgcc_s.so.1 must be installed" and exits non-zero AFTER writing its .meg,
#       so ParaFly reports every command as failed; outputs are still valid (same as the 2024 runs)
export PATH=path/to/home/anaconda3/bin:path/to/home/bin:$PATH
B=path/to/project/17.TE_landscape_2
R=$B/my_run
N=path/to/projects_all/pan_genome/X_genome_M46_2/5.genome_annotation_repeat_predict/output/EDTA3
O=$B/output_LTR_time_M46_2
rm -rf $O/Mfi_46.intact_LTR
mkdir -p $O/Mfi_46.intact_LTR
cd $O
ln -sf $N/M46_genome_Chr_reorder.genome.ctg.fa Mfi_46.genome.fa
cd Mfi_46.intact_LTR
ln -sf $N/M46_genome_Chr_reorder.genome.ctg.fa.mod.EDTA.raw/M46_genome_Chr_reorder.genome.ctg.fa.mod.LTR.intact.raw.gff3 M46_2.mod.LTR.intact.gff3
# new EDTA writes "classification=" (lowercase); old p12 matches "Classification=" -> patched copy, otherwise identical
sed 's/Classification=(\.\*?);",lines)/[Cc]lassification=(.*?);",lines)/g' $R/p12.get_intact_LTR_type.py > p12.get_intact_LTR_type.M46_2.py
python p12.get_intact_LTR_type.M46_2.py M46_2.mod.LTR.intact.gff3 Mfi_46
for t in Copia Gypsy unknown; do
  echo "[$(date)] $t start"
  python $R/p13.get_shared_intact_TE_fa_type.py ../Mfi_46.genome.fa Mfi_46_LTR_${t}.lLTR.bed Mfi_46_LTR_${t}.rLTR.bed $t
  python $R/p10.cal_K_value.py LTR_$t 1; ParaFly -c LTR_$t.cmd1.sh -CPU 10 > /dev/null 2>&1
  python $R/p10.cal_K_value.py LTR_$t 2; ParaFly -c LTR_$t.cmd2.sh -CPU 10 > /dev/null 2>&1
  python $R/p11.get_k_Value_from_result.py LTR_${t}_res2 > Mfi_46.intact_LTR_${t}.shared.txt
  cut -f 3 Mfi_46.intact_LTR_${t}.shared.txt | sed 's/$/\tMfi_46/' > Mfi_46.intact_LTR_${t}.shared.time.txt
  echo "[$(date)] $t: $(ls LTR_$t/*.fa | wc -l) LTR pairs, $(ls LTR_${t}_res1/*.meg | wc -l) aligned, $(wc -l < Mfi_46.intact_LTR_${t}.shared.txt) with K"
done
echo DONE
