# 过滤后的Chip seq数据比对到参考基因组
# -k 10：每个 read pair 最多报告 10 个 valid alignments

import re
import os
import sys

def main(clean_fq1, clean_fq2, ref_genome, sample, threads):
	cmd = (f"bowtie2-build --threads 10 {ref_genome} {sample} && "
		   f"bowtie2 -p {threads} -x {sample} -1 {clean_fq1} -2 {clean_fq2} -k 10 -X 2000 | "
	       f"samtools view -@ {threads} -b -o {sample}.k10.bam && "
	       f"samtools sort -@ {threads} -o {sample}.k10.sort.bam {sample}.k10.bam && " 
	       f"samtools index -@ {threads} {sample}.k10.sort.bam")

	print(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 6:
		print("python %s clean_fq1 clean_fq2 ref_genome sample_name threads"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])


