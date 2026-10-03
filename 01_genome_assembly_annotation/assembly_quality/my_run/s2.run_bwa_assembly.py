# 评估基因组组装的完整性 使用bwa-mem将二代基因组reads比对组装序列上

import re
import os 
import sys

script = "path/to/software/bwa/bwa-0.7.17/bwa"

def main(assembly_fa, paired_read1, paired_read2):
	cmd1 = f"{script} index -a bwtsw {assembly_fa} -p genome"
	print(cmd1)

	cmd2 = f"{script} mem genome {paired_read1} {paired_read2} -o out.sam -t 20"
	print(cmd2)

	cmd3 = f"samtools flagstat out.sam > flagstat.txt"
	print(cmd3)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s assembly_fa paired_read1 paired_read2"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2], sys.argv[3])