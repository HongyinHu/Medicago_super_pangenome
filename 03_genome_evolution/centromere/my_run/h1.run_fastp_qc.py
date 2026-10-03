#fastp 质控和接头去除, fastp 的 -q 是“合格碱基的 Phred 阈值

import re
import os
import sys

def main(fq1, fq2, threads):
	o_fq1 = re.sub("fq.gz","clean.fq.gz",fq1)
	o_fq2 = re.sub("fq.gz","clean.fq.gz",fq2)
	cmd = f"fastp -w {threads} -i {fq1} -I {fq2} -o {o_fq1} -O {o_fq2} -q 20"
	print(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s fq1 fq2 threads"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3])