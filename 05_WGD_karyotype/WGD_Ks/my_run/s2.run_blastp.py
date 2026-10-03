# 获取蛋白序列的blastp比对结果

import re
import os 
import sys

script1 = "path/to/home/sofeware/CEGMA_v2/depend_sofe/blast+/ncbi-blast-2.10.0+/bin/makeblastdb"
script2 = "path/to/home/sofeware/CEGMA_v2/depend_sofe/blast+/ncbi-blast-2.10.0+/bin/blastp"

def main(pep_fa,prefix):
	pep_fa = os.path.abspath(pep_fa)
	# cds_fa = os.path.abspath(cds_fa)

	# cmd1 = 'seqkit grep -f <(cut -f 7 {gff} ) {cds_fa} | seqkit seq --id-regexp \"^(.*?)\\.\\d\" -i > {prefix}.cds'.format(gff=gff, cds_fa=cds_fa,prefix=prefix)
	# cmd2 = 'seqkit grep -f <(cut -f 7 {gff} ) {pep_fa} | seqkit seq --id-regexp \"^(.*?)\\.\\d\" -i > {prefix}.pep'.format(gff=gff, pep_fa=pep_fa,prefix=prefix)
	cmd3 = "{script1} -in {prefix}.pep -dbtype prot".format(script1=script1, prefix=prefix)
	cmd4 = "{script2} -num_threads 20 -db {prefix}.pep -query {prefix}.pep -outfmt 6 -evalue 1e-5 -num_alignments 20  -out {prefix}.blastp.txt".format(
		script2=script2, prefix=prefix)

	# print("\n".join([cmd1, cmd2, cmd3, cmd4]))
	print("\n".join([cmd3, cmd4]))

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s pep_fa prefix"%sys.argv[0])
		sys.exit(0)

	main(sys.argv[1], sys.argv[2])
