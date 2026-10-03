# 使用mcscan绘制共线性图

import re
import os 
import sys

def get_bed(gff3, prefix):
	with open(gff3) as inf, open("%s.bed"%prefix,"w") as ouf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split("\t")
				if line[2] == "mRNA":
					seqid = re.search(r"ID=(.*?);",line[-1]).group(1)
					outline = "\t".join(map(str,[line[0],line[3],line[4], seqid, line[-2], line[-3]]))
					print(outline,file=ouf)


def main(gff3, cds_fa, pep_fa, prefix):

	# get bed from gff3
	# print("python -m jcvi.formats.gff bed --type=mRNA --key=ID {gff3} -o {prefix}.bed".format(gff3=gff3, prefix=prefix))
	get_bed(gff3,prefix)

	# get uniq bed
	print("python -m jcvi.formats.bed uniq {prefix}.bed".format(prefix=prefix))

	# get cds and pep fasta file
	print("seqkit grep -f <(cut -f4 {prefix}.uniq.bed) {cds_fa}  | seqkit seq -i > {prefix}.cds".format(prefix=prefix, cds_fa=cds_fa))
	print("seqkit grep -f <(cut -f4 {prefix}.uniq.bed) {pep_fa}  | seqkit seq -i > {prefix}.pep".format(prefix=prefix, pep_fa=pep_fa))

if __name__ == '__main__':
	if len(sys.argv) != 5:
		print("python %s gff3 cds_fa pep_fa prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
