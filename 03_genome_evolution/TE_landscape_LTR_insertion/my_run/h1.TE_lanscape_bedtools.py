# 蛋白编码基因前后10kp附近,各个类型TE的分布情况
# HHY 20241017
# bedtools closest -a LTR_Copia.TE.bed -b genename.gene.bed -D a  > LTR_Copia.TE.bed.txt 

import os
import sys
import re
from Bio import SeqIO

def main(gene_annotation_file, TE_annotation_file, genome_fa):
	gene_annotation = list()
	TE_annotation = dict()

	with open("TE/genome.info.bed",'w') as ouf:
		for seq in SeqIO.parse(genome_fa,"fasta"):
			print("\t".join([seq.id,str(len(seq.seq))]),file=ouf)

	TE_flag = 1
	with open(TE_annotation_file) as inf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split("\t")
				if not re.search("scaffold",line[0]) and not re.search("Chr0c",line[0]):
					start = line[3]
					end = line[4]
					stand = line[6]
					Class = re.search("Classification=(.*?);",line[-1]).group(1)
					Class = re.sub(r"/","_",Class)
					if re.search("Helitro",Class):
						Class = "Helitro"
					if re.search("DNA|MITE",Class):
						Class = "DNA_TIRs"
					if re.search("Low_complexity",Class):
						Class = "Simple_repeat"
					TE_name = f"{Class}{TE_flag}"
					TE_flag += 1
					TE_annotation.setdefault(Class,[]).append((line[0],start,end,TE_name,stand))

	TE_file_list = list()
	for k,v in TE_annotation.items():
		TE_file_list.append(f"{k}.TE.bed")
		with open(f"TE/{k}.TE.bed","w") as ouf:
			for i in v:
				print("\t".join(i),file=ouf)

	gene_flag = 1
	with open(gene_annotation_file) as inf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split("\t")
				if line[2] == "mRNA":
					start = line[3]
					end = line[4]
					stand = line[6]
					genename = f"gene{gene_flag}"
					gene_flag += 1
					if not re.search("scaffold",line[0])  and not re.search("Chr0c",line[0]):
						gene_annotation.append((line[0],start,end,genename,stand))
					
	with open(f"TE/gene.info.bed","w") as ouf:
		for i in gene_annotation:
			print("\t".join(i),file=ouf)

	with open("TE/bedtools.run.sh","w") as ouf:
		for i in TE_file_list:
			print(f"bedtools closest -a {i} -b gene.info.bed -D a > {i}.out",file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s gene_annotation_file TE_annotation_file genome_fa")
	main(sys.argv[1],sys.argv[2],sys.argv[3])