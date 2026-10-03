# 使用orthofinder聚类的同源基因家族簇，计算四类范基因组四类基因集的ka/ks
# perl path/to/home/sofeware/ParaAT2.0/ParaAT.pl -h ./res_homologs.txt -n ./res_cds.fa -a ./res_pep.fa  -p proc -m mafft -f axt -g -k -o result_dir


import re
import os 
import sys
from Bio import SeqIO

def main(all_pep,all_cds,group_csv,genefam):
	pan_gene_fam = list()
	with open(genefam) as inf:
		for lines in inf.readlines():
			pan_gene_fam.append(lines.strip())

	ortho_genes = dict()
	with open(group_csv) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index > 0:
				line = re.split(r"\s+|,",lines)
				line_new = [i for i in line if i != ""]
				gene_lis = list()
				for i in line_new[1::]:
					# gene_name = re.split(r"\|",i)[1]
					gene_name = i
					gene_lis.append(gene_name)

				ortho_genes.setdefault(line_new[0],[])
				ortho_genes[line_new[0]] += gene_lis

	

	res_dict = dict()
	for i in pan_gene_fam:
		res_dict[i] = ortho_genes[i]

	with open("res_homologs.txt","w") as ouf:
		for k,v in res_dict.items():
			print("\t".join(v),file=ouf)

	pep_dict = dict()
	cds_dict = dict()
	for seq in SeqIO.parse(all_pep,"fasta"):
		pep_dict[seq.id] = seq.seq

	for seq in SeqIO.parse(all_cds,"fasta"):
		cds_dict[seq.id] = seq.seq

	with open("res_pep.fa","w") as ouf1, open("res_cds.fa","w") as ouf2:
		for k,v in res_dict.items():
			for i in v:
				print(">"+i+"\n"+pep_dict[i],file=ouf1)
				print(">"+i+"\n"+cds_dict[i],file=ouf2)


if __name__ == '__main__':
	if len(sys.argv) != 5:
		print("python %s all_pep all_cds orthogroup_tsv pan_class_genelist"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])



