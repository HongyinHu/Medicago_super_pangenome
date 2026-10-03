# 以下是一个Python脚本，可以从OrthoFinder结果中区分泛基因组中的四类基因（core, soft-core, shell和specific）：
import os
import re
import sys

def write_file(name,lis):
	with open("%s.txt"%name,"w") as ouf:
		for i in lis:
			print(i,file=ouf)

def species_class(OrthoTsv,lis,file_name):
	species_dict = dict()
	with open(OrthoTsv) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index > 0:
				line = lines.strip().split("\t")
				
				if line[0] in lis:
					items = re.split(r"[,\s]",lines)
					for item in items:
						if re.search(r"\|",item):
							spec_name, gene_name = re.split(r"\|",item)
							species_dict.setdefault(spec_name,[]).append(gene_name)

	if not os.path.exists("out_class"):
		os.mkdir("out_class")

	for k,v in species_dict.items():
		with open("%s/%s_%s.txt"%("out_class",k,file_name),"w") as ouf:
			for i in v:
				print(i,file=ouf)

		print("\t".join(map(str,[k,file_name,len(v)])))


def main(geneCount,OrthoTsv):
	core_list = list()
	softcore_list = list()
	shell_list = list()
	specific_list = list()


	with open(geneCount) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index > 0:
				line = lines.strip().split()
				line_num = 0
				
				for i in line[1:-1]:
					if int(i) > 0:
						line_num += 1

				if line_num == 18:
					core_list.append(line[0])
				elif 15 <= line_num <= 17:
					softcore_list.append(line[0])
				elif 2 <= line_num <= 14:
					shell_list.append(line[0])
				else:
					specific_list.append(line[0])

	write_file("pangenome_Core",core_list)
	write_file("pangenome_SoftCore",softcore_list)
	write_file("pangenome_Shell",shell_list)
	write_file("pangenome_Specific",specific_list)

	# 对物种分类
	species_class(OrthoTsv,core_list,"Core")
	species_class(OrthoTsv,softcore_list,"SoftCore")
	species_class(OrthoTsv,shell_list,"Shell")
	species_class(OrthoTsv,specific_list,"Specific")

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s geneCount_tsv Orthogroup_tsv"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])





