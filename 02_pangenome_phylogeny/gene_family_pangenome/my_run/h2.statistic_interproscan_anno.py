# 统计interproscan注释结果

import re
import os 
import sys

def main(indir,file_class):
	anno_gene_list = set()
	indir = os.path.abspath(indir)
	for dirpath,dirnames,files in os.walk(indir):
		for file in files:
			if file.endswith("fa.out"):
				fullfile = os.path.join(dirpath,file)

				with open(fullfile) as inf:
					for lines in inf.readlines():
						line = lines.strip().split("\t")
						line[0] = re.sub(":","_",line[0])
						anno_gene_list.add(line[0])

	class_genes = set()
	with open(file_class) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			for i in line:
				class_genes.add(i)

	anno_list = set()
	unanno_list = set()
	for i in class_genes:
		if i in anno_gene_list:
			anno_list.add(i)
		else:
			unanno_list.add(i)

	anno_num = len(anno_list)
	unanno_num = len(unanno_list)

	with open("%s.anno.txt"%file_class,"w") as ouf1, open("%s.unanno.txt"%file_class,"w") as ouf2:
		for i in anno_list:
			print(i,file=ouf1)

		for ii in unanno_list:
			print(ii,file=ouf2)


	print("all gene in class_file %d"%len(class_genes))
	print("anno numer %d"%anno_num)
	print("unanno number %d"%unanno_num)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s indir_interpro class_gene_file"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])

