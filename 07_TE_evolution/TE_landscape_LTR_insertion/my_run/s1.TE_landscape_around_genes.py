# 蛋白编码基因前后10kp附近，各个类型TE的分布情况
import os
import sys
import re
from Bio import SeqIO
from multiprocessing import Pool

def count_TE_types(lis):
	lis1,lis2 = lis
	counts = {}
	
	for i in lis1:
		upstream = (i[0]-10000,i[0])
		downstream = (i[1],i[1]+10000)
		for ii in lis2:
			if ii[0] >= upstream[0] and ii[1] <= upstream[1]:
				ave_te = ii[0] + (ii[1]-ii[0]) / 2 
				dis = ii[1] - upstream[1]
				counts.setdefault("upstream",{}).setdefault(ii[2],[]).append((dis,i[-1]))
			# elif ii[0] < upstream[0] and ii[1] >= upstream[0]  and ii[1] <= upstream[1]:
			# 	ave_te = ii[0] + (ii[1]-ii[0]) / 2 
			# 	dis = ave_te - upstream[1]
			# 	counts.setdefault("upstream",{}).setdefault(ii[2],[]).append((dis,i[-1]))

			elif ii[0] >= downstream[0] and ii[1] <= downstream[1]:
				ave_te = ii[0] + (ii[1]-ii[0]) / 2 
				dis = ii[0] - downstream[0]
				counts.setdefault("downstream",{}).setdefault(ii[2],[]).append((dis,ii[-1]))
			# elif ii[0] >= downstream[0] and ii[0] <= downstream[1]  and ii[1] >= downstream[1]:
			# 	ave_te = ii[0] + (ii[1]-ii[0]) / 2 
			# 	dis = ave_te - downstream[0]
			# 	counts.setdefault("upstream",{}).setdefault(ii[2],[]).append((dis,i[-1]))
	
	return counts

def main(gene_annotation_file, TE_annotation_file,genome_fa):
	gene_annotation = dict()
	TE_annotation = dict()

	genome_dict = dict()
	for seq in SeqIO.parse(genome_fa,"fasta"):
		genome_dict[seq.id] = len(seq.seq)

	with open(gene_annotation_file) as inf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split("\t")
				line[3] = int(line[3])
				line[4] = int(line[4])
				if line[2] == "mRNA":
					# if line[6] == "-":
					# 	start = genome_dict[line[0]] - line[4] + 1
					# 	end = genome_dict[line[0]] - line[3] + 1
					# else:
					start = line[3]
					end = line[4]
					gene_annotation.setdefault(line[0],[]).append((start,end,line[0]))

	with open(TE_annotation_file) as inf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split("\t")
				line[3] = int(line[3])
				line[4] = int(line[4])
				Class = re.search("Classification=(.*?);",line[-1]).group(1)
				Class = re.sub(r"/","_",Class)

				# if line[6] == "-":
				# 	start = genome_dict[line[0]] - line[4] + 1
				# 	end = genome_dict[line[0]] - line[3] + 1
				# else:
				start = line[3]
				end = line[4]
				
				TE_annotation.setdefault(line[0],[]).append((start,end,Class,line[0]))

	pool = Pool(20)
	para = list()
	for k in genome_dict.keys():
		if k in gene_annotation and k in TE_annotation:
			para.append((gene_annotation[k],TE_annotation[k]))

	# para.append((gene_annotation["Chr1"],TE_annotation["Chr1"]))
	
	print(len(para))
	res = pool.map(count_TE_types,para)

	for i in res:
		for k,v in i.items():
			for k2,v2 in v.items():
				if re.search("Helitro",k2):
					k2 = "Helitro"
				if re.search("DNA|MITE",k2):
					k2 = "DNA_TIRs"
				if re.search("Low_complexity",k2):
					k2 = "Simple_repeat"
				
				if k == "upstream":
					with open(f"upstream/{k2}.txt","a+") as ouf:
						for i in v2:
							print("\t".join([k2,str(i[0]/1000)]),file=ouf)
				else:
					with open(f"downstream/{k2}.txt","a+") as ouf:
						for i in v2:
							print("\t".join([k2,str(i[0]/1000)]),file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s gene_annotation_file TE_annotation_file genome_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3])				




















