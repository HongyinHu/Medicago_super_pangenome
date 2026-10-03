#获取每个两两成组结果的共线性比例
#time:20241125
#author:HHY

import re
import os
import sys
from Bio import SeqIO

def main(indir):
	indir = os.path.abspath(indir)

	for dirpath,dirnames,files in os.walk(indir):
		for dirname in dirnames:
			if dirname == "70":
				all_junction_dict = dict()
				all_match = dict()

				full_dirname = os.path.join(dirpath,dirname)

				name1 = os.path.basename(os.path.dirname(full_dirname))
				name2 = "_".join(re.split("_",name1)[0:2])
				name3 = "_".join(re.split("_",name1)[2:4])

				full_file = os.path.join(dirpath,"ltr_juction.fa")
				for seq in SeqIO.parse(full_file,"fasta"):
					all_junction_dict[seq.id] = str(seq.seq)

				for file in os.listdir(full_dirname):
					if file.endswith("match"):
						newfile = re.sub("match","fna",file)
						full_file_name = os.path.join(full_dirname, newfile)
						all_seq_id = list()
						for seq in SeqIO.parse(full_file_name,"fasta"):
							all_seq_id.append(seq.id)
						
						if re.search(name2,"".join(all_seq_id)) and re.search(name3,"".join(all_seq_id)):
							for seq in SeqIO.parse(full_file_name,"fasta"):
								all_match[seq.id]=str(seq.seq)
				
				name2_junction_num, name3_junction_num = 0, 0
				name2_match_num, name3_match_num = 0, 0

				for k in all_junction_dict:
					if re.search(name2,k):
						name2_junction_num += 1
					elif  re.search(name3,k):
						name3_junction_num += 1

				for k in all_match:
					if re.search(name2,k):
						name2_match_num += 1
					elif  re.search(name3,k):
						name3_match_num += 1

				res_name2 = "%.2f"%((name2_match_num/name2_junction_num)*100)
				res_name3 = "%.2f"%((name3_match_num/name3_junction_num)*100)
				print("\t".join([name2,name3,str(res_name2)]))
				print("\t".join([name3,name2,str(res_name3)]))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s indir_groups"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])







