# 过滤溯祖树比对序列
# 删除基因名，保留物种名称

import re
import os
import sys
from Bio import SeqIO

def main(indir,outdir):
	for file in os.listdir(indir):
		if file.endswith("aln-cln"):
			fullfile = os.path.join(indir,file)

			content_list=dict()
			for seq in SeqIO.parse(fullfile,"fasta"):
				if str(seq.seq).count("-")/len(seq.seq) < 0.2:
					newid = re.sub(r"\|.*","",seq.id)
					content_list[newid]=str(seq.seq)

			if len(content_list) == 17:
				with open(f"{outdir}/{file}","w") as ouf:
					for k,v in content_list.items():
						print(f">{k}\n{v}",file=ouf)


if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s indir outdir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])
