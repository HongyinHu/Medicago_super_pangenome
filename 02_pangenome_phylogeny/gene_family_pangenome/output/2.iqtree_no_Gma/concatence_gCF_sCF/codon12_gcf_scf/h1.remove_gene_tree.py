# 删除基因树中的基因，保留物种名称

import re
import os 
import sys

def dels(i):
	res = re.search(r"(|.*?):",i).group(1)
	return res

def main(indir,outdir):
	for file in os.listdir(indir):
		if file.endswith("treefile"):
			fullfile = os.path.join(indir,file)
			with open(fullfile) as inf, open(f"{outdir}/{file}","w") as ouf:
				for lines in inf.readlines():
					line = re.sub(r"\|\w+[.]*[-]*\w+[.]*\w+[.]*\w+","",lines)
					print(line,end="",file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s indir outdir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])

