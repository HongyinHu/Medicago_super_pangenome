# 提取ka/ks结果

import re
import os 
import sys

def main(indir,prefix):
	with open("kaks.list","w") as ouf:
		for dirpath,dirnames,files in os.walk(indir):
			for file in files:
				if file.endswith("kaks") and not file.startswith(".") and not file.startswith("msg"):
					fullfile = os.path.join(dirpath,file)
					if os.path.getsize(fullfile) != 0:
						with open(fullfile) as inf:
							for index,lines in enumerate(inf.readlines()):
								if index == 1:
									line = lines.strip().split()
									if line[4] == "NA": continue
									if float(line[4]) < 10:
										print(prefix,line[4],sep="\t",file=ouf)
									# try:
									# 	line = lines.strip().split()
									# 	if int(line[4]) == 50:
									# 		print(fullfile)
									# 		exit(0)
									# 	print(prefix,line[4],sep="\t",file=ouf)
									# except:
									# 	pass
					# else:
					# 	print(fullfile)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s res_indir prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])
