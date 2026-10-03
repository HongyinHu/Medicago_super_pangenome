#从cactus多序列比对结果提取物种物种特异性片段的bed文件
import re
import os
import sys

def main(maf_file,prefix):
	with open(maf_file) as inf:
		content = inf.read().split("\n\n")
		for i in content:
			species = set()

			if re.search(prefix,i):
				for lines in i.split("\n"):
					if lines.startswith("s"):
						line = lines.strip().split()
						name = re.split(r"\.",line[1])[0]
						if not re.search("Anc",name):
							species.add(name)

			if len(species) == 1:
				# print(i)
				# print("##"*20)
				for lines in i.split("\n"):
					if re.search(prefix,lines):
						line = lines.strip().split()
						Chr = re.split(r"\.",line[1])[1]
						END = int(line[2]) + int(line[3])
						print("\t".join([Chr,line[2],str(END),line[3],line[4]]))


if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s maf_file prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])