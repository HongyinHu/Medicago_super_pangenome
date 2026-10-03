# 修改格式Orthotsv适应PanGP软件

import re
import os 
import sys

def main(OrthoTsv):

	key_list = list()
	with open(OrthoTsv) as inf:
		for index,lines in enumerate(inf.readlines()):
			
			
			if index == 0:
				line = lines.strip().split()
				for i in line[1:]:
					key_list.append(i.strip())

				print("\t".join(line))

			else:
				line = re.split(r"[\s+,]",lines)
				
				line_dict = dict()

				for i in key_list:
					line_dict.setdefault(i,[])

				for i in line[1:]:
					if i:
						key_name = i.strip().split("|")[0]
						if key_name in key_list:
							line_dict[key_name].append(i)

				print(line[0],end="")
				for k,v in line_dict.items():
					if len(v) == 0:
						v = "-"
						print("\t"+v,end="")
					else:
						print("\t"+",".join(v),end="")

				print()


if __name__ == '__main__':
	main(sys.argv[1])


