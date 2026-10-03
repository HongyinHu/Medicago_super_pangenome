"""
1. 将gblock过滤后的多序列比对结果，进一步使用trimal进行过滤
"""
import re
import os 
import sys

def main(indir,outdir):
	indir = os.path.abspath(indir)
	outdir = os.path.abspath(outdir)

	for file in os.listdir(indir):
		if file.endswith("align"):
			fullfile = os.path.join(indir,file)
			newfile = re.sub("align","aln-cln",file)
			print("trimal -in %s -out %s/%s -gt 0.8"%(fullfile,outdir,newfile))

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("python %s gblock_indir trimal_outdir\ngblock_indir--gblock文件夹\ntrimal_outdir--trimal输出文件夹")
		sys.exit(0)
	
	main(sys.argv[1],sys.argv[2])
