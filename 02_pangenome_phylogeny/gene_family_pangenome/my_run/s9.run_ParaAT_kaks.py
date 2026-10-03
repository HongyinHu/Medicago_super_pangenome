# 使用ParaAT批量计算kaks
# perl path/to/home/sofeware/ParaAT2.0/ParaAT.pl -h ./res_homologs.txt -n ./res_cds.fa -a ./res_pep.fa  -p proc -m mafft -f axt -g -k -o result_dir
import re
import os 
import sys

script = "perl path/to/home/sofeware/ParaAT2.0/ParaAT.pl"

def full_dir(file):
	fullfile = os.path.abspath(file)
	return fullfile

def main(homologs,cds,pep):
	cds = full_dir(cds)
	pep = full_dir(pep)
	# proc = full_dir(proc)

	if not os.path.exists("homologs_split"):
		os.mkdir("homologs_split")

	if not os.path.exists("homologs_result"):
		os.mkdir("homologs_result")

	os.system("cp proc homologs_result")

	with open(homologs) as inf:
		lines = inf.readlines()

		n = 500
		for i in range(0, len(lines), n):
			chunk = lines[i:i+n]
			with open(f'homologs_split/homologs_{i//n}.txt', 'w') as f:
				f.writelines(chunk)

		if len(lines) % n != 0:
			chunk = lines[-(len(lines) % n):]
			with open(f'homologs_split/homologs_{len(lines)//n}.txt', 'w') as f:
				f.writelines(chunk)

	for file in os.listdir("homologs_split"):
		fullfile = os.path.join("homologs_split",file)
		res_out = os.path.join("homologs_result",re.sub("txt","out",file))
		cmd = "%s -h %s -n %s -a %s -p proc -m mafft -f axt -k -o %s"%(script,fullfile,cds,pep,res_out)
		print(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s homologs cds pep"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3])