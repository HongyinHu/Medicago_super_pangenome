# 获取wgdi要求的gff格式

import re
import os 
import sys
from Bio import SeqIO

# 不同基因组需要修改，去除contig
# chr_prefix = "Chr"
# chr_prefix = "chr"
chr_prefix = "MtrunA17Chr"

chr_name = list()
for i in range(1,9):
	chr_name.append(chr_prefix+str(i))


# 基本不用修改
def main(genome_fa,genome_gff,genome_cds,genome_pep,prefix):
	chr_len = dict()
	for seq in SeqIO.parse(genome_fa,"fasta"):
		if seq.id in chr_name:
			if re.search("Chr",seq.id):
				seqid = re.search(r"Chr(\d+)",seq.id).group(1)
				chr_len[seqid] = len(seq.seq)

	chr_gene = dict()
	flag = 1
	flag_name = None

	gene_name_change = dict()
	with open(genome_gff) as inf, open(f"{prefix}.gff","w") as ouf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split("\t")
				# if line[0] in chr_name:
				if not re.search("Chr",line[0]): continue
				
				if flag_name != line[0]:
					flag_name = line[0]
					flag = 1

				if line[2] == "mRNA":
					# mraname = re.split(':|=|;',line[8])[1]
					mraname = re.search(r"ID=(.*?);",line[8]).group(1)
					genename = re.split(':|=|;',line[8])[3]

					line_chr = re.search(r"Chr(\d+)",line[0]).group(1)
					chr_gene.setdefault(line_chr,[]).append(line[8])
					line[0] = re.sub('Chr|chr|MtrunA17Chr',"",line[0])

					# dict oldname = newname
					newname = f"{prefix}{line[0]}g%05d"%flag
					gene_name_change[mraname] = newname
					
					res_lis = [line_chr,newname,line[3],line[4],line[6],flag,genename]
					print("\t".join(map(str,res_lis)),file=ouf)
					flag += 1

	with open(f"{prefix}.lens","w") as ouf:
		for k,v in chr_len.items():
			ks = re.sub('Chr|chr|MtrunA17Chr',"",k)
			print("\t".join(map(str,[k,v,len(chr_gene[k])])),file=ouf)

	with open(f"{prefix}.cds","w") as ouf:
		for seq in SeqIO.parse(genome_cds,"fasta"):
			if seq.id in gene_name_change.keys():
				print(f">{gene_name_change[seq.id]}\n{seq.seq}",file=ouf)

	with open(f"{prefix}.pep","w") as ouf:
		for seq in SeqIO.parse(genome_pep,"fasta"):
			if seq.id in gene_name_change.keys():
				seq_seq = re.sub(r"\.","",str(seq.seq))
				print(f">{gene_name_change[seq.id]}\n{seq_seq}",file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 6:
		print("python %s genome_fa genome_gff genome_cds genome_pep prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])



