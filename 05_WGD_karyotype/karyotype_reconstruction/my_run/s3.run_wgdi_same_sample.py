# 书写wgdi配置文件和运行各个命令

import re
import os 
import sys

def run_dotplot(blastp_result, gff1, gff2, lens1, lens2, genome1_name, genome2_name, prefix):
	conf1 = """
[dotplot]
blast = {blastp_result}
gff1 =  {gff1}
gff2 =  {gff2}
lens1 = {lens1}
lens2 = {lens2}
genome1_name =  {genome1_name}
genome2_name =  {genome2_name}
multiple  = 1
score = 100
evalue = 1e-5
repeat_number = 30
position = order
blast_reverse = false
ancestor_left = none
ancestor_top = none
markersize = 0.5
figsize = 10,10
savefig = {prefix}.png
""".format(blastp_result=blastp_result, gff1=gff1, gff2=gff2, lens1=lens1, lens2=lens2, genome1_name=genome1_name, 
	genome2_name=genome2_name, prefix=prefix)
	
	cmd1 = "wgdi -d total.conf"
	return (conf1, cmd1)

def run_collinearity(blastp_result, gff1, gff2, lens1, lens2, prefix):
	conf2="""
[collinearity]
gff1 = {gff1}
gff2 = {gff2}
lens1 = {lens1}
lens2 = {lens2}
blast = {blastp_result}
blast_reverse = false
multiple  = 1
process = 8
evalue = 1e-5
score = 100
grading = 50,40,25
mg = 40,40
pvalue = 0.2
repeat_number = 10
position = order
savefile = {prefix}.collinearity.txt
""".format(blastp_result=blastp_result, gff1=gff1, gff2=gff2, lens1=lens1, lens2=lens2, prefix=prefix)
	cmd2 = "wgdi -icl total.conf"
	return (conf2, cmd2)

def run_ks(cds_fa, pep_fa, prefix):
	conf3="""
[ks]
cds_file = {cds_fa}
pep_file = {pep_fa}
align_software = muscle
pairs_file = {prefix}.collinearity.txt
ks_file = {prefix}.ks.txt
""".format(cds_fa=cds_fa, pep_fa=pep_fa, prefix=prefix)
	cmd3 = "wgdi -ks total.conf"
	return (conf3, cmd3)

def run_blockinfo(blastp_result, gff1, gff2, lens1, lens2, prefix):
	conf4 = """
[blockinfo]
blast = {blastp_result}
gff1 =  {gff1}
gff2 =  {gff2}
lens1 = {lens1}
lens2 = {lens2}
collinearity = {prefix}.collinearity.txt
score = 100
evalue = 1e-5
repeat_number = 30
position = order
ks = {prefix}.ks.txt
ks_col = ks_NG86
savefile = {prefix}_block_information.csv
""".format(blastp_result=blastp_result, gff1=gff1, gff2=gff2, lens1=lens1, lens2=lens2, prefix=prefix)
	cmd4 = "wgdi -bi total.conf"
	return (conf4, cmd4)

def run_correspondence(prefix, lens1, lens2):
	conf5 = """
[correspondence]
blockinfo =  {prefix}_block_information.csv
lens1 = {lens1}
lens2 = {lens2}
tandem = true
tandem_length = 200
pvalue = 0.05
block_length = 5
multiple  = 1
homo = 0,1
savefile = {prefix}_correspondence_information.csv
""".format(prefix=prefix, lens1=lens1, lens2=lens2)
	cmd5 = "wgdi -c total.conf"
	return (conf5, cmd5)

def run_blockks(lens1, lens2, genome1_name, genome2_name, prefix):
	conf6 = """
[blockks]
lens1 = {lens1}
lens2 = {lens2}
genome1_name =  {genome1_name}
genome2_name =  {genome2_name}
blockinfo = {prefix}_block_information.csv
pvalue = 0.05
tandem = true
tandem_length = 200
markersize = 1
area = 0,2
block_length = 10
figsize = 8,8
savefig = {prefix}.ks.dotplot.pdf
""".format(lens1=lens1, lens2=lens2, genome1_name=genome1_name, genome2_name=genome2_name, prefix=prefix)
	cmd6 = "wgdi -bk total.conf"
	return (conf6, cmd6)

def run_kspeak(prefix):
	conf7 = """
[kspeaks]
 blockinfo = {prefix}_block_information.csv
 pvalue = 0.05
 tandem = true
 block_length = 10
 ks_area = 0,10
 multiple = 1
 homo = 0,1
 fontsize = 9
 area = 0,3
 figsize = 10,6.18
 savefig = {prefix}.kspeak.pdf
 savefile = {prefix}.kspeak.csv
""".format(prefix=prefix)
	cmd7 = "wgdi -kp total.conf"
	return (conf7, cmd7)

def run_peakfit(prefix):
	conf8="""
[peaksfit] 
blockinfo = {prefix}_block_information.csv
mode = median
tandem_length = 200
bins_number = 200
fontsize = 9 
area = -1,3
figsize = 10,6.18 
savefig = {prefix}.peakfit.pdf
""".format(prefix=prefix)
	cmd8 = "wgdi -pf total.conf"
	return (conf8, cmd8)


def run_km(blastp_result, prefix, gff1, gff2, aak_color, lens1):
	conf9=f"""
[karyotype_mapping]
blast = {blastp_result}
blast_reverse = false
gff1 = {gff1}
gff2 = {gff2}
score = 100
evalue = 1e-5
repeat_number = 5
ancestor_top = {aak_color}
the_other_lens = {lens1}
blockinfo = {prefix}_block_information.csv
blockinfo_reverse = false
limit_length = 5
the_other_ancestor_file =  km_result.txt
	"""
	cmd9 = "wgdi -km total.conf"
	return(conf9, cmd9)

def run_k():
	conf10=f"""
[karyotype]
ancestor = km_result.txt
width = 0.5
figsize = 10,6.18
savefig = karyotype.png
	"""
	cmd10 = "wgdi -k total.conf"
	return(conf10, cmd10)



def main(blastp_result, gff1, gff2, lens1, lens2, genome1_name, genome2_name, aak_color, gff1_cds, gff1_pep):
	blastp_result = os.path.abspath(blastp_result)
	gff1 = os.path.abspath(gff1)
	gff2 = os.path.abspath(gff2)
	lens1 = os.path.abspath(lens1)
	lens2 = os.path.abspath(lens2)
	cds_fa = os.path.abspath(gff1_cds)
	pep_fa = os.path.abspath(gff1_pep)

	prefix = "_".join([genome1_name, genome2_name])
	a1 = run_dotplot(blastp_result, gff1, gff2, lens1, lens2, genome1_name, genome2_name, prefix)
	a2 = run_collinearity(blastp_result, gff1, gff2, lens1, lens2, prefix)
	a3 = run_ks(gff1_cds, gff1_pep, prefix)
	
	a4 = run_blockinfo(blastp_result, gff1, gff2, lens1, lens2, prefix)
	a5 = run_correspondence(prefix, lens1, lens2)
	a6 = run_blockks(lens1, lens2, genome1_name, genome2_name, prefix)
	#a7 = run_kspeak(prefix)
	#a8 = run_peakfit(prefix)
	a9 = run_km(blastp_result, prefix, gff1, gff2, aak_color, lens1)
	a10 = run_k()



	with open("total.conf","w") as ouf:
		print("\n".join([a1[0], a2[0], a3[0], a4[0], a5[0], a6[0], a9[0], a10[0]]),file=ouf)

	print("\n".join([a1[1], a2[1], a3[1], a4[1], a5[1], a6[1], a9[1], a10[1]]))

if __name__ == '__main__':
	if len(sys.argv) != 11:
		print("python %s blastp_result gff1 gff2_anc lens1 lens2_anc genome1_name genome2_name aak_color all_cds all_pep"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9],sys.argv[10])


