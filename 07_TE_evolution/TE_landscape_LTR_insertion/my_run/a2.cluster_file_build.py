#获取文件两两组合
import itertools
import os
import sys

def main(indir):
	# path = '/usr_storage/stt/Populus_pangenome/TE_analysis/fl_LTR/'
	

	#combines = list(itertools.combinations('Psze_LTR.fa','Pkor_LTR.fa','Ptri_LTR.fa','Pdel_LTR.fa','Psim_LTR.fa',
	#                                       'Pwua_LTR.fa','Pyun_LTR.fa','Ppse_LTR.fa','Plas_LTR.fa','Pdav_LTR.fa','Prot_LTR.fa',
	#                                       'Ptre_LTR.fa','Palby_LTR.fa','Palb_LTR.fa','Pade_LTR.fa','Pqio_LTR.fa','Ppru_LTR.fa',
	#                                       'Pili_LTR.fa','Peup_LTR.fa'))

	files = os.listdir(indir) 
	combines = list(itertools.combinations(files,2))

	for combine in combines:
		s1 = combine[0].split(r'.')[0]
		s2 = combine[1].split(r'.')[0]
		s = s1+'_'+s2
		mkdir_cmd = f'mkdir {s}'
		cat_cmd = f"cat {combine[0]} {combine[1]} > {s}/ltr_juction.fa"

		# cat_cmd = 'cat /usr_storage/stt/Populus_pangenome/TE_analysis/fl_LTR/'+combine[0]+' /usr_storage/stt/Populus_pangenome/TE_analysis/fl_LTR/'+combine[1]+' > '+s+'/ltr_juction.fa'
		# print(mkdir_cmd)
		# print(cat_cmd)
		os.system(mkdir_cmd)
		os.system(cat_cmd)

main(sys.argv[1])
 
