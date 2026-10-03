

import re
import os
import sys
import pandas as pd
from plotnine import *
import skmisc
import numpy as np

def main(simu_res):
	data = pd.read_csv(simu_res,sep=" ",names=["genome_num","pan_fam","core_fam"])
	data = data.dropna()
	p1 = ggplot(data) + \
		geom_point(aes(x="genome_num",y="pan_fam",),color="#2C8EC4",alpha=0.7,shape=6,size=2.5) + \
		geom_smooth(aes(x="genome_num",y="pan_fam"),method="loess",color="#3F80C2",se=False,alpha=0.3,size=0.8) + \
		geom_point(aes("genome_num","core_fam"),color="#EC645F",alpha=0.7,shape="o",size=2.5) + \
		geom_smooth(aes(x="genome_num",y="core_fam"),method="loess",color="#C25E63",se=False,alpha=0.3,size=0.8) + \
		theme_classic()

	# p1.save("p1.simu_pan.",dpi=100,format='svg',width=10,height=10)
	p1.save("p1.simu_pan.pdf",dpi=100,format='pdf',width=6,height=6)
	
	#ggsave(filename="p1.simu_pan.png",plot=p1)
	#save_as_pdf_pages(p1,"p1.simu_pan.png")

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s simulate_pan_core_file"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])