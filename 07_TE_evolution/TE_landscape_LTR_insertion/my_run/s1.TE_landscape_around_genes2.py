# 蛋白编码基因前后10kp附近，各个类型TE的分布情况
import os
import sys
import re
from Bio import SeqIO
from BCBio import GFF

# def filter_gff(gff_file, output_file):
# 	with open(gff_file, 'r') as fin, open(output_file, 'w') as fout:
# 		for line in fin:
# 			if not line.startswith('#'):
# 				fout.write(line)

def get_upstream_downstream(record, distance=10000): 
	if record.strand == 1:
		upstream = (record.location.start - distance, record.location.start)
		downstream = (record.location.end, record.location.end + distance)
	else:
		upstream = (record.location.end, record.location.end + distance)
		downstream = (record.location.start - distance, record.location.start)
	return upstream, downstream

def overlap(region, feature):
	feature_start = min(feature.location.start, feature.location.end)
	feature_end = max(feature.location.start, feature.location.end)
	return feature_start < region[1] and feature_end > region[0]

def count_TE_types(gene_annotation_file, TE_annotation_file, genome_file):
	genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))  
	gene_annotation_data = GFF.parse(gene_annotation_file, base_dict=genome)
	TE_annotation_data = GFF.parse(TE_annotation_file, base_dict=genome)

	counts = {}
	i = 1
	for record in gene_annotation_data:
		for feature in record.features:
			for sub_feature in feature.sub_features:
				if sub_feature.type == "mRNA":
					i += 1
					upstream, downstream = get_upstream_downstream(sub_feature)
				
					for TE_record in TE_annotation_data:
						if TE_record.id == record.id:
							for TE_sub_feature in TE_record.features:
								if overlap(upstream, TE_sub_feature) or overlap(downstream, TE_sub_feature):
									if overlap(upstream, TE_sub_feature):
										dis = TE_sub_feature.location.end - sub_feature.location.start
										counts.setdefault(TE_sub_feature.type,[]).append(dis)
									if overlap(downstream, TE_sub_feature):
										dis = TE_sub_feature.location.start - sub_feature.location.end
										counts.setdefault(TE_sub_feature.type,[]).append(dis)

								# TE_info = {"strand":TE_sub_feature.strand, 
								# 			"start":TE_sub_feature.location.start, 
								# 			"end":TE_sub_feature.location.end, 
								# 			"region": "upstream" if overlap(upstream, TE_sub_feature) else "downstream"}
								# counts[TE_sub_feature.type].setdefault("TE_info", []).append(TE_info)
								# counts[TE_sub_feature.type].setdefault("counts", 0)
								# counts[TE_sub_feature.type]["counts"] += 1
	print(i)
	return counts

def main(gene_annotation_file, TE_annotation_file, genome_file):
	res = count_TE_types(gene_annotation_file, TE_annotation_file, genome_file)
	for k,v in res.items():
		with open(f"TE_{k}.txt","w") as ouf:
			for i in v:
				print(i/1000,file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s gene_annotation_file TE_annotation_file genome_file"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3])
