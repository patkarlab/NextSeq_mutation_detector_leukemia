#!/usr/bin/env python
# This script will write a chrwise list file given a bed file
import sys

in_bedfile = sys.argv[1]	# input bedfile

chr_id=""
region_list = []
region_name_list = []
with open(in_bedfile) as bedfile:
	for lines in bedfile:
		region_data = lines.strip().split('\t')
		chr = region_data[0]
		region = chr + ':' + str(region_data[1]) + '-' + str(region_data[2])
		if chr != chr_id:
			chr_id = chr
			print ((',').join(region_list),"\t",(',').join(region_name_list))
			region_list.clear()
			region_name_list.clear()
		else:
			pass
		region_list.append(region)
		region_name_list.append(region_data[3])
print ((',').join(region_list),"\t",(',').join(region_name_list))