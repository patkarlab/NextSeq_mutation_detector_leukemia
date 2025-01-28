#!/usr/bin/env python
# This script will write a chrwise list file given a bed file
import sys

in_bedfile = sys.argv[1]	# input bedfile


def region_check(name):
	"""
	This function will modify the region name 
	"""
	mod_name = ""
	split_name = name.split('_')
	if split_name[0].isnumeric and split_name[1].isnumeric():
		mod_name = '_'.join(split_name[2:])
	else:
		mod_name = name
	return mod_name

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
		region_name_list.append(region_check(region_data[3]))	# Removing the numeric values at the start of region name
		if (len(region_list) == 100):	# print 100 values per line
			print ((',').join(region_list),"\t",(',').join(region_name_list))
			region_list.clear()
			region_name_list.clear()
print ((',').join(region_list),"\t",(',').join(region_name_list))