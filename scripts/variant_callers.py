#!/usr/bin/env python3
# This script will take the somaticseq combined .vcf and the output .tiers.tsv file 

import sys
import re
import pandas as pd
import csv

tiers_csv_file = sys.argv[1]	# 
combined_vcf_file = sys.argv[2]	# 
outfile = sys.argv[3]			# Output file

output_file = open (outfile,'w')
#print ("The sequence of variant callers is obtained from the combined vcf file")
with open (combined_vcf_file) as vcf:
	for vcf_lines in vcf:
		if 'ID=MVDLKFP' in vcf_lines:
			info_tag_split = vcf_lines.split(':')[1]
			variant_caller_ids = info_tag_split.split(',')
			filtered_ids = [re.sub("[^A-Za-z0-9]", "", i) for i in variant_caller_ids]	# Removing non alphabetic/numeric characters
			print (filtered_ids)
		else:
			pass

with open (tiers_csv_file,'r') as tsv:
	tsv_handle = csv.reader(tsv, delimiter="\t")
	header = next (tsv_handle)
	print (*header, "VARIANT_CALLERS", "NO_OF_VARIANTS", file=output_file, sep="\t")
	variant_caller_column = header.index('CALL_CONFIDENCE')	# Index of the required column
	for lines in tsv_handle:
		variant_list = list()
		variant_callers = list(lines[variant_caller_column].split(','))
		variant_callers = list(map(int, variant_callers))
		no_of_variant_callers =  sum(variant_callers)
		for j in range(len(variant_callers)):
			if variant_callers[j] == 1:
				variant_list.append(filtered_ids[j])
				
		print (*lines,variant_list,no_of_variant_callers, file=output_file, sep="\t")
