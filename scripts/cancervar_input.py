#!/usr/bin/env python3

import sys
import csv
import re
import pandas as pd

csv_file = sys.argv[1]	# csv file from the annovar output
output_file = sys.argv[2] # tsv file for filtered output
outfile = open (output_file,'w')

df = pd.read_csv(csv_file)
for row_nums in df.index:
	variant_sum = 0
	#print (df['Start'][row_nums], df['End'][row_nums], df['Ref'][row_nums], df['Alt'][row_nums], df['Otherinfo1'][row_nums], file=output_file, sep="\t")
	if df['PopFreqMax'][row_nums] < 0.01:
		other_info_column_list = df['Otherinfo1'][row_nums].split(';')
		for items in other_info_column_list:
			value_list = items.split()
			for values in value_list:
				if 'MDLKFPGS' in values:
					variants = values.split('=')[1].replace(',','')
					for callers in variants:
						variant_sum = variant_sum + int(callers)
					if variant_sum > 1:
						print (df['Chr'][row_nums],df['Start'][row_nums], df['End'][row_nums], df['Ref'][row_nums], df['Alt'][row_nums], file=outfile,sep="\t")

outfile.close()
