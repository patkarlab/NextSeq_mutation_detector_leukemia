#!/bin/env python3

import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(description="reformat flt3_itd_ext vcf to tsv ")
parser.add_argument("-v", "--flt3_itd_ext_vcf", required=True, help="Path to flt3_itd_ext vcf")
parser.add_argument("-o", "--output", required=True, help="Path to out tsv file")
args = parser.parse_args()

#read the vcf file
df = pd.read_csv(args.flt3_itd_ext_vcf, sep="\t",  comment = "#", names = ['CHROM','POS','ID',	'REF', 'ALT','QUAL','FILTER','INFO','FORMAT','Sample'])

if df.empty:
	df = pd.DataFrame(columns=['POS', 'REF', 'ALT', 'Strand', 'length','VAF%', 'CDS', 'AA',
	   'Raw Allelic Ratio(RAR)', 'Raw Total depth(RDP)',
	   'Raw Variant Depth(RVD)', 'Coverage Radius(CR)'])

	df.to_csv(args.output, sep= "\t", index=False)

else:	

	#split 'INFO' column and add the split columns to dataframe
	#split_info = df['INFO'].str.split(';', expand=True)
	#df = pd.concat([df, split_info], axis=1)

	split_info = df['INFO'].apply(lambda x: re.split(r';(?=\w+=)', x))
	max_len = split_info.map(len).max()
	split_info = pd.DataFrame(split_info.tolist(), columns=range(max_len))

	# Concatenate split columns back to df
	df = pd.concat([df, split_info], axis=1)

	#drop unnecessary columns
	columns_to_keep = [ 'POS','REF','ALT','Sample',1,2,3,4,11,12,13,14]
	df = df[columns_to_keep]

	#split selected columns from  previously split columns on '='
	cols_to_split = [1,2,3,4,11,12,13,14]
	delimiter = '='

	for col in cols_to_split:
		# Create new columns with names as '_1' and '_2' added to original name
		df[[f'{col}_1', f'{col}_2']] = df[col].str.split(delimiter, expand=True)

	#define the new column names for renaming selected columns
	column_mapping = {
		'1_2':'Strand',
		'2_2': 'length',
		'3_2': 'CDS',
		'4_2':'AA',
		'11_2':'Raw Allelic Ratio(RAR)',
		'12_2':'VAF',
		'13_2':'Raw Total depth(RDP)',
		'14_2':'Raw Variant Depth(RVD)'
	}

	#rename the columns
	df.rename(columns=column_mapping, inplace=True)

	columns_to_drop = [1,2,3,4,11,12,13,14,'1_1','2_1','3_1','4_1','11_1','12_1','13_1','14_1' ]
	df = df.drop(columns=columns_to_drop)

	#split the 'Sample' column
	split_sample = df['Sample'].str.split(':', expand=True)
	df = pd.concat([df, split_sample], axis=1)

	#multiply vaf with 100 to get percentage
	df['VAF'] = pd.to_numeric(df['VAF'])
	df['VAF%'] = df['VAF'] * 100

	#drop unnecessary columns
	columns_to_drop = ['Sample', 0,1,2,3,4,5,6,7,8,9,'VAF' ]
	df = df.drop(columns=columns_to_drop)

	#rename column
	df = df.rename(columns={10:'Coverage Radius(CR)'})

	#rearrange columns
	new_order = [ 'POS', 'REF', 'ALT', 'Strand', 'length','VAF%', 'CDS', 'AA',
	   'Raw Allelic Ratio(RAR)', 'Raw Total depth(RDP)',
	   'Raw Variant Depth(RVD)', 'Coverage Radius(CR)' ]
	df = df[new_order]

	df.to_csv(args.output, index=False)

