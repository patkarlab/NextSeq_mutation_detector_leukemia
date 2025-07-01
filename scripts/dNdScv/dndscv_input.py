#!/usr/bin/env python3

import sys
import os
import pandas as pd

args = sys.argv
input_excel = args[1]
sample_id = args[2]
sheet_name = "append_final_concat"
output_name = sample_id + ''.join('_dndscv.tsv')
output_gene_list = sample_id + ''.join('_genelist.tsv')
genes = []

if os.path.getsize(input_excel) != 0:
	df = pd.read_excel(input_excel, sheet_name=sheet_name)
	dndscv_df = df [["Chr", "Start", "Ref", "Alt"]]
	dndscv_df.insert(0, 'sampleID', sample_id)	# Inserting a column with sample name at the start
	dndscv_df_rename = dndscv_df.rename(columns={'Chr': 'chr', 'Start' : 'pos', 'Ref' : 'ref' , 'Alt' : 'alt' })	# Renaming the columns
	dndscv_df_rename['chr'] = dndscv_df_rename['chr'].str.replace('chr', '')
	# Remove ';' from gene names
	for values in df ["Gene.refGene"]:
		gene_list = values.split(';')
		for name in gene_list:
			genes.append(name)
	gene_names = pd.DataFrame(genes).drop_duplicates()
else:
	dndscv_df_rename = pd.DataFrame()
	gene_names = pd.DataFrame()

dndscv_df_rename.to_csv (output_name, sep='\t', index=False)
gene_names.to_csv (output_gene_list, sep='\t', header=False, index=False)