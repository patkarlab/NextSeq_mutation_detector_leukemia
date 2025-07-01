#!/usr/bin/env python
#This script will take gene info from *_varout.tsv and extract the required dnds values from *_geneout.tsv 

import os, sys
import pandas as pd
import openpyxl

args = sys.argv
variants = args[1]
dnds_info = args[2]
input_excel = args[3]
sheet_name = "append_final_concat"
mutation_acronyms = {'Missense' : 'wmis_cv', 'Essential_Splice' : 'wspl_cv' , 'no-SNV' : 'wind_cv', 'Nonsense' : 'wnon_cv' }
remove_synon = "Synonymous"
default_val = int(-1)
final_col_name = "dNdScv"

if os.path.getsize(variants) != 0:
	df = pd.read_csv(variants, sep="\t")
	# print (df[["chr", "pos", "ref" , "mut", "gene", "impact"]])
	gene_list = df["gene"].drop_duplicates().to_list()
	# print (gene_list)
	if os.path.getsize(dnds_info) != 0:
		df_annot = pd.read_csv(dnds_info, sep="\t")
		# mapping gene names and dndscv values for variants
		merged_df = pd.merge(df, df_annot, left_on ="gene", right_on = "gene_name" ,how='left') 

		# impact_column_index = merged_df.columns.get_loc('impact')
		# for index, row in merged_df.iterrows():
		# 	print(index, merged_df.iloc[index,impact_column_index])
		# print (merged_df.shape)
		
		# print (df_filtered[["chr", "chr", "ref" , "mut", "gene", "impact", "wmis_cv", "wnon_cv", "wspl_cv", "wind_cv"]])
		mutatns_mapped_df = merged_df[["chr", "pos", "ref" , "mut"]].copy() # Making a copy is required
		# Adding a new column for output and obtaining its index
		mutatns_mapped_df.insert(len(mutatns_mapped_df.columns), final_col_name, default_val)
		dNdScv_col_index = mutatns_mapped_df.columns.get_loc(final_col_name)
		chr_column_index = mutatns_mapped_df.columns.get_loc('chr')

		impact_column_index = merged_df.columns.get_loc('impact')
		# wnon_column_index = merged_df.columns.get_loc('wnon_cv')
		# wspl_column_index = merged_df.columns.get_loc('wspl_cv')
		# wind_column_index = merged_df.columns.get_loc('wind_cv')
		for index, row in mutatns_mapped_df.iterrows():
			impact_col_val = merged_df.iloc[index, impact_column_index]
			chr_col_val = "chr" + mutatns_mapped_df.iloc[index, chr_column_index]
			if impact_col_val in mutation_acronyms:
				dnds_type = mutation_acronyms[impact_col_val]
				# Obtain index for the dnds types
				windex = merged_df.columns.get_loc(dnds_type)
				dnds_val = merged_df.iloc[index, windex]
			else:
				dnds_val = default_val
			mutatns_mapped_df.iloc[index, dNdScv_col_index] = dnds_val
			mutatns_mapped_df.iloc[index, chr_column_index] = chr_col_val
		# print (merged_df.shape, mutatns_mapped_df.shape)
	else:
		mutatns_mapped_df = pd.DataFrame()
else:
	mutatns_mapped_df = pd.DataFrame()

# print(mutatns_mapped_df)
if os.path.getsize(input_excel) != 0:
	df_excel = pd.read_excel(input_excel, sheet_name=sheet_name)
	excel_mutatns_mergedf = pd.merge(df_excel, mutatns_mapped_df, left_on = ["Chr", "Start", "Ref", "Alt"] , right_on = ["chr", "pos", "ref", "mut"] , how='left')
	excel_mutatns_mergedf = excel_mutatns_mergedf.drop(columns=['chr', 'pos', 'ref', 'mut'])

	# Writing the output to a new sheet
	with pd.ExcelWriter(input_excel, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
		excel_mutatns_mergedf.to_excel(writer, sheet_name=final_col_name, index=False)

else:
	print ("Excel not modified")
