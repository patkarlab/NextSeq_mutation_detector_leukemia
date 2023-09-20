#!/usr/bin/env python3

import sys
import csv
import re
import os
import openpyxl
import pandas as pd

acmg_output = sys.argv[1]	# output file from acmg (tsv)
csv_file = sys.argv[2] 		# cancervar input file
acmgcsv = sys.argv[3] 		# Output csv file with acmg mapping
#input_xlsx_file = sys.argv[4]	# Excel file to append the sheet

output_file = open (acmgcsv,'w')
# Extract the information from acmg output
acmg_class_dict = dict()
acmg_rules_dict = dict()
acmg_classify_dict = dict()
if os.path.getsize(acmg_output) != 0:
	with open (acmg_output, 'r') as tsv_file:
		tsv_handle = csv.reader(tsv_file, delimiter = '\t')
		header = next (tsv_handle)
		for tsv_values in tsv_handle:
			chrom = tsv_values[0]
			pos = tsv_values[1]
			ref = tsv_values[2]
			alt = tsv_values[3]
			acmg_class = tsv_values[4]
			acmg_rules = tsv_values[5]
			classificatn = tsv_values[6]

			if re.search('[a-zA-Z]', str(chrom)):
				chrom = re.sub ("chr","", chrom, flags = re.IGNORECASE)
				chrom = re.sub ("X","23", chrom, flags = re.IGNORECASE)
				chrom = re.sub ("Y","y", chrom, flags = re.IGNORECASE)
			

			acmg_class = re.sub (",",";", acmg_class, flags = re.IGNORECASE)
			acmg_rules = re.sub (",",";", acmg_rules, flags = re.IGNORECASE)
			classificatn = re.sub (",",";", classificatn, flags = re.IGNORECASE)

			map_id = str(chrom) + ':' +''.join(str(pos)) + ':' +''.join(str(ref)) + ':' +''.join(str(alt))
			acmg_class_dict[map_id] = acmg_class
			acmg_rules_dict[map_id] = acmg_rules
			acmg_classify_dict[map_id] = classificatn
			#print (acmg_class_dict[map_id], acmg_rules_dict[map_id], acmg_classify_dict[map_id], sep ="\t")
			#print (map_id)

# Map the acmg output to cancervar file
if os.path.getsize(csv_file) != 0:
	with open (csv_file,'r') as f:
		csv_handle = csv.reader(f)
		header = next (csv_handle)
		print (*header,"acmg classification", "acmg rules", "classifications", file=output_file, sep=",")
		for lines in csv_handle:
			chrom_csv = lines[0]
			pos_csv = lines[1]
			ref_csv = lines[3]
			alt_csv = lines[4]
			verdict = -1
			rules = -1
			classify = -1

			if re.search('[a-zA-Z]', str(chrom_csv)):
				chrom_csv = re.sub ("chr","", chrom_csv, flags = re.IGNORECASE)
				chrom_csv = re.sub ("X","23", chrom_csv, flags = re.IGNORECASE)
				chrom_csv = re.sub ("Y","y", chrom_csv, flags = re.IGNORECASE)
			
			csv_id = str(chrom_csv) + ':' +''.join(str(pos_csv)) + ':' +''.join(str(ref_csv)) + ':' +''.join(str(alt_csv))
			if csv_id in acmg_class_dict:
				verdict = acmg_class_dict[csv_id]
				rules = acmg_rules_dict[csv_id]
				classify = acmg_classify_dict[csv_id]

			print (*lines, verdict, rules, classify, file=output_file, sep=",")
			#print (csv_id)

# Add the csv file as a sheet to the final excel sheet
#wb = openpyxl.load_workbook(input_xlsx_file)
#if os.path.getsize(acmgcsv) != 0:
#	worksheet = wb.create_sheet('acmg')
#	read_csv = csv.reader(open(acmgcsv, 'r', encoding='utf-8'), delimiter = ',')
#	ind = 1 
#	for rows in read_csv:
#		worksheet.append(rows)
#		print (ind)
#		ind = ind + 1

#wb.save(input_xlsx_file)
#xl_file = pd.read_excel(input_xlsx_file, sheet_name=None)
#with pd.ExcelWriter('output_temp.xlsx') as writer:
#	for sheet_names in xl_file.keys():
#		df = pd.read_excel(input_xlsx_file, sheet_name=sheet_names)
#		df.fillna(value= -1 , inplace=True)
#		
#		if 'VAF' in df:
#			df["VAF"] = df["VAF"].astype(str).str.rstrip("%").astype(float)
#		df = df.rename(columns={"VAF":"VAF (%)"})
#
#		if 'PopFreqMax' in df:
#			df["PopFreqMax"] = df["PopFreqMax"].astype(str).astype(float)
#			if 'acmg' in sheet_names:
#				df = df[df.PopFreqMax < 0.01]           # Removing varinats with Population frequency less than 0.01
#
#		if 'ALT_COUNT' in df:
#			df["ALT_COUNT"] = df["ALT_COUNT"].astype(str).astype(float)
#
#		if 'VariantCaller_Count' in df and 'acmg' in sheet_names:
#			df = df[df.VariantCaller_Count > 1]		# Removing variants which are called by only one variant caller
#
#		df.to_excel(writer, sheet_name=sheet_names, index=False)
