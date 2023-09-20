#!/usr/bin/env python3

import sys
import csv
import re
import os
import openpyxl
import pandas as pd

acmgcsv = sys.argv[1] 		# Output csv file with acmg mapping
input_xlsx_file = sys.argv[2]	# Excel file to append the sheet

# Add the csv file as a sheet to the final excel sheet
wb = openpyxl.load_workbook(input_xlsx_file)
if os.path.getsize(acmgcsv) != 0:
	worksheet = wb.create_sheet('acmg')
	read_csv = csv.reader(open(acmgcsv, 'r', encoding='utf-8'), delimiter = ',')
	for rows in read_csv:
		worksheet.append(rows)

wb.save(input_xlsx_file)

xl_file = pd.read_excel(input_xlsx_file, sheet_name=None)
with pd.ExcelWriter('output_temp.xlsx') as writer:
	for sheet_names in xl_file.keys():
		df = pd.read_excel(input_xlsx_file, sheet_name=sheet_names)
		df.fillna(value= -1 , inplace=True)
		
		if 'VAF' in df:
			df["VAF"] = df["VAF"].astype(str).str.rstrip("%").astype(float)
		df = df.rename(columns={"VAF":"VAF (%)"})

		if 'PopFreqMax' in df:
			df["PopFreqMax"] = df["PopFreqMax"].astype(str).astype(float)
			if 'acmg' in sheet_names:
				df = df[df.PopFreqMax < 0.01]           # Removing varinats with Population frequency less than 0.01

		if 'ALT_COUNT' in df:
			df["ALT_COUNT"] = df["ALT_COUNT"].astype(str).astype(float)

		if 'VariantCaller_Count' in df and 'acmg' in sheet_names:
			df = df[df.VariantCaller_Count > 1]		# Removing variants which are called by only one variant caller

		df.to_excel(writer, sheet_name=sheet_names, index=False)
