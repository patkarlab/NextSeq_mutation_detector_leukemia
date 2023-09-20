#!/usr/bin/env python3
# This script will take the *cancervar.csv as input and will write the variants.txt 

import sys
import csv
import os
import re

infile = sys.argv[1]	# Input csv file
outfile = sys.argv[2]	# Output file with variants "variants.txt"
# Columns to screen
PopFreqMax = "PopFreqMax"
popfreqlimit = 0.01
FuncrefGene = "Func.refGene"
ExonicFuncrefGene = "ExonicFunc.refGene"
AltCount = "ALT_COUNT"
altcountlimit = 20
Vaf = "VAF"
Vaflimit = 5
# Strings to search
exonic_splice = "exonic|splic"
synonymous = r"\bsynonymous"

output_file = open (outfile,'w')
if os.path.getsize(infile) != 0:
	with open (infile,'r') as f:
		csv_handle = csv.reader(f)
		header = next (csv_handle)	
		for lines in csv_handle:
			popfreq = float(lines[header.index(PopFreqMax)])
			funcrefgene = str(lines[header.index(FuncrefGene)])
			exonicfuncrefgene = str(lines[header.index(ExonicFuncrefGene)])
			alt_count = int(lines[header.index(AltCount)])
			vaf = float(lines[header.index(Vaf)].rstrip("%"))

			chrm = lines[0]
			pos = lines[1]
			ref = lines[3]
			alt = lines[4]
			if ref == '-':
				ref=""
			if alt == '-':
				alt=""

			if popfreq < popfreqlimit:														# Remove variants with freq < popfreqlimit
				if re.search(exonic_splice, funcrefgene, flags = re.IGNORECASE):			# Select variants which have 'exonic_splice' enteries
#					if chrm == "chr17":
#						print (chrm, pos, ref, alt, file=output_file, sep='-')
					if not re.search(synonymous, exonicfuncrefgene, flags = re.IGNORECASE):	# Remove variants other than 'synonymous'
						if alt_count >= altcountlimit and vaf >= Vaflimit:
							print (chrm, pos, ref, alt, file=output_file, sep='-')

output_file.close ()
