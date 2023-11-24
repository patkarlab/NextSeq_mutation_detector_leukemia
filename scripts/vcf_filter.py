#!/usr/bin/env python3
#This script will filter the vcf for exonic, synonymous and mutation variants

import sys, os, re

invcf = sys.argv[1]		# Input vcf file
outfile = sys.argv[2]	# Ouput	vcf file

info_column = "INFO"
FuncrefGene = r'^Func.refGene'
exonic = "exonic"

ExonicFuncrefGene = r'^ExonicFunc.refGene'
synonymous = r"\bsynonymous"

PopFreqMax = "PopFreqMax"
popfreqlimit = 0.01

output_file = open(outfile,'w')
if os.path.getsize(invcf) != 0:
	with open (invcf,'r') as vcf:
		for lines in vcf:
			lines = lines.rstrip()
			if re.search(r'^#', lines):	# Printing the entire header of the input vcf file
				print (lines, file=output_file)

				if re.search(r'^#CHROM', lines):	# Extracting the column no. for info column 
					info_column_index = lines.split('\t').index(info_column)
					#print (info_column_index)
			else:
				info_column_data = lines.split('\t')[info_column_index]
				info_column_list = info_column_data.split(';')

				exonic_var = 0
				non_synonymous_var = 0
				maf_filter = 0

				for values in info_column_list:
					if re.search(FuncrefGene, values):
						if exonic in values.split('=')[1]:		# Extracting only the exonic variants
							exonic_var = 1

					if re.search(ExonicFuncrefGene, values):	# Removing the synonymous
						if not re.search(synonymous, values.split('=')[1], flags = re.IGNORECASE):
							non_synonymous_var = 1
					
					if re.search(PopFreqMax, values, flags = re.IGNORECASE):	# Screening for PopFreq
						if values.split('=')[1] == '.':
							pop_freq = -1
						else:
							pop_freq = float(values.split('=')[1])

						if pop_freq < popfreqlimit:
							maf_filter = 1

				if exonic_var == 1 and non_synonymous_var == 1 and maf_filter == 1:
					print (lines, file=output_file)	
else:				
	pass

output_file.close()
