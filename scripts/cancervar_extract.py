#! /usr/bin/env python3

import sys
import csv
import re

cancervar_file = sys.argv[1]# pred file from the CancerVar
vep_file = sys.argv[2]		# file from VEP
output_file = sys.argv[3]	# output file

map_evidence = dict()		# Dictionary for the cancervar evidence
map_ensembl_score = dict()	# Dictionary for score
with open (cancervar_file, 'r') as tsv_file:
	cancervar_handle = csv.reader(tsv_file, delimiter = '\t')
	header = next (cancervar_handle)
	for vep_values in cancervar_handle:
		chrom = vep_values[0]
		start = vep_values[1]
		end = vep_values[2]
		ref = vep_values[3]
		alt = vep_values[4]
		evidence = str(vep_values[13]).replace(',',' ')
		ensembl_score = str(vep_values[44]).replace(',',' ')

		if re.search('[a-zA-Z]', str(chrom)):
			chrom = re.sub ("chr","", chrom, flags = re.IGNORECASE)
			chrom = re.sub ("X","23", chrom, flags = re.IGNORECASE)
			chrom = re.sub ("Y","y", chrom, flags = re.IGNORECASE)

		evidence = re.sub (r'EVS=.*',"", evidence, flags = re.IGNORECASE)
		evidence = re.sub ("CancerVar","CancerVar Score", evidence, flags = re.IGNORECASE)
		evidence = re.sub ("#"," ", evidence, flags = re.IGNORECASE)
		cancervar_id = str(chrom) + ':' +''.join(str(start)) + ':' +''.join(str(end)) + ':' +''.join(str(ref)) + ':' +''.join(str(alt))
		map_evidence[cancervar_id] = evidence
		map_ensembl_score[cancervar_id] = ensembl_score
	
outfile = open (output_file, 'w')
with open (vep_file,'r') as vcf:
	csv_handle = csv.reader(vcf, delimiter = ',')
	header = next (csv_handle)      # Removing the header
	print (*header,"CancerVar evidence", "AI score", file=outfile, sep=",")
	for str_lines in csv_handle:
		chrom_vep = str_lines[0]
		start_vep = str_lines[1]
		end_vep = str_lines[2]
		ref_vep = str_lines[3]
		alt_vep = str_lines[4]
		evidence_val = -1
		ensembl_score_val = -1		
		
		if re.search('[a-zA-Z]', str(chrom_vep)):
			chrom_vep = re.sub ("chr","", chrom_vep, flags = re.IGNORECASE)
			chrom_vep = re.sub ("X","23", chrom_vep, flags = re.IGNORECASE)
			chrom_vep = re.sub ("Y","y", chrom_vep, flags = re.IGNORECASE)
		
		vep_id = str(chrom_vep) + ':' +''.join(str(start_vep)) + ':' +''.join(str(end_vep)) + ':' +''.join(str(ref_vep)) + ':' +''.join(str(alt_vep))
		if vep_id in map_evidence:
			evidence_val = map_evidence[vep_id]
			ensembl_score_val = map_ensembl_score[vep_id]
			print (*str_lines, evidence_val, ensembl_score_val, file=outfile, sep=",")
		else:
			print (*str_lines, evidence_val, ensembl_score_val, file=outfile, sep=",")

outfile.close ()
