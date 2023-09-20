#!/usr/bin/env python3
# This script will check the final.concat.csv files and extract the correponding info from pcgr and cpsr outputs
# The output of this code 
import sys
import re
import csv

sample = sys.argv[1]
filepath = sys.argv[2]
pcgr_tsv_file = sys.argv[3]
cpsr_tsv_file = sys.argv[4]
csv_file = filepath + sample + '.final.concat.csv'
output = sample + '.final.concat_append.csv'

pcgr_map1 = dict()                    # Dictionary for markers
pcgr_map2 = dict()             # Dictionary for the marker occurence
pcgr_map3 = dict()   		# Dictionary for comments
pcgr_map4 = dict()
pcgr_map5 = dict()
pcgr_map6 = dict()
pcgr_map7 = dict()
with open (pcgr_tsv_file,'r') as tsv:
	tsv_handle = csv.reader(tsv, delimiter="\t")
	header = next (tsv_handle)
	column1 = header.index('GENOMIC_CHANGE')
	column2 = header.index('HGVSp')
	column3 = header.index('HGVSc')
	column4 = header.index('DBSNPRSID')
	column5 = header.index('CLINVAR')
	column6 = header.index('CLINVAR_CLNSIG')
	column7 = header.index('TIER')

	for pcgr_lines in tsv_handle:
		pcgr_chr = pcgr_lines[0]
		pcgr_pos = pcgr_lines[1]
		pcgr_ref = pcgr_lines[2]
		pcgr_alt = pcgr_lines[3]

		genomic_change = str(pcgr_lines[column1]).replace(",", "&")
		hgvsp = str(pcgr_lines[column2]).replace(",", "&")
		hgvsc = str(pcgr_lines[column3]).replace(",", "&")
		dbsnprsid = str(pcgr_lines[column4]).replace(",", "&")
		clinvar = str(pcgr_lines[column5]).replace(",", "&")
		clinvar_clnsig = str(pcgr_lines[column6]).replace(",", "&")
		tier = str(pcgr_lines[column7]).replace(",", "&")
		
		if re.search('[a-zA-Z]', str(pcgr_chr)):
			pcgr_chr = re.sub ("chr","", pcgr_chr, flags = re.IGNORECASE)
			pcgr_chr = re.sub ("X","23", pcgr_chr, flags = re.IGNORECASE)
			pcgr_chr = re.sub ("Y","y", pcgr_chr, flags = re.IGNORECASE)

		pcgr_id = str(pcgr_chr) + ':' +''.join(pcgr_pos) + ':' +''.join(pcgr_ref) + ':' +''.join(pcgr_alt)
		#print (pcgr_id)
		pcgr_map1[pcgr_id] = genomic_change
		pcgr_map2[pcgr_id] = hgvsp
		pcgr_map3[pcgr_id] = hgvsc
		pcgr_map4[pcgr_id] = dbsnprsid
		pcgr_map5[pcgr_id] = clinvar
		pcgr_map6[pcgr_id] = clinvar_clnsig
		pcgr_map7[pcgr_id] = tier


cpsr_map1 = dict ()
cpsr_map2 = dict ()
cpsr_map3 = dict ()
cpsr_map4 = dict ()
with open (cpsr_tsv_file,'r') as tsv1:
	tsv_handle1 = csv.reader(tsv1, delimiter="\t")
	header1 = next (tsv_handle1)
	col1 = header1.index('VAR_ID')
	col2 = header1.index('CLINVAR_CLASSIFICATION')
	col3 = header1.index('CLINVAR_VARIANT_ORIGIN')
	col4 = header1.index('FINAL_CLASSIFICATION')
	col5 = header1.index('CPSR_CLASSIFICATION')

	for cpsr_lines in tsv_handle1:
		chr_pos_ref_alt = cpsr_lines[col1].split('_')
		cpsr_chr = chr_pos_ref_alt[0]
		cpsr_pos = chr_pos_ref_alt[1]
		cpsr_ref = chr_pos_ref_alt[2]
		cpsr_alt = chr_pos_ref_alt[3]
		
		clinvar_class = str(cpsr_lines[col2]).replace(",", "&")
		clinvar_varor = str(cpsr_lines[col3]).replace(",", "&")
		final_classif = str(cpsr_lines[col4]).replace(",", "&")
		cpsr_classifi = str(cpsr_lines[col5]).replace(",", "&")

		if re.search('[a-zA-Z]', str(cpsr_chr)):
			cpsr_chr = re.sub ("chr","", cpsr_chr, flags = re.IGNORECASE)
			cpsr_chr = re.sub ("X","23", cpsr_chr, flags = re.IGNORECASE)
			cpsr_chr = re.sub ("Y","y", cpsr_chr, flags = re.IGNORECASE)

		cpsr_id = str(cpsr_chr) + ':' +''.join(cpsr_pos) + ':' +''.join(cpsr_ref) + ':' +''.join(cpsr_alt)
		cpsr_map1[cpsr_id] = clinvar_class
		cpsr_map2[cpsr_id] = clinvar_varor
		cpsr_map3[cpsr_id] = final_classif
		cpsr_map4[cpsr_id] = cpsr_classifi


outfile = open (output,'w')
with open (csv_file,'r') as vcf:
	csv_handle = csv.reader(vcf)
	header2 = next (csv_handle)      # Removing the header
	print (*header2,"GENOMIC_CHANGE","HGVSp","HGVSc","DBSNPRSID","CLINVAR","CLINVAR_CLNSIG","PCGR_TIER","CLINVAR_CLASSIFICATION","CLINVAR_VARIANT_ORIGIN","	CLINVAR+CPSR_CLASSIFICATION","CPSR_CLASSIFICATION",file=outfile, sep=",")
	val1 = val2 = val3 = val4 = val5 = val6 = val7 = -1
	val8 = val9 = val10 = val11 = -1

	for str_lines in csv_handle:
		csv_chr = str_lines[0]
		csv_pos = str_lines[1]
		csv_ref = str_lines[3]
		csv_alt = str_lines[4]
		
		if re.search('[a-zA-Z]', str(csv_chr)):
			csv_chr = re.sub ("chr","", csv_chr, flags = re.IGNORECASE)
			csv_chr = re.sub ("X","23", csv_chr, flags = re.IGNORECASE)
			csv_chr = re.sub ("Y","y", csv_chr, flags = re.IGNORECASE)

		csv_id = str(csv_chr) + ':' +''.join(csv_pos) + ':' +''.join(csv_ref) + ':' +''.join(csv_alt)
		#print (csv_id)
		if csv_id in pcgr_map1:
			val1 = pcgr_map1[csv_id]

		if csv_id in pcgr_map2:
			val2 = pcgr_map2[csv_id]

		if csv_id in pcgr_map3:
			val3 = pcgr_map3[csv_id]

		if csv_id in pcgr_map4:
			val4 = pcgr_map4[csv_id]

		if csv_id in pcgr_map5:
			val5 = pcgr_map5[csv_id]

		if csv_id in pcgr_map6:
			val6 = pcgr_map6[csv_id]

		if csv_id in pcgr_map7:
			val7 = pcgr_map7[csv_id]

		if csv_id in cpsr_map1:
			val8 = cpsr_map1[csv_id]

		if csv_id in cpsr_map2:
			val9 = cpsr_map2[csv_id]

		if csv_id in cpsr_map3:
			val10 = cpsr_map3[csv_id]

		if csv_id in cpsr_map4:
			val11 = cpsr_map4[csv_id]

		print (*str_lines, val1 , val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, file=outfile, sep=",")


outfile.close()
