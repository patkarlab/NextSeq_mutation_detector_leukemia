#!/usr/bin/env python3
##############
# The analysis of json proceeds as follows
# Step 1 : Open it and load , use python library json
# Step 2 : check the type and length of individual objects
# Step 3 : If the object is a list => loop over it, If its a dictionary => print out the needed key's value 
# Iterate steps 2 and 3
##############
import json
import sys
import os

infile = sys.argv[1]	# Input is file in json format
outfile = sys.argv[2]	# Output is a tsv file

output = open (outfile,'w')
print ("chrom", "pos", "ref", "alt", "acmg classification", "acmg rules", "classifications", file=output, sep="\t")

if os.path.getsize(infile) != 0:
	with open (infile) as f:
		data = json.load(f)
else:
	sys.exit()		
	
def extract_data(values):
	chrm_pos_ref_alt_verdict = list()

	if (values.get('alt')):
		alt=values['alt']
	else:
		alt="-"

	if (values.get('pos')):
		pos=values['pos']
	else:
		pos="-"

	if (values.get('chromosome')):
		chrm=values['chromosome']
	else:
		chrm="-"

	if (values.get('ref')):
		ref=values['ref']
	else:
		ref="-"
	
	if (values.get('acmg_annotation')):
		verdict=values['acmg_annotation']['verdict']['ACMG_rules']['verdict']
		acmg_rules=values['acmg_annotation']['verdict']['ACMG_rules']
		classification=values['acmg_annotation']['verdict']['classifications']
	else:
		verdict="NA"
		acmg_rules="NA"
		classification="NA"

	if (values.get('original_variant')):
		original_variant=values['original_variant']
	else:
		original_variant="NA"

	chrm_pos_ref_alt_verdict.append(chrm)
	chrm_pos_ref_alt_verdict.append(pos)
	chrm_pos_ref_alt_verdict.append(ref)
	chrm_pos_ref_alt_verdict.append(alt)
	chrm_pos_ref_alt_verdict.append(original_variant)
	chrm_pos_ref_alt_verdict.append(verdict)
	chrm_pos_ref_alt_verdict.append(acmg_rules)
	chrm_pos_ref_alt_verdict.append(classification)

	return chrm_pos_ref_alt_verdict


def parse_data (orig_variant_list):
	chrm_pos_ref_alt = list()
	chr_ex = orig_variant_list[0]
	pos_ex = orig_variant_list[1]
	if orig_variant_list[2] == "":
		ref_ex = "-"
	else:
		ref_ex = orig_variant_list[2]
	
	if orig_variant_list[3] == "":
		alt_ex = "-"
	else:
		alt_ex = orig_variant_list[3]

	chrm_pos_ref_alt.append(chr_ex)
	chrm_pos_ref_alt.append(pos_ex)
	chrm_pos_ref_alt.append(ref_ex)
	chrm_pos_ref_alt.append(alt_ex)

	return chrm_pos_ref_alt


if data.__class__ == list:
	#print ("The is a list", "The length is", len(data))
	for values in data:										# data is a list of dictionaries
		results = extract_data(values)
		#print (results[0], results[1], results[2], results[3], results[5], results[6], results[7], file=output, sep="\t")
		# Extracting the chrom pos ref alt from the original variant key
		orig_variant = parse_data(results[4].split('-'))
		print (orig_variant[0], orig_variant[1], orig_variant[2], orig_variant[3], results[5], results[6], results[7], file=output, sep="\t")
else:
	results = extract_data(data)
	#print (results[0], results[1], results[2], results[3], results[5], results[6], results[7], file=output, sep="\t")
	# Extracting the chrom pos ref alt from the original variant key
	orig_variant = parse_data(results[4].split('-'))
	print (orig_variant[0], orig_variant[1], orig_variant[2], orig_variant[3], results[5], results[6], results[7], file=output, sep="\t")

output.close ()
