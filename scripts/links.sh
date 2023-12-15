#!/usr/bin/bash
# This script will make links for the .bam and .bai files

output_folder=$1
sample_list=$2

for samples in `cat ${sample_list}`
do 
	ln -s ${output_folder}/${samples}/${samples}.final.bam ./
	#touch ${output_folder}/${samples}/${samples}.final.bam.bai
	ln -s ${output_folder}/${samples}/${samples}.final.bam.bai ./
	#touch ${samples}.final.bam
	touch ${samples}.final.bam.bai
done		
