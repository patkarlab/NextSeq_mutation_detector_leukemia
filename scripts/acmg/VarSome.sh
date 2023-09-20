#!/usr/bin/env bash

variants_file=$1
output_file=$2
acmg_outfile=$3

#varsome_api_run.py -g hg19 -k U1t?AcAG@A6nB@VLB1eTaAiATRZ0u3eX6#9qlTup -q 'chr17-7579472-G-A' -p add-all-data=1 -u https://staging-api.varsome.com/ -o annotations.txt
# The above search is translated to url: //lookup/chr17-7579472-G-A/hg19?add-all-data=1
#varsome_api_run.py -g hg19 -k 'U1t?AcAG@A6nB@VLB1eTaAiATRZ0u3eX6#9qlTup' -i ${variants_file}  -p add-ACMG-annotation=1 add-source-databases="none"  -u https://staging-api.varsome.com/ -o ${output_file}

source /home/programs/VarSomeApi/bin/activate
touch ${output_file}
if [ -s ${variants_file} ]; then
	signal=1
	while [ ${signal} -ne 0 ]; do
		rm ${output_file}
		varsome_api_run.py -g hg19 -k 'U1t?AcAG@A6nB@VLB1eTaAiATRZ0u3eX6#9qlTup' -i ${variants_file} -p add-ACMG-annotation=1 -u https://staging-api.varsome.com/ -o ${output_file}
		grep -i 'Connection failure' ${output_file}
		if [ $? -eq 0 ]; then
			signal=1
		else
			signal=0
		fi
	done		
fi
deactivate

# Extract the classification data to acmg_outfile
/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/acmg/json_file.py ${output_file} ${acmg_outfile}
