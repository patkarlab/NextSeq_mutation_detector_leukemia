#!/usr/bin/bash

cns_file=$1
samples=$2

export PATH=$PATH:/home/programs/MoChA/bcftools-1.16
#AnnotSV -SVinputFile ${input} -outputFile ${output_tsvfile} -svtBEDcol 4 -genomeBuild GRCh37	# for bedfile as input
#AnnotSV -SVinputFile ${input} -outputFile ./${output_tsvfile} -genomeBuild GRCh37	# for vcf as input

echo -e "# chrom\tStart\tEnd" > ${samples}.bed
awk 'BEGIN{OFS="\t"}NR>1{if($5 > 0.4 || $5 < -0.4) print $1,$2,$3}' ${cns_file} >> ${samples}.bed
no_of_line=$(wc -l ${samples}.bed | awk '{print $1}')
#echo ${no_of_line}

if [ ${no_of_line} -gt 1 ];then
	AnnotSV -SVinputFile ${samples}.bed -outputFile ./${samples}_annotsv.tsv -genomeBuild GRCh37
	rm ${samples}.bed
else
	touch ./${samples}_annotsv.tsv
fi
