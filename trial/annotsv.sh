#!/usr/bin/bash

input=$1
#output_tsvfile=$2

export PATH=$PATH:/home/programs/MoChA/bcftools-1.16
#AnnotSV -SVinputFile ${input} -outputFile ${output_tsvfile} -svtBEDcol 4 -genomeBuild GRCh37	# for bedfile as input
#AnnotSV -SVinputFile ${input} -outputFile ./${output_tsvfile} -genomeBuild GRCh37	# for vcf as input

for samples in `cat ${input}`
do
	#ls ../${samples}/cnvkit/${samples}.final.cns
	echo -e "# chrom\tStart\tEnd" > ${samples}.bed
	awk 'BEGIN{OFS="\t"}NR>1{if($5 > 0.4 || $5 < -0.4) print $1,$2,$3}' ../${samples}/cnvkit/${samples}.final.cns >> ${samples}.bed
	AnnotSV -SVinputFile ${samples}.bed -outputFile ./${samples}.tsv -genomeBuild GRCh37
	rm ${samples}.bed
done
