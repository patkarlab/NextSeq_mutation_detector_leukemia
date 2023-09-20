#! /usr/bin/bash 

for i in `cat samplesheet.csv`
do 
	vep -i Final_Output/${i}/${i}.somaticseq.vcf --cache -o Final_Output/${i}/${i}_vep.txt --offline --tab --force_overwrite --af_1kg --af --af_gnomadg --pubmed --sift b --canonical --hgvs --shift_hgvs 1
	filter_vep -i Final_Output/${i}/${i}_vep.txt -o Final_Output/${i}/${i}_filtered.txt --filter "(CANONICAL is YES) and (AF < 0.01 or not AF)" --force_overwrite
	grep -v "##" Final_Output/${i}/${i}_filtered.txt > Final_Output/${i}/${i}_vep_delheaders.txt
done
