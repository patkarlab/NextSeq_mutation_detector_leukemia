#!/usr/bin/bash

#genome=$1	# hg19_all.fasta
#vcf_file=$2	
#bam_file=$3	
#outfile=$4	# sample.html

#source activate cnvkit

#create_report ${vcf_file} --fasta ${genome} --ideogram /home/programs/IGV_reports/igv-reports/test/data/hg19/cytoBandIdeo.txt --info-columns GENE TISSUE TUMOR COSMIC_ID GENE SOMATIC --tracks ${vcf_file} ${bam_file} --output ${outfile}

#conda deactivate


source activate new_base

perl /home/programs/annovar_latest/annovar//table_annovar.pl ${vcf_file} --out temp.annovar --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring . --otherinfo --thread 10 /home/programs/annovar_latest/annovar//humandb/ --xreffile /home/programs/annovar_latest/annovar//example/gene_fullxref.txt -vcfinput

conda deactivate

source activate cnvkit

create_report ${vcf_file} --fasta ${genome} --ideogram /home/programs/IGV_reports/igv-reports/test/data/hg19/cytoBandIdeo.txt --info-columns Gene.refGene cosmic84 SOMATIC --tracks ${vcf_file} ${bam_file} /home/reference_genomes/refgene_hg37/refGene.txt.gz --output ${outfile}

conda deactivate
