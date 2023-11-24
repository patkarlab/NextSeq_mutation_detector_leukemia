#!/usr/bin/bash

genome=$1	# hg19_all.fasta
vcf_file=$2	
bam_file=$3	
outfile=$4	# sample.html

touch ${bam_file}
touch ${bam_file}.bai

#awk '{if ($1 ~ /^#/ ) print $0 ; else if ($0  ~ /Func.refGene=exonic/) print $0}' ${vcf_file} > filtered.vcf
/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/vcf_filter.py ${vcf_file} filtered.vcf

source activate cnvkit

create_report filtered.vcf --fasta ${genome} --ideogram /home/programs/IGV_reports/igv-reports/test/data/hg19/cytoBandIdeo.txt --info-columns Gene.refGene cosmic84 SOMATIC --tracks ${vcf_file} ${bam_file} /home/reference_genomes/refgene_hg37/refGene.txt.gz --output ${outfile}

conda deactivate
