#! /usr/bin/bash

samplesheet=$1

source activate new_base
for samples in `cat ${samplesheet}`
do
	fastqc -f bam Final_Output/${samples}/${samples}.final.bam --extract --outdir Final_Output/${samples}/
	fastqc -f fastq sequences/${samples}_*R1_*.fastq.gz sequences/${samples}_*R2_*.fastq.gz --extract --outdir Final_Output/${samples}/
done
