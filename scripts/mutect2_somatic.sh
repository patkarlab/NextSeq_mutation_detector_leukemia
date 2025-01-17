#!/usr/bin/env bash

java_path=$1
gatk_path=$2
genome=$3
normal_bam=$4
tumor_bam=$5
output_file=$6
dbsnp=$7
bedfile=$8

source activate new_base
${java_path}/java -Xmx10G -jar ${gatk_path} -T MuTect2 -R ${genome} -I:tumor ${tumor_bam} -I:normal ${normal_bam} -o ${output_file} --dbsnp ${dbsnp} -L ${bedfile} -nct 64 -mbq 25
#conda deactivate
