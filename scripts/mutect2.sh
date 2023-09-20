#!/usr/bin/env bash

java_path=$1
gatk_path=$2
genome=$3
bam_file=$4
output_file=$5
dbsnp=$6
bedfile=$7

source activate new_base
${java_path}/java -Xmx10G -jar ${gatk_path} -T MuTect2 -R ${genome} -I:tumor ${bam_file} -o ${output_file} --dbsnp ${dbsnp} -L ${bedfile} -nct 25 -mbq 25
#conda deactivate
