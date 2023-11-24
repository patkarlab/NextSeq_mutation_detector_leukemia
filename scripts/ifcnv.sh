#!/usr/bin/bash
# This script will run ifCNV on all the bam files 

input_dir=$1
bedfile=$2
output_dir=$3

source activate ifcnv

ifCNV -i ${input_dir} -b ${bedfile} -o ${output_dir} -sv True

conda deactivate
