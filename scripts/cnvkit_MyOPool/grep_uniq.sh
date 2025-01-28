#!/usr/bin/bash


while IFS='' read -r line ; do
	grep -m1 -w "$line" /home/pipelines/mutation_detector_nextflow/bedfile/MYOPOOL_231224.bed 
done < "$1"
