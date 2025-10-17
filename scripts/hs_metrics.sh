#! /usr/bin/bash

for i in `cat samplesheet.csv` 
do 
	#java -jar /home/programs/picard/build/libs/picard.jar CollectHsMetrics I=/home/pipelines/NextSeq_mutation_detector_leukemia/Final_Output/$i/$i".final.bam" O=/home/pipelines/NextSeq_mutation_detector_leukemia/Final_Output/$i".hsmetrics.txt" R=/home/reference_genomes/hg19_broad/hg19_all.fasta BAIT_INTERVALS=/home/pipelines/mutation_detector_nextflow/bedfile/06112021_Leukemia_Panel_sorted.intervals.list TARGET_INTERVALS=/home/pipelines/mutation_detector_nextflow/bedfile/06112021_Leukemia_Panel_sorted.intervals.list VALIDATION_STRINGENCY=LENIENT
	sample=${i%%-*}
	echo -ne $i'\t'; grep -v '#' /home/pipelines/NextSeq_mutation_detector_leukemia/Final_Output/$sample/$sample"_hsmetrics.txt" | awk 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print $7,$8}'

done
