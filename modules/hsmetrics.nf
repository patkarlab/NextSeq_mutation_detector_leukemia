#!/usr/bin/env nextflow

process HSMETRICS {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_hsmetrics*.txt'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		path ("${Sample}_hsmetrics.txt"), emit : genewise
		path ("${Sample}_hsmetrics_exonwise.txt"), emit : exonwise
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics_exonwise.txt BAIT_INTERVALS= ${params.bedfile_exonwise}.interval_list TARGET_INTERVALS= ${params.bedfile_exonwise}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	"""
}

process HSMETRICS_COLLECT {
	publishDir "${params.output}/", mode: 'copy'
	input:
		file (GeneWise) 
		file (ExonWise)
	output:
		file("hsmetrics_probewise.txt") 
		file("hsmetrics_exonwise.txt")
	script:
	"""
	echo -e "Sample name\tOn target\tOff target" > hsmetrics_probewise.txt
	for i in ${GeneWise}
	do
		samp_name=\$(basename -s .txt \${i})
		grep -v '#' \${i} | awk -v name=\${samp_name} 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print name,\$7,\$8}' >> hsmetrics_probewise.txt
	done

	echo -e "Sample name\tOn target\tOff target" > hsmetrics_exonwise.txt
	for i in ${ExonWise}
	do
		samp_name=\$(basename -s .txt \${i})
		grep -v '#' \${i} | awk -v name=\${samp_name} 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print name,\$7,\$8}' >> hsmetrics_exonwise.txt
	done
	"""
}