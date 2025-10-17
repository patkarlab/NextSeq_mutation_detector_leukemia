#!/usr/bin/env nextflow

process ICHOR_CNA {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/${Sample}-WGS/" , mode: 'copy'
	input:
	tuple val(Sample), file(final_bam), file(final_bai)

	output:
	tuple val(Sample), path("${Sample}_ichorCNA")

	script:
	"""
	run_ichorCNA.sh ${Sample} ${final_bam} ${Sample}_ichorCNA
	"""
}
