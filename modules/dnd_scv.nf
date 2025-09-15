#!/usr/bin/env nextflow

process DNDSCV {
	publishDir "$PWD/Final_Output/", mode : 'copy'
	input:
		path all_samples 
	output:
		path ("*_dndscv.xlsx")
	script:
	"""
	all_excels="${all_samples.join(' ')}"
	
	for excel in \$all_excels; do
		sampleID=`basename \$excel .xlsx`
		cp \$excel \${sampleID}_dndscv.xlsx
		run_dndscv.sh \${sampleID}_dndscv.xlsx \${sampleID}
	done
	"""
}
