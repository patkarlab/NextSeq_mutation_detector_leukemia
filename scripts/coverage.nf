#!/usr/bin/env nextflow

process COVERAGE {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed"), file ("${Sample}_pindel.counts.bed"), file("${Sample}_pindel_ubtf.counts.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBams[0]} > ${Sample}.bed
	${params.bedtools} coverage -counts -a ${params.bedfile_exonwise}.bed -b ${Sample}.bed > ${Sample}.counts.bed
	${params.bedtools} coverage -counts -a ${params.flt3_bedfile}.bed -b ${Sample}.bed > ${Sample}_pindel.counts.bed
	grep 'UBTF' ${params.bedfile}.bed > ubtf_pindel.bed
	${params.bedtools} coverage -counts -a ubtf_pindel.bed -b ${Sample}.bed > ${Sample}_pindel_ubtf.counts.bed
	"""
}

process COVERVIEW {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: "*.coverview_regions.csv"
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.coverview_regions.csv")
	script:
	"""
	${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile_exonwise}.bed -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	"""
}