#!/usr/bin/env nextflow

process ICHOR_CNA {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/${Sample}-WGS/" , mode: 'copy'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: "*_genomeWide_all_sols.pdf"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: "*_genomeWide.pdf"
	input:
		tuple val(Sample), file(final_bam), file(final_bai)
	output:
		tuple val(Sample), path("${Sample}_ichorCNA"), file("*_genomeWide_all_sols.pdf"), file("*_genomeWide.pdf")
	script:
	"""
	run_ichorCNA.sh ${Sample} ${final_bam} ${Sample}_ichorCNA
	cp ${Sample}_ichorCNA/${Sample}/${Sample}_genomeWide_all_sols.pdf ./
	cp ${Sample}_ichorCNA/${Sample}/${Sample}_genomeWide.pdf ./
	"""
}

process OFFTARGET_BAM_GEN {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	
	output:
		tuple val (Sample), file ("${Sample}_offtarget.bam"), file("${Sample}_offtarget.bam.bai")
	script:
	"""
	bedtools intersect -abam ${finalBam} -b ${params.bedfile}.bed -v > ${Sample}_offtarget.bam
	samtools index -@16 ${Sample}_offtarget.bam
	"""
}

process ICHORCNA_OFFTARGET {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/" , mode: 'copy'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: "*_genomeWide_all_sols.pdf"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: "*_genomeWide.pdf"
	input:
		tuple val(Sample), file(offtarget_bam), file(offtarget_bai)
	output:
		tuple val(Sample), path("${Sample}_offtarget"), file("*_genomeWide_all_sols.pdf"), file("*_genomeWide.pdf")
	script:
	"""
	awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,\$4,"-"}' ${params.bedfile}.bed > probes.bed
	awk 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,\$4,"-"}' ${params.bedfile_exonwise}.bed > exons.bed
	offtarget_ichorCNA.sh ${Sample} ${offtarget_bam} probes.bed exons.bed
	cp ${Sample}_offtarget/${Sample}/${Sample}_genomeWide_all_sols.pdf ./
	cp ${Sample}_offtarget/${Sample}/${Sample}_genomeWide.pdf ./
	"""
}	

