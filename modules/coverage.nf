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

process COVERAGE_WGS {
    tag "${Sample}"
    publishDir "${params.output}/${Sample}/" , mode: 'copy', pattern: '*txt'

    input:
    tuple val(Sample), file(final_bam), file(final_bai)

    output:
    file("${Sample}_WGS_coverage_metrics.txt")

    script:
    """
    ${params.java_path}/java -jar ${params.picard_path} CollectWgsMetrics \
    I=${final_bam} \
    R=${params.genome} \
    O=${Sample}_WGS_coverage_metrics.txt
    """
}

process COVERAGE_WGS_COLLECT {
	publishDir "${params.output}/", mode: 'copy'
	input:
		file (coverage_wgs)
	output:
		file("WGS_metrics.txt")
	script:
	"""
	echo -e "Sample name\tMean coverage" > WGS_metrics.txt
	for i in ${coverage_wgs}
	do
		samp_name=\$(basename -s -WGS_WGS_coverage_metrics.txt \${i})
		grep -v '#' \${i} | awk -v name=\${samp_name} 'BEGIN{FS="\t"; OFS="\t"} /^GENOME_TERRITORY/{getline; print name, \$2}' >> WGS_metrics.txt
	done
	"""
}
