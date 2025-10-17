#!/usr/bin/nextflow

process TRIM {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val(Sample), file(read1), file(read2)
	output:
		tuple val(Sample), file("${Sample}_trim_R1.fastq"), file("${Sample}_trim_R2.fastq")
	script:
	"""
	${params.fastp} -i ${read1} -I ${read2} -o ${Sample}_trim_R1.fastq -O ${Sample}_trim_R2.fastq --adapter_fasta ${params.adaptors} -w 4
	"""
}

process MAPBAM {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val(Sample), file(trim1), file(trim2)
	output:
		tuple val(Sample), file ("${Sample}.bam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" \
	-M -t 8 ${params.genome} ${trim1} ${trim2} | ${params.samtools} sort -@ 8 -o ${Sample}.bam -
	"""
}

// process SORT {
// 	tag "${Sample}"
// 	input:
// 		tuple val(Sample), file(bamfile)
// 	output:
// 		tuple val(Sample), file ("${Sample}_sortd.bam")
// 	script:
// 	"""
// 	samtools sort -o ${Sample}_sortd.bam  ${bamfile}
// 	"""
// }

process MARK_DUPS {
	tag "${Sample}"
	input:
		tuple val(Sample), file(sortd_bam)
	output:
		tuple val(Sample), file("${Sample}_markdups.bam"), file("${Sample}_marked_dup_metrics.txt")
	script:
	"""
	${params.java_path}/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx72G \
		-jar ${params.picard_path} MarkDuplicates \
		I=${sortd_bam} \
		O=${Sample}_markdups.bam \
		M=${Sample}_marked_dup_metrics.txt \
		TMP_DIR=${Sample}_tmp
	"""
}

process BQSR {
	tag "${Sample}"
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics)
	output:
		tuple val(Sample), file("${Sample}_recal.table")
	script:
	"""
	${params.gatk} BaseRecalibrator \
		-I ${markdups_bam} \
		-R ${params.genome} \
		--known-sites ${params.site1} \
		--known-sites ${params.site2} \
		--bqsr-baq-gap-open-penalty 30.0 \
		-O ${Sample}_recal.table
	"""
}

process APPLY_BQSR {
	// publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_final.bam'
	// publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_final.bam.bai'
	tag "${Sample}"
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics), file(recal_table)
	output:
		tuple val(Sample), file("${Sample}_final.bam"), file("${Sample}_final.bam.bai")
	script:
	"""
	${params.gatk} ApplyBQSR \
		-R ${params.genome} \
		-I ${markdups_bam} \
		--bqsr-recal-file ${recal_table} \
		-O ${Sample}_final.bam

	mv ${Sample}_final.bai ${Sample}_final.bam.bai
	"""
}

process ALIGNMENT_METRICS {
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_alignment_summary_metrics.txt'
	tag "${Sample}"
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai)
	output:
		tuple val(Sample), file("${Sample}_alignment_summary_metrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectAlignmentSummaryMetrics \
		R=${params.genome} \
		I=${final_bam} \
		O=${Sample}_alignment_summary_metrics.txt
	"""
}

process INSERT_SIZE_METRICS {
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_insert_size_metrics.txt'
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_insert_size_metrics.pdf'
	tag "${Sample}"
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai)
	output:
		tuple val(Sample), file("${Sample}_insert_size_metrics.txt"), file("${Sample}_insert_size_metrics.pdf")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectInsertSizeMetrics \
		I=${final_bam} \
		O=${Sample}_insert_size_metrics.txt \
		H=${Sample}_insert_size_metrics.pdf \
		HISTOGRAM_WIDTH=500 \
		TMP_DIR=${Sample}_tmp
	"""
}

process ALIGNMENT_METRICS_WGS {
	publishDir "${params.output}/${Sample}/${Sample}-WGS/", mode: 'copy', pattern: '*_alignment_summary_metrics.txt'
	tag "${Sample}"
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai)
	output:
		tuple val(Sample), file("${Sample}_alignment_summary_metrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectAlignmentSummaryMetrics \
		R=${params.genome} \
		I=${final_bam} \
		O=${Sample}_alignment_summary_metrics.txt
	"""
}

process INSERT_SIZE_METRICS_WGS {
	publishDir "${params.output}/${Sample}/${Sample}-WGS/", mode: 'copy', pattern: '*_insert_size_metrics.txt'
	publishDir "${params.output}/${Sample}/${Sample}-WGS/", mode: 'copy', pattern: '*_insert_size_metrics.pdf'
	tag "${Sample}"
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai)
	output:
		tuple val(Sample), file("${Sample}_insert_size_metrics.txt"), file("${Sample}_insert_size_metrics.pdf")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectInsertSizeMetrics \
		I=${final_bam} \
		O=${Sample}_insert_size_metrics.txt \
		H=${Sample}_insert_size_metrics.pdf \
		HISTOGRAM_WIDTH=500 \
		TMP_DIR=${Sample}_tmp
	"""
}

process ABRA_BAM { 
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	input:
		tuple val(Sample), file (bamin), file(baiin)
	output:
		tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai"), file ("*.old_final.bam"), file ("*.old_final.bam.bai")
	script:
	"""
	${params.bedtools} sort -i ${params.bedfile}.bed > sorted.bed

	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${bamin} --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --targets sorted.bed --tmpdir ./ > abra.log

	${params.samtools} sort ${bamin} > ${Sample}.old_final.bam
	${params.samtools} index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	${params.samtools} sort ${Sample}.abra.bam > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

workflow FASTQTOBAM {
	take:
		samples_ch
	main:
	 TRIM(samples_ch)
	 MAPBAM(TRIM.out)
	//  SORT(MAPBAM.out)
	 MARK_DUPS(MAPBAM.out)
	 BQSR(MARK_DUPS.out)
	 APPLY_BQSR(MARK_DUPS.out.join(BQSR.out))
	 ALIGNMENT_METRICS(APPLY_BQSR.out)
	 INSERT_SIZE_METRICS(APPLY_BQSR.out)
	emit:
		final_bams_ch = APPLY_BQSR.out	
}

workflow FASTQTOBAM_WGS {
	take:
		samples_ch
	main:
	 TRIM(samples_ch)
	 MAPBAM(TRIM.out)
	 MARK_DUPS(MAPBAM.out)
	 BQSR(MARK_DUPS.out)
	 APPLY_BQSR(MARK_DUPS.out.join(BQSR.out))
	 ALIGNMENT_METRICS_WGS(APPLY_BQSR.out)
	 INSERT_SIZE_METRICS_WGS(APPLY_BQSR.out)
	emit:
		wgs_bam_ch = APPLY_BQSR.out
}