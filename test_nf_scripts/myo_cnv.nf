#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MARKDUPLICATES {
	tag "${Sample}"
	publishDir "/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_myeloidpanel_240425_markduplicates/", mode: 'copy', pattern: '*'
	input:
		val (Sample)
	output:
		tuple val(Sample), file ("${Sample}_MD.bam"), file ("${Sample}_MD.bai")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} MarkDuplicates \
	INPUT=${params.sequences}/${Sample}_tumor.bam OUTPUT=${Sample}_MD.bam REMOVE_DUPLICATES=false \
	METRICS_FILE=BNC1.txt CREATE_INDEX=true
	"""
}

process cnvkit {
	publishDir "$PWD/Final_Output/${Sample}/", mode : 'copy'
	input:
		val(Sample)
	output:
		tuple val(Sample), file("*scatter.png"), file("*chr_gene_scatter.pdf")
	script:
	"""
	#${params.java_path}/java -jar ${params.picard_path} MarkDuplicates \
	#INPUT=${params.sequences}/${Sample}_tumor.bam OUTPUT=${Sample}_MD.bam REMOVE_DUPLICATES=false \
	#METRICS_FILE=BNC1.txt CREATE_INDEX=true

	${params.cnvkit_path} ${params.sequences}/${Sample}_tumor.bam ${params.cnvkitRef} ./
	#${params.cnvkit_path} ${Sample}_MD.bam ${params.cnvkitRef} ./
	${params.gene_scatter}/custom_scatter_chrwise.py ${params.gene_scatter_list}/chrwise_list.txt ./${Sample}_tumor.cnr ./${Sample}_tumor.cns ${Sample}_chr_
	"""
}

process hsmetrics_run{
	// publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_hsmetrics*.txt'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai)
	output:
		path ("${Sample}_hsmetrics.txt") , emit : genewise
		path ("${Sample}_hsmetrics_exonwise.txt"), emit : exonwise
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics_exonwise.txt BAIT_INTERVALS= ${params.bedfile_exonwise}.interval_list TARGET_INTERVALS= ${params.bedfile_exonwise}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT

	"""
}

process generatefinalbam{
	// publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	// publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	input:
		val Sample
	output:
		tuple val(Sample), file ("*.old_final.bam"), file ("*.old_final.bam.bai")
		
	script:
	"""
	#${params.bedtools} sort -i ${params.bedfile}.bed > sorted.bed
	#${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${params.sequences}/${Sample}_tumor.bam --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --targets sorted.bed --tmpdir ./ > abra.log
	${params.samtools} sort ${params.sequences}/${Sample}_tumor.bam > ${Sample}.old_final.bam
	${params.samtools} index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	#${params.samtools} sort ${Sample}.tumor.bam > ${Sample}.final.bam
	#${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

process hsmetrics_collect {
	input:
		file (GeneWise) 
		file (ExonWise)
	script:
	"""
	echo -e "Sample name\tOn target\tOff target" > hsmetrics_genwise.txt
	for i in ${GeneWise}
	do
		samp_name=\$(basename -s .txt \${i})
		grep -v '#' \${i} | awk -v name=\${samp_name} 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print name,\$7,\$8}' >> hsmetrics_genwise.txt
	done
	#ls ${ExonWise} 
	"""
}

process DNDSCV {
	publishDir "$PWD/Final_Output/", mode : 'copy'
	input:
		val (Sample)
	output:
		tuple val(Sample), file ("${Sample}_dndscv.xlsx")
	script:
	"""
	cp $PWD/Final_Output/${Sample}/${Sample}.xlsx ${Sample}_dndscv.xlsx
	/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/dNdScv/run_dndscv.sh ${Sample}_dndscv.xlsx ${Sample}
	"""
}

workflow CNV {
	samples_ch = Channel.fromPath(params.input).splitCsv().flatten()
	// MARKDUPLICATES(samples_ch)
	// cnvkit(samples_ch)
	// generatefinalbam(samples_ch)
	// hsmetrics_run(generatefinalbam.out)
	// Run the perSampleProcess for each sample
	// ch_results_genewise = hsmetrics_run.out.genewise.collect()
	// ch_results_exon = hsmetrics_run.out.exonwise.collect()

	// hsmetrics_collect(ch_results_genewise, ch_results_exon)
	//hsmetrics_collect(ch_results_exon, "Exonwise")

	// Collect all output files before processing them together
	DNDSCV(samples_ch)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
