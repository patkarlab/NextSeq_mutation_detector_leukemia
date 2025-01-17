#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process cnvkit {
	publishDir "$PWD/Final_Output/${Sample}/", mode : 'copy'
	input:
		val(Sample)
	output:
		tuple val(Sample), file("*scatter.png"), file("*chr_gene_scatter.pdf")
	script:
	"""
	${params.cnvkit_path} ${params.sequences}/${Sample}_tumor.bam ${params.cnvkitRef} ./
	/${params.gene_scatter}/custom_scatter_chrwise.py ${params.gene_scatter_list}/chrwise_list.txt ./${Sample}_tumor.cnr ./${Sample}_tumor.cns ${Sample}_chr_
	"""
}

workflow {
	samples_ch = Channel.fromPath(params.input).splitCsv().flatten()
	// samples_ch.view()
	cnvkit(samples_ch)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}