#!/usr/bin/env nextflow

process CNVKIT {
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*gene_scatter.pdf'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final-scatter.png'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final-diagram.pdf'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.final.cns"), file ("*.final.cnr"), file ("*gene_scatter.pdf"), file ("*.final-scatter.png"), file ("*.final-diagram.pdf") 
	script:
	"""
	
	#cnvkit.py batch ${finalBam} -r ${params.cnvkitRef} -m hybrid --drop-low-coverage --output-dir ${PWD}/${Sample}/cnvkit/ --diagram --scatter
	cnvkit.sh ${finalBam} ${params.cnvkitRef} ./
	# Commenting the following command for CNV myeloid panel, uncomment for the usual AL panel
	#/${params.gene_scatter}/custom_scatter_v3.py ${params.gene_scatter}/chr_list_all.txt ./${Sample}.final.cnr ./${Sample}.final.cns ${Sample}

	custom_scatter_chrwise.py ${params.gene_scatter_list}/chrwise_list.txt ./${Sample}.final.cnr ./${Sample}.final.cns ${Sample}_chr_
	"""
}

process ANNOT_SV {
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_AnnotSV.tsv'
	input:
		tuple val (Sample), file (finalCns), file (finalCnr), file (geneScatter), file (finalScatter), file (finalDiagram)
	output:
		tuple val (Sample), file ("*_AnnotSV.tsv")
	script:
	"""
	annotsv.sh ${finalCns} ${Sample}
	substitute_null.py ${Sample}_annotsv.tsv ${Sample}_AnnotSV.tsv
	"""
}

process IFCNV {
	input:
		val Sample
	output:
		//tuple val (Sample), file ("*.html"), file ("*.tsv")
		val (Sample)
	script:
	"""
	links.sh $PWD/Final_Output/ ${params.input}
	mkdir ifCNV
	ifcnv.sh ./ ${params.bedfile}.bed ifCNV

	# Making ifCNV's output directory for each sample
	for i in `cat ${params.input}`
	do
		if [ ! -d $PWD/Final_Output/\${i}/ifCNV ]; then
			mkdir $PWD/Final_Output/\${i}/ifCNV
		fi		
	done

	# Copying output of ifCNV to respective samples
	if [ -f ifCNV/ifCNV.tsv ]; then
		for i in `awk 'NR>1{print \$3}' ifCNV/ifCNV.tsv | awk 'BEGIN{FS="."}{print \$1}' | sort | uniq`	
		do
			cp ifCNV/\${i}.*.html $PWD/Final_Output/\${i}/ifCNV/	
		done
	fi
	"""
}