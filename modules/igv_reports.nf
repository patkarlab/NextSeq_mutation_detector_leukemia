#!/usr/bin/env nextflow

process IGV_REPORTS {
	tag "${Sample}"
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno), file (cancervarMultianno)	
	output:
		val(Sample)
	script:
	"""
	perl ${params.annovarLatest_path}/table_annovar.pl ${somaticVcf} --out ${Sample}.annovar --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring . --otherinfo --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt -vcfinput 

	igv_reports.sh ${params.genome} ${Sample}.annovar.hg19_multianno.vcf $PWD/Final_Output/${Sample}/${Sample}.final.bam $PWD/Final_Output/${Sample}/${Sample}_igv.html
	"""		
}