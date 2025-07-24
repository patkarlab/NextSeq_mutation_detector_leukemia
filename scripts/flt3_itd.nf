#!/usr/bin/nextflow

process FILT3R {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_filt3r_out.csv'
	input:
		tuple val(Sample), file (bamin), file(baiin)
	output:
		tuple val (Sample), file("*_filt3r_out.csv")
	script:
	"""
	${params.bedtools} bamtofastq -i ${bamin} -fq ${Sample}_R1.fastq -fq2 ${Sample}_R2.fastq
	filt3r -k 12 --ref ${params.filt3r_ref} --sequences ${Sample}_R1.fastq,${Sample}_R2.fastq --nb-threads 64 --vcf --out ${Sample}_filt3r.json
	python3 /home/diagnostics/pipelines/Validation/scripts/convert_json_to_csv.py ${Sample}_filt3r.json ${Sample}_filt3r_json.csv
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}_filt3r.vcf --outfile ${Sample}.filt3r.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.filt3r.avinput --out ${Sample}.filt3r --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	# Check if the multianno file is empty
	if [[ ! -s ${Sample}.filt3r.hg19_multianno.csv ]]; then
		touch ${Sample}.filt3r__final.csv
		touch ${Sample}_filt3r_json_filtered.csv
		touch ${Sample}_filt3r_out.csv
	else
		python3 /home/diagnostics/pipelines/Validation/scripts/somaticseqoutput-format_filt3r.py ${Sample}.filt3r.hg19_multianno.csv ${Sample}.filt3r__final.csv
		python3 /home/diagnostics/pipelines/Validation/scripts/filter_json.py ${Sample}_filt3r_json.csv ${Sample}_filt3r_json_filtered.csv
		python3 /home/diagnostics/pipelines/Validation/scripts/merge_filt3r_csvs.py ${Sample}.filt3r__final.csv ${Sample}_filt3r_json_filtered.csv ${Sample}_filt3r_out.csv
	fi
	"""
}

process GETITD {
	tag "${Sample}"
	publishDir "${params.output}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		path "*_getitd"
	script:
	"""
	${params.samtools} sort ${finalBam} -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}