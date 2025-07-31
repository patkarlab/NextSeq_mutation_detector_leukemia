#!/usr/bin/env nextflow

process SOMATICSEQ {
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.somaticseq.vcf'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file(oldfinalBam), file(oldfinalBamBai) , file(mutectVcf), file(vardictVcf), file(DeepSomaticVcf), file(lofreqVcf), file(strelkaVcf), file(freebayesVcf), file(platypusVcf)
	output:
		tuple val (Sample), file("*.somaticseq.vcf"), file("*.hg19_multianno.csv"), file("*.hg19_multianno.txt.cancervar.ensemble.pred")
	script:
	"""
	${params.vcf_sorter_path} ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${platypusVcf} ${Sample}.platypus.sorted.vcf
	${params.vcf_sorter_path} ${DeepSomaticVcf} ${Sample}.deepsomatic.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.deepsomatic.sorted.vcf -snv ${Sample}_deepsomatic_snvs.vcf -indel ${Sample}_deepsomatic_indels.vcf

	${params.vcf_sorter_path} ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_deepsomatic_snvs.vcf ${Sample}_deepsomatic_snvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_deepsomatic_indels.vcf ${Sample}_deepsomatic_indels_sort.vcf

	somaticseq_parallel.py --output-directory ./${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --algorithm xgboost  --dbsnp-vcf  /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${mutectVcf} --vardict-vcf ${vardictVcf} --lofreq-vcf ${lofreqVcf} --strelka-vcf ${strelkaVcf} --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypus_cnvs_sort.vcf ${Sample}_deepsomatic_snvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf ${Sample}_deepsomatic_indels_sort.vcf
	
	${params.vcf_sorter_path} ./${Sample}.somaticseq/Consensus.sSNV.vcf ./${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ./${Sample}.somaticseq/somaticseq_snv.vcf > ./${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ./${Sample}.somaticseq/somaticseq_snv.vcf.gz
	
	${params.vcf_sorter_path} ./${Sample}.somaticseq/Consensus.sINDEL.vcf ./${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ./${Sample}.somaticseq/somaticseq_indel.vcf > ./${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ./${Sample}.somaticseq/somaticseq_indel.vcf.gz
	
	${params.bcftools_path} concat -a ./${Sample}.somaticseq/somaticseq_snv.vcf.gz ./${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ./${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=MDLK012,Number=7,Type=Integer,Description="Calling decision of the 7 algorithms: MuTect, VarDict, LoFreq, Strelka, SnvCaller_0, SnvCaller_1, SnvCaller_2">/##INFO=<ID=MDLKFPGS,Number=7,Type=String,Description="Calling decision of the 7 algorithms: MuTect, VarDict, LoFreq, Strelka, Freebayes, Platypus, DeepSomatic">/g' ${Sample}.somaticseq.vcf

	sed -i 's/MDLK012/MDLKFPS/g' ${Sample}.somaticseq.vcf
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf  --outfile ${Sample}.somaticseq.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	${params.cancervar} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}
	"""
}

process COMBINE_VARIANTS{
	tag "${Sample}"
	input:
		tuple val (Sample), file(mutectVcf), file(vardictVcf), file(DeepSomaticVcf), file(lofreqVcf), file(strelkaVcf), file(freebayesVcf), file(platypusVcf)
	output:
		tuple val(Sample), file("${Sample}.combined.vcf")
	script:
	"""
	${params.vcf_sorter_path} ${mutectVcf} ${Sample}.mutect.sorted.vcf
	${params.vcf_sorter_path} ${vardictVcf} ${Sample}.vardict.sorted.vcf
	${params.vcf_sorter_path} ${DeepSomaticVcf} ${Sample}.DeepSomatic.sorted.vcf
	${params.vcf_sorter_path} ${lofreqVcf} ${Sample}.lofreq.sorted.vcf
	${params.vcf_sorter_path} ${strelkaVcf} ${Sample}.strelka.sorted.vcf
	${params.vcf_sorter_path} ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${platypusVcf} ${Sample}.platypus.sorted.vcf

	${params.java_path}/java -jar ${params.GATK38_path} -T CombineVariants -R ${params.genome} --variant ${Sample}.mutect.sorted.vcf --variant ${Sample}.vardict.sorted.vcf --variant ${Sample}.DeepSomatic.sorted.vcf --variant ${Sample}.lofreq.sorted.vcf --variant ${Sample}.strelka.sorted.vcf --variant ${Sample}.freebayes.sorted.vcf --variant ${Sample}.platypus.sorted.vcf -o ${Sample}.combined.vcf -genotypeMergeOptions UNIQUIFY
	"""
}