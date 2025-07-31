#!/usr/bin/env nextflow


process PLATYPUS {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file ("*.platypus.vcf")
	script:
	"""
	python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBams[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --minMapQual=50 --maxVariants=6 --minReads=6 --regions=${params.bedfile}_regions.txt
	"""
}

process FREEBAYES {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")
	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBam} -t ${params.bedfile}.bed > ${Sample}.freebayes.vcf 	
	"""
}

process MUTECT2 {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.mutect2.vcf")
	script:
	"""
	${params.samtools} view -bs 40.1 ${finalBam} > subsampled_01.bam
	${params.samtools} index subsampled_01.bam
	${params.mutect2} ${params.java_path} ${params.GATK38_path} ${params.genome} subsampled_01.bam ${Sample}.mutect2.vcf ${params.site2} ${params.bedfile}.bed
	"""
}

process VARDICT {
	tag "${Sample}"	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${finalBam} -O 50 -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf
	"""
}

process DEEPSOMATIC {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file("*_DS.vcf")
	script:
	""" 
	mkdir output
	outpath=`realpath output`
	bam_path=`realpath ${finalBam} | awk 'BEGIN{OFS=FS="/"} {\$NF=""; print \$0}'`
	vcf_output=${Sample}_DS.vcf
	control_bam_path=`realpath /home/programs/DeepSomatic/reads/OCIAML3.bam`
	echo \$bam_path \$outpath \$vcf_output ${finalBam} \${control_bam_path} ${params.genome} ${params.bedfile}.bed
	${params.deepsomatic} \$bam_path \$outpath \$vcf_output ${finalBam} \${control_bam_path} ${params.genome} ${params.bedfile}.bed
	"""
}

process LOFREQ {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file ("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${oldfinalBam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.005 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process STRELKA {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.strelka.vcf")
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBam} --referenceFasta ${params.genome} --callRegions  ${params.bedfile}.bed.gz --targeted --runDir ./
	./runWorkflow.py -m local -j 20
	gunzip -f ./results/variants/variants.vcf.gz
	mv ./results/variants/variants.vcf ./${Sample}.strelka.vcf
	"""
}

process PINDEL {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*_pindel.hg19_multianno.csv")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=${params.pindel}/Adaptor.pm
	printf '%s\t%s\t%s' ${finalBam} "300" ${Sample} > ./config.txt
	${params.pindel}/pindel -f ${params.genome} -i ./config.txt -c chr13 -o ${Sample}_pindel
	${params.pindel}/pindel2vcf -r ${params.genome} -P ${Sample}_pindel -R hg19 -d 07102019 -v ${Sample}_pindel_SI.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}_pindel_SI.vcf --outfile ${Sample}_pindel.avinput --withzyg --includeinfo

	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}_pindel.avinput ${params.annovarLatest_path}/humandb/ -buildver hg19 -out ${Sample}_pindel --remove -protocol refGene,cytoBand,cosmic84 --operation g,r,f -nastring '.' --otherinfo --csvout --thread 10 --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
 	"""
}

process PINDEL_UBTF {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*_pindel_ubtf.hg19_multianno.csv")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=${params.pindel}/Adaptor.pm
	printf '%s\t%s\t%s' ${finalBam} "300" ${Sample} > ./config.txt
	${params.pindel}/pindel -f ${params.genome} -i ./config.txt -c chr17 -o ${Sample}_pindel_ubtf
	${params.pindel}/pindel2vcf -r ${params.genome} -P ${Sample}_pindel_ubtf -R hg19 -d 07102019 -v ${Sample}_pindel_ubtf_SI.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}_pindel_ubtf_SI.vcf --outfile ${Sample}_pindel_ubtf.avinput --withzyg --includeinfo

	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}_pindel_ubtf.avinput ${params.annovarLatest_path}/humandb/ -buildver hg19 -out ${Sample}_pindel_ubtf --remove -protocol refGene,cytoBand,cosmic84 --operation g,r,f -nastring '.' --otherinfo --csvout --thread 10 --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
 	"""	
}