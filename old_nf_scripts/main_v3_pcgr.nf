#!/usr/bin/env nextflow
nextflow.enable.dsl=2

"mkdir Coverview".execute()


log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=

Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}

"""


process trimming_fastq_mcf {
	maxForks 15
	publishDir "${PWD}/${Sample}/processed_reads/", mode: 'copy'
	input:
		val (Sample)
	output:
		tuple val (Sample), file ("*.fastq")
	script:
	"""
	${params.ea_utils_path}/fastq-mcf -o ${Sample}.R1.trimmed.fastq -o ${Sample}.R2.trimmed.fastq -l 53 -k 0 -q 0 ${params.adaptors} ${params.sequences}/${Sample}*_R1_*.fastq.gz ${params.sequences}/${Sample}*_R2_*.fastq.gz
	sleep 5s
	"""
}

process gzip{
	publishDir "${PWD}/${Sample}/processed_reads/", mode: 'copy'
	input:
		tuple val (Sample), file(trimmedFiles)
	output:
		val Sample
	script:
	"""
	mkdir "$PWD/${Sample}/Annovar_Modified/" 
	gzip -f ${PWD}/${Sample}/processed_reads/${trimmedFiles[0]}
	gzip -f ${PWD}/${Sample}/processed_reads/${trimmedFiles[1]}
	sleep 5s
	"""
}

process trimming_trimmomatic {
	maxForks 10
	publishDir "$PWD/${Sample}/trimmed", mode: 'copy'
	input:
			val Sample
	output:
			tuple val (Sample), file("*.fq.gz")
	script:
	"""
	trimmomatic PE \
	${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${params.adaptors}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
	sleep 5s
	"""
}

process minimap_getitd {
	publishDir "$PWD/${Sample}/", mode: 'copy', pattern: '*_getitd'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
			val Sample
	output:
			path "*_getitd"
	script:
	"""
	minimap2 -ax sr ${params.genome_minimap_getitd} ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz > ${Sample}.sam
	${params.samtools} view -b -h ${Sample}.sam -o ${Sample}.bam
	${params.samtools} sort ${Sample}.bam -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process pair_assembly_pear {
	memory '7.0 GB'
	publishDir "${PWD}/${Sample}/assembled_reads/", mode: 'copy'
	publishDir "${PWD}/${Sample}/Annovar_Modified/", mode: 'copy'
	input:
		tuple val (Sample), file(trimmedFiles)
	output:
		tuple val (Sample), file("*") 
	script:
	"""
	sleep 5s
	${params.pear_path} -f ${PWD}/${Sample}/trimmed/${trimmedFiles[0]} -r ${PWD}/${Sample}/trimmed/${trimmedFiles[2]} -o ${Sample} -n 53 -j 25
	"""
}

process mapping_reads{
	maxForks 15
	publishDir "${PWD}/${Sample}/mapped_reads/", mode: 'copy'
	input:
		tuple val (Sample), file (pairAssembled)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled[0]} > ${Sample}.sam
	"""
} 



process sam_conversion{
	maxForks 15
	publishDir "$PWD/${Sample}/mapped_reads/", mode: 'copy', pattern: '*.fxd_sorted.bam'
	publishDir "$PWD/${Sample}/mapped_reads/", mode: 'copy', pattern: '*.fxd_sorted.bam.bai'

	input:
		tuple val (Sample), file(samFile)
	output:
		tuple val(Sample), file ("*.fxd_sorted.bam"), file ("*.fxd_sorted.bam.bai")
	
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} FixMateInformation I= ${samFile} O= ${Sample}.fxd.sam VALIDATION_STRINGENCY=SILENT
	cp ${samFile} ${Sample}.fxd.sam
	${params.samtools} view -bT ${params.genome} ${Sample}.fxd.sam > ${Sample}.fxd.bam
	${params.samtools} sort ${Sample}.fxd.bam > ${Sample}.fxd_sorted.bam
	${params.samtools} index ${Sample}.fxd_sorted.bam > ${Sample}.fxd_sorted.bam.bai
	"""
}

process RealignerTargetCreator {
	publishDir "${PWD}/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.intervals'
	
	input:
		tuple val (Sample), file (bamFile), file(bamBai)
	output:
		tuple val(Sample), file ("*.intervals")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFile} --known ${params.site1} -o ${Sample}.intervals
	"""
}

process IndelRealigner{
	publishDir "${PWD}/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.realigned.bam'
	input:
		tuple val(Sample), file (targetIntervals), file(bamFile), file(bamBai)
	output:
		tuple val(Sample), file ("*.realigned.bam")
	script:
	"""
	echo ${Sample} ${targetIntervals} ${bamFile}
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T IndelRealigner -R ${params.genome} -I ${bamFile} -known ${params.site1} --targetIntervals ${targetIntervals} -o ${Sample}.realigned.bam
	"""
}

process BaseRecalibrator{
	publishDir "${PWD}/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.recal_data.table'
	input:
		tuple val (Sample), file (realignedBam)
	output:
		tuple val(Sample), file ("*.recal_data.table")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T BaseRecalibrator -R ${params.genome} -I ${realignedBam} -knownSites ${params.site2} -knownSites ${params.site3} -maxCycle 600 -o ${Sample}.recal_data.table
	"""
}

process PrintReads{
	publishDir "${PWD}/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.aligned.recalibrated.bam'
	input:
		tuple val (Sample), file (realignedBam), file (recal_dataTable)
	output:
		tuple val (Sample), file ("*.aligned.recalibrated.bam")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T PrintReads -R ${params.genome} -I ${realignedBam} --BQSR ${recal_dataTable} -o ${Sample}.aligned.recalibrated.bam
	"""
}

process generatefinalbam{
	publishDir "$PWD/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/${Sample}/gatk38_processing/", mode: 'copy', pattern: '*.final.bam.bai'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	
	input:
		tuple val (Sample), file(alignedRecalibratedBam)
	output:
		tuple val(Sample), file ("*.final.bam"),  file ("*.final.bam.bai")
		
	script:
	"""
	${params.samtools} sort ${alignedRecalibratedBam} > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

process hsmetrics_run{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_hsmetrics.txt'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.picard_interval} TARGET_INTERVALS= ${params.picard_interval} R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	"""
}

process mutect2_run{
	maxForks 10
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.mutect2.vcf'
	
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*.mutect2.vcf")

	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${finalBam} -o ${Sample}.mutect2.vcf --dbsnp ${params.site2} -L ${params.bedfile}.bed -nct 25 -contamination 0.02 -mbq 30
	"""
}

process freebayes_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.freebayes.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")

	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBam} -t ${params.bedfile}.bed > ${Sample}.freebayes.vcf 	
	"""
}


process vardict_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.vardict.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${finalBam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf
	"""
}

process varscan_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_snp.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_indel.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_snp.vcf.gz'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan_indel.vcf.gz'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.varscan.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val(Sample), file ("*.varscan_snp.vcf"),  file ("*.varscan_indel.vcf"), file("*.varscan.vcf")
		
	script:
	"""
	${params.samtools} mpileup -f ${params.genome} ${finalBam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process lofreq_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.lofreq.filtered.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val(Sample), file ("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${finalBam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.005 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process strelka_run{
	publishDir "$PWD/${Sample}/variants/strelka", mode: 'copy', pattern: '*.strelka.vcf'
	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		val (Sample)
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${finalBam} --referenceFasta ${params.genome} --callRegions  ${params.bedfile}.bed.gz --targeted --runDir ${PWD}/${Sample}/variants/strelka/
	${PWD}/${Sample}/variants/strelka/runWorkflow.py -m local -j 20
	gunzip -f ${PWD}/${Sample}/variants/strelka/results/variants/variants.vcf.gz
	mv ${PWD}/${Sample}/variants/strelka/results/variants/variants.vcf $PWD/${Sample}/variants/${Sample}.strelka.vcf
	
	${params.strelka_path}/configureStrelkaSomaticWorkflow.py --normalBam ${params.NA12878_bam}  --tumorBam ${finalBam} --referenceFasta ${params.genome} --callRegions ${params.bedfile}.bed.gz --targeted --runDir ${PWD}/${Sample}/variants/strelka-somatic/
	${PWD}/${Sample}/variants/strelka-somatic/runWorkflow.py -m local -j 20
	
	${params.bcftools_path} concat -a ${PWD}/${Sample}/variants/strelka-somatic/results/variants/somatic.indels.vcf.gz ${PWD}/${Sample}/variants/strelka-somatic/results/variants/somatic.snvs.vcf.gz -o ${PWD}/${Sample}/variants/${Sample}.strelka-somatic.vcf
	"""
}

process somaticSeq_run {
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.somaticseq.vcf'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample}/ANNOVAR/", mode: 'copy', pattern: '*.hg19_multianno.csv'
	input:
		tuple val (Sample), file(mutectVcf), file(vardictVcf), file(varscanVcf), file(lofreqVcf)
	output:
		tuple val (Sample), file ("*.somaticseq.vcf"), file("*.hg19_multianno.csv")
	script:
	"""
	${params.vcf_sorter_path} ${PWD}/${Sample}/variants/${Sample}.freebayes.vcf ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${PWD}/${Sample}/variants/${Sample}.platypus.vcf ${Sample}.platypus.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf

	${params.vcf_sorter_path} ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf

	somaticseq_parallel.py --output-directory ${PWD}/${Sample}/variants/${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --algorithm xgboost  --dbsnp-vcf  /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${PWD}/${Sample}/gatk38_processing/${Sample}.final.bam --mutect2-vcf ${PWD}/${Sample}/variants/${Sample}.mutect2.vcf --vardict-vcf ${PWD}/${Sample}/variants/${Sample}.vardict.vcf --varscan-vcf ${PWD}/${Sample}/variants/${Sample}.varscan.vcf --lofreq-vcf ${PWD}/${Sample}/variants/${Sample}.lofreq.filtered.vcf --strelka-vcf ${PWD}/${Sample}/variants/${Sample}.strelka.vcf --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypus_cnvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf
	
	${params.vcf_sorter_path} ${PWD}/${Sample}/variants/${Sample}.somaticseq/Consensus.sSNV.vcf ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf > ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz
	
	${params.vcf_sorter_path} ${PWD}/${Sample}/variants/${Sample}.somaticseq/Consensus.sINDEL.vcf ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf > ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz
	
	${params.bcftools_path} concat -a ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_snv.vcf.gz ${PWD}/${Sample}/variants/${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=MVDLK01,Number=7,Type=Integer,Description="Calling decision of the 7 algorithms: MuTect, VarScan2, VarDict, LoFreq, Strelka, SnvCaller_0, SnvCaller_1">/##INFO=<ID=MVDLKFP,Number=7,Type=String,Description="Calling decision of the 7 algorithms: MuTect, VarScan2, VarDict, LoFreq, Strelka, Freebayes, Platypus">/g' ${Sample}.somaticseq.vcf

	sed -i 's/MVDLK01/MVDLKFP/g' ${Sample}.somaticseq.vcf

	mkdir -p "$PWD/${Sample}/PCGR"
	mkdir -p "$PWD/Final_Output/${Sample}/PCGR"
	mkdir -p "$PWD/${Sample}/CPSR"
	mkdir -p "$PWD/Final_Output/${Sample}/CPSR"
	${params.pcgr_script_path} ${Sample}.somaticseq.vcf $PWD/${Sample}/PCGR/ ${Sample} $PWD/${Sample}/CPSR/ ${params.ensemblid_path}
	${params.variant_call_path} $PWD/${Sample}/PCGR/${Sample}*.tiers.tsv ${Sample}.somaticseq.vcf $PWD/${Sample}/PCGR/${Sample}_output.tsv
	cp $PWD/${Sample}/PCGR/*.html $PWD/${Sample}/PCGR/${Sample}*.tiers.tsv $PWD/${Sample}/PCGR/${Sample}_output.tsv $PWD/Final_Output/${Sample}/PCGR/
	cp $PWD/${Sample}/CPSR/*.html $PWD/${Sample}/CPSR/${Sample}*.tiers.tsv $PWD/Final_Output/${Sample}/CPSR/

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf  --outfile ${Sample}.somaticseq.avinput --withzyg --includeinfo
	cp ${Sample}.somaticseq.vcf ${PWD}/Final_Output/${Sample}/	
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	"""
}

process platypus_run{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.platypus.vcf'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
	output:
		tuple val(Sample), file ("*.platypus.vcf")
	script:
	"""
	python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBams[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6 --regions=${params.bedfile}_regions.txt
	"""
}

process coverage {
	publishDir "$PWD/${Sample}/coverage/", mode: 'copy'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBams[0]} > ${Sample}.bed
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}.bed > ${Sample}.counts.bed
	${params.bedtools} coverage -counts -a ${params.flt3_bedfile}.bed -b ${Sample}.bed > ${Sample}_pindel.counts.bed	
	"""
}

process pindel {
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*pindel_SI.vcf'
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*_pindel.hg19_multianno.csv'

	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	export BAM_2_PINDEL_ADAPT=${params.pindel}/Adaptor.pm
	sh ${params.pindel_config_script} -s ${Sample}
	${params.pindel}/pindel -f ${params.genome} -i $PWD/config.txt -c chr13 -o ${Sample}_pindel
	${params.pindel}/pindel2vcf -r ${params.genome} -p ${Sample}_pindel_SI -R hg19 -d 07102019 -v ${Sample}_pindel_SI.vcf


	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}_pindel_SI.vcf --outfile ${Sample}_pindel.avinput --withzyg --includeinfo

	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}_pindel.avinput ${params.annovarLatest_path}/humandb/ -buildver hg19 -out ${Sample}_pindel --remove -protocol refGene,cytoBand,cosmic84 --operation g,r,f -nastring '.' --otherinfo --csvout --thread 10 --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
 	"""
}

process format_pindel {
	publishDir "$PWD/${Sample}/pindel/", mode: 'copy', pattern: '*_final.pindel.csv'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_final.pindel.csv'
	
	input:
		tuple val (Sample), file(bedfile), file (multianno)
	output:
		val Sample
	script:
	"""
	python3 ${params.format_pindel_script} ${PWD}/${Sample}/coverage/${Sample}_pindel.counts.bed ${PWD}/${Sample}/pindel/${Sample}_pindel.hg19_multianno.csv ${PWD}/${Sample}/pindel/${Sample}_final.pindel.csv
	"""
}

process cnvkit_run {
	publishDir "$PWD/${Sample}/cnvkit/", mode: 'copy'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai)
	output:
		val Sample
	script:
	"""
	#cnvkit.py batch ${finalBam} -r ${params.cnvkitRef} -m hybrid --drop-low-coverage --output-dir ${PWD}/${Sample}/cnvkit/ --diagram --scatter
	${params.cnvkit_path} ${finalBam} ${params.cnvkitRef} ${PWD}/${Sample}/cnvkit/
	/${params.gene_scatter}/custom_scatter_v3.py ${params.gene_scatter}/chr_list_all.txt ${PWD}/${Sample}/cnvkit/${Sample}.final.cnr ${PWD}/${Sample}/cnvkit/${Sample}.final.cns ${Sample}
	cp *gene_scatter.pdf $PWD/${Sample}/cnvkit/
	cp *gene_scatter.pdf $PWD/Final_Output/${Sample}/
	"""

}


process coverview_run {
	executor="local"
	publishDir "$PWD/${Sample}/Coverview/", mode: 'copy'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile}.bed -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	cp ${Sample}.coverview_regions.csv ${PWD}/Coverview/${Sample}.coverview_regions.csv
	"""
}

process coverview_report {
	errorStrategy 'ignore'
	executor="local"
	input:
		val (Sample)
	output:
		val Sample
	script:
	"""
	python3 ${params.coverview_report_path} ${PWD}/Coverview/ ${PWD}/Final_Output/
	"""
}

process combine_variants{
	publishDir "$PWD/${Sample}/variants/", mode: 'copy'
	publishDir "$PWD/${Sample}/variants/", mode: 'copy', pattern: '*.avinput'
	publishDir "$PWD/${Sample}/ANNOVAR/", mode: 'copy', pattern: '*.hg19_multianno.csv'
	
	input:
		tuple val (Sample), file(freebayesVcf), file(platypusVcf)
	output:
		tuple val(Sample), file ("*.combined.vcf"),  file ("*.hg19_multianno.csv")
	script:
	"""
	grep "^#" ${PWD}/${Sample}/variants/${Sample}.freebayes.vcf > ${Sample}.freebayes.sorted.vcf
	grep -v "^#" ${PWD}/${Sample}/variants/${Sample}.freebayes.vcf | sort -k1,1V -k2,2g >> ${Sample}.freebayes.sorted.vcf
	
	grep "^#" ${PWD}/${Sample}/variants/${Sample}.platypus.vcf > ${Sample}.platypus.sorted.vcf
	grep -v "^#" ${PWD}/${Sample}/variants/${Sample}.platypus.vcf | sort -k1,1V -k2,2g >> ${Sample}.platypus.sorted.vcf
	
	${params.java_path}/java -jar ${params.GATK38_path} -T CombineVariants -R ${params.genome} --variant ${Sample}.freebayes.sorted.vcf --variant ${Sample}.platypus.sorted.vcf -o ${Sample}.combined.vcf -genotypeMergeOptions UNIQUIFY
	
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.combined.vcf  --outfile ${Sample}.combined.avinput --withzyg --includeinfo
	
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.combined.avinput --out ${Sample}.combined --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	"""
}

process cava {
	publishDir "$PWD/${Sample}/CAVA/", mode: 'copy'
	
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno), file(combinedVcf)
	
	output:
		tuple val(Sample), file ("*.cava.csv")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i $PWD/${Sample}/variants/${Sample}.somaticseq.vcf -o ${Sample}.somaticseq
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i $PWD/${Sample}/variants/${Sample}.combined.vcf -o ${Sample}.combined
	python3 ${params.cava_script_path} ${Sample}.somaticseq.txt ${Sample}.combined.txt ${Sample}.cava.csv
	"""
}

process format_somaticseq_combined {
	input:
		tuple val (Sample), file(somaticseqVcf), file (multianno), file (combinedVcf),  file (hg19_multianno)
	output:
		val Sample
	script:
	"""
	python3 ${params.format_somaticseq_script} ${PWD}/${Sample}/ANNOVAR/${Sample}.somaticseq.hg19_multianno.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.somaticseq.csv
	python3 ${params.format_combined_script} ${PWD}/${Sample}/ANNOVAR/${Sample}.combined.hg19_multianno.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.combined.csv
	"""
}

process format_concat_combine_somaticseq {
	input:
		tuple val (Sample), file ("*")
	output:
		val Sample
	script:
	"""
	sed -i '1d' ${PWD}/${Sample}/Annovar_Modified/${Sample}.combined.csv
	sed -i '1d' ${PWD}/${Sample}/Annovar_Modified/${Sample}.somaticseq.csv
	#python3 ${params.format_concat_script} ${PWD}/${Sample}/Annovar_Modified/${Sample}.combined.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.somaticseq.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.concat.csv
	cp ${PWD}/${Sample}/Annovar_Modified/${Sample}.somaticseq.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.concat.csv
	python3 ${params.format_remove_artefact_script} ${PWD}/${Sample}/Annovar_Modified/${Sample}.concat.csv ${params.artefactFile} ${PWD}/${Sample}/Annovar_Modified/${Sample}.final.concat.csv ${PWD}/${Sample}/Annovar_Modified/${Sample}.artefacts.csv
	sed -i '1iChr,Start,End,Ref,Alt,Variant_Callers,FILTER,SOMATIC_FLAG,VariantCaller_Count,REF_COUNT,ALT_COUNT,VAF,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,Gene_full_name.refGene,Function_description.refGene,Disease_description.refGene,cosmic84,PopFreqMax,1000G_ALL,ExAC_ALL,CG46,ESP6500siv2_ALL,InterVar_automated' ${PWD}/${Sample}/Annovar_Modified/${Sample}.final.concat.csv
	sed -i '1iChr,Start,End,Ref,Alt,Variant_Callers,FILTER,SOMATIC_FLAG,VariantCaller_Count,REF_COUNT,ALT_COUNT,VAF,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,Gene_full_name.refGene,Function_description.refGene,Disease_description.refGene,cosmic84,PopFreqMax,1000G_ALL,ExAC_ALL,CG46,ESP6500siv2_ALL,InterVar_automated' ${PWD}/${Sample}/Annovar_Modified/${Sample}.artefacts.csv
	"""
}

process merge_csv {
	input:
		tuple val (Sample), file (cava_csv)
	output:
		val Sample
	script:
	"""
	sed -i 's/\t/,/g' ${PWD}/${Sample}/cnvkit/${Sample}.final.cnr
	python3 ${params.pharma_marker_script} ${Sample} ${PWD}/${Sample}/Annovar_Modified/ ${params.pharma_input_xlxs} ${PWD}/${Sample}/${Sample}_pharma.csv
	python3 ${params.merge_csvs_script} ${Sample} ${PWD}/${Sample}/Annovar_Modified/ ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ${PWD}/${Sample}/CAVA/ ${PWD}/${Sample}/Coverview/${Sample}.coverview_regions.csv ${PWD}/${Sample}/pindel/${Sample}_final.pindel.csv ${PWD}/${Sample}/cnvkit/${Sample}.final.cnr ${PWD}/${Sample}/${Sample}_pharma.csv

	${params.concat_script_path} ${Sample} ${PWD}/${Sample}/Annovar_Modified/ $PWD/${Sample}/PCGR/${Sample}*.tiers.tsv $PWD/${Sample}/CPSR/${Sample}*.tiers.tsv
	${params.pcgr_cpsr_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ${Sample}.final.concat_append.csv $PWD/${Sample}/PCGR/${Sample}*.tiers.tsv $PWD/${Sample}/CPSR/${Sample}*.tiers.tsv
	mv output_temp.xlsx ${PWD}/Final_Output/${Sample}/${Sample}.xlsx
	"""
}

process Final_Output {
	input:
		tuple val (Sample), file ("*")
	output:
		val Sample
	script:
	"""
	python3 ${params.coveragePlot_script} ${Sample} $PWD/${Sample}/coverage/${Sample}.counts.bed $PWD/${Sample}/coverage/
	cp ${PWD}/${Sample}/coverage/${Sample}.Low_Coverage.png ${PWD}/Final_Output/${Sample}/
	cp ${PWD}/${Sample}/cnvkit/${Sample}.final-scatter.png ${PWD}/${Sample}/cnvkit/${Sample}.final-diagram.pdf ${PWD}/Final_Output/${Sample}/
	"""
}

process remove_files{
	errorStrategy 'ignore'
	input:
		tuple val (Sample), file ("*")
	script:
	"""
	rm -rf ${PWD}/${Sample}/
	"""
}

workflow MIPS {

    Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	
	main:
	trimming_trimmomatic(samples_ch) | pair_assembly_pear | mapping_reads | sam_conversion
	minimap_getitd(samples_ch)
	RealignerTargetCreator(sam_conversion.out)
	IndelRealigner(RealignerTargetCreator.out.join(sam_conversion.out)) | BaseRecalibrator
	PrintReads(IndelRealigner.out.join(BaseRecalibrator.out)) | generatefinalbam
	hsmetrics_run(generatefinalbam.out)
	platypus_run(generatefinalbam.out)
	coverage(generatefinalbam.out)
	freebayes_run(generatefinalbam.out)
	mutect2_run(generatefinalbam.out)
	vardict_run(generatefinalbam.out)
	varscan_run(generatefinalbam.out)
	lofreq_run(generatefinalbam.out)
	strelka_run(generatefinalbam.out)
	somaticSeq_run(mutect2_run.out.join(vardict_run.out.join(varscan_run.out.join(lofreq_run.out.join(strelka_run.out)))))
	pindel(generatefinalbam.out)
	cnvkit_run(generatefinalbam.out)
	coverview_run(generatefinalbam.out)
	coverview_report(coverview_run.out.toList())
	combine_variants(freebayes_run.out.join(platypus_run.out))
	cava(somaticSeq_run.out.join(combine_variants.out))
	format_somaticseq_combined(somaticSeq_run.out.join(combine_variants.out))
	format_concat_combine_somaticseq(format_somaticseq_combined.out)
	format_pindel(pindel.out.join(coverage.out))
	merge_csv(format_concat_combine_somaticseq.out.join(cava.out.join(format_pindel.out)))
	Final_Output(coverage.out.join(cnvkit_run.out))
	remove_files(merge_csv.out.join(coverview_run.out.join(Final_Output.out)))
}
