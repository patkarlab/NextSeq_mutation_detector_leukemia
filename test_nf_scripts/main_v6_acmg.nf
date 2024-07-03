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

process trimming_trimmomatic {
	maxForks 10
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

process pair_assembly_pear {
	memory '7.0 GB'
	input:
		tuple val (Sample), file(trimmedFiles)
	output:
		tuple val (Sample), file("*assembled.fastq")
	script:
	"""
	sleep 5s
	${params.pear_path} -f ${trimmedFiles[0]} -r ${trimmedFiles[2]} -o ${Sample} -n 53 -j 25
	"""
}

process mapping_reads{
	maxForks 15
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
	input:
		tuple val (Sample), file(samFile)
	output:
		tuple val(Sample), file ("*.fxd_sorted.bam"), file ("*.fxd_sorted.bam.bai")

	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} FixMateInformation I= ${samFile} O= ${Sample}.fxd.sam VALIDATION_STRINGENCY=SILENT
	mv ${samFile} ${Sample}.fxd.sam
	${params.samtools} view -bT ${params.genome} ${Sample}.fxd.sam > ${Sample}.fxd.bam
	${params.samtools} sort ${Sample}.fxd.bam > ${Sample}.fxd_sorted.bam
	${params.samtools} index ${Sample}.fxd_sorted.bam > ${Sample}.fxd_sorted.bam.bai
	"""
}

process minimap_getitd {
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

process RealignerTargetCreator {
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
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	input:
		tuple val (Sample), file(alignedRecalibratedBam)
	output:
		tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai"), file ("*.old_final.bam"), file ("*.old_final.bam.bai")
		
	script:
	"""
	${params.bedtools} sort -i ${params.bedfile}.bed > sorted.bed

	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${params.sequences}/${Sample}*.bam --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --targets sorted.bed --tmpdir ./ > abra.log

	${params.samtools} sort ${params.sequences}/${Sample}*.bam > ${Sample}.old_final.bam
	${params.samtools} index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	${params.samtools} sort ${Sample}.abra.bam > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

process hsmetrics_run{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_hsmetrics.txt'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	${params.hsmetrics_all} $PWD/Final_Output/hsmetrics.tsv ${Sample} ${Sample}_hsmetrics.txt
	"""
}

process mutect2_run{
	maxForks 10	
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.mutect2.vcf")
	script:
	"""
	#${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${finalBam} -o ${Sample}.mutect2.vcf --dbsnp ${params.site2} -L ${params.bedfile}.bed -nct 25 -contamination 0.02 -mbq 30
	${params.samtools} view -bs 40.1 ${finalBam} > subsampled_01.bam
	${params.samtools} index subsampled_01.bam
	${params.mutect2} ${params.java_path} ${params.GATK38_path} ${params.genome} subsampled_01.bam ${Sample}.mutect2.vcf ${params.site2} ${params.bedfile}.bed 
	"""
}

process freebayes_run{	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.freebayes.vcf")
	script:
	"""
	${params.freebayes_path} -f ${params.genome} -b ${finalBam} -t ${params.bedfile}.bed > ${Sample}.freebayes.vcf 	
	"""
}

process vardict_run{	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${finalBam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf
	"""
}

process varscan_run{	
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file("*.varscan.vcf")
		
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

process strelka_run{
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

process somaticSeq_run {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.somaticseq.vcf'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file(oldfinalBam), file(oldfinalBamBai), file(mutectVcf), file(vardictVcf), file(varscanVcf), file(lofreqVcf), file(strelkaVcf), file(freebayesVcf), file(platypusVcf)
	output:
		tuple val (Sample), file("*.somaticseq.vcf"), file("*.hg19_multianno.csv"), file("*.hg19_multianno.txt.cancervar.ensemble.pred")
	script:
	"""
	${params.vcf_sorter_path} ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${platypusVcf} ${Sample}.platypus.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf

	${params.vcf_sorter_path} ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf

	somaticseq_parallel.py --output-directory ./${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --algorithm xgboost  --dbsnp-vcf  /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${mutectVcf} --vardict-vcf ${vardictVcf} --varscan-vcf ${varscanVcf} --lofreq-vcf ${lofreqVcf} --strelka-vcf ${strelkaVcf} --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypus_cnvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf 

	${params.vcf_sorter_path} ./${Sample}.somaticseq/Consensus.sSNV.vcf ./${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ./${Sample}.somaticseq/somaticseq_snv.vcf > ./${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ./${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ./${Sample}.somaticseq/Consensus.sINDEL.vcf ./${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ./${Sample}.somaticseq/somaticseq_indel.vcf > ./${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ./${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ./${Sample}.somaticseq/somaticseq_snv.vcf.gz ./${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ./${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=MVDLK01,Number=7,Type=Integer,Description="Calling decision of the 7 algorithms: MuTect, VarScan2, VarDict, LoFreq, Strelka, SnvCaller_0, SnvCaller_1">/##INFO=<ID=MVDLKFP,Number=7,Type=String,Description="Calling decision of the 7 algorithms: MuTect, VarScan2, VarDict, LoFreq, Strelka, Freebayes, Platypus">/g' ${Sample}.somaticseq.vcf
	sed -i 's/MVDLK01/MVDLKFP/g' ${Sample}.somaticseq.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf  --outfile ${Sample}.somaticseq.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	${params.cancervar} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}
	"""
}

process CombineVariants {
	input:
		tuple val (Sample), file(mutectVcf), file(vardictVcf), file(varscanVcf), file(lofreqVcf), file(strelkaVcf), file(freebayesVcf), file(platypusVcf)
	output:
		tuple val (Sample), file("${Sample}.combined.vcf")
	script:
	"""
	${params.vcf_sorter_path} ${mutectVcf} ${Sample}.mutect.sorted.vcf
	${params.vcf_sorter_path} ${vardictVcf} ${Sample}.vardict.sorted.vcf
	${params.vcf_sorter_path} ${varscanVcf} ${Sample}.varscan.sorted.vcf
	${params.vcf_sorter_path} ${lofreqVcf} ${Sample}.lofreq.sorted.vcf
	${params.vcf_sorter_path} ${strelkaVcf} ${Sample}.strelka.sorted.vcf
	${params.vcf_sorter_path} ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${platypusVcf} ${Sample}.platypus.sorted.vcf
	
	${params.java_path}/java -jar ${params.GATK38_path} -T CombineVariants -R ${params.genome} --variant ${Sample}.mutect.sorted.vcf --variant ${Sample}.vardict.sorted.vcf \
	--variant ${Sample}.varscan.sorted.vcf --variant ${Sample}.lofreq.sorted.vcf --variant ${Sample}.strelka.sorted.vcf --variant ${Sample}.freebayes.sorted.vcf --variant ${Sample}.platypus.sorted.vcf \
	-o ${Sample}.combined.vcf -genotypeMergeOptions UNIQUIFY
	"""
}

process platypus_run{
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file ("*.platypus.vcf")
	script:
	"""
	python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBams[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --maxVariants=6 --minReads=6 --regions=${params.bedfile}_regions.txt
	"""
}

process coverage {
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed"), file ("${Sample}_pindel.counts.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBams[0]} > ${Sample}.bed
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}.bed > ${Sample}.counts.bed
	${params.bedtools} coverage -counts -a ${params.flt3_bedfile}.bed -b ${Sample}.bed > ${Sample}_pindel.counts.bed	
	"""
}

process pindel {
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

process format_pindel {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_final.pindel.csv'
	input:
		tuple val (Sample), file (pindelMultianno), file (countsBed), file (pindelCountsBed)
	output:
		tuple val (Sample), file("*_final.pindel.csv")
	script:
	"""
	python3 ${params.format_pindel_script} ${pindelCountsBed} ${pindelMultianno} ${Sample}_final.pindel.csv
	"""
}

process cnvkit_run {
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
	${params.cnvkit_path} ${finalBam} ${params.cnvkitRef} ./
	/${params.gene_scatter}/custom_scatter_v3.py ${params.gene_scatter}/chr_list_all.txt ./${Sample}.final.cnr ./${Sample}.final.cns ${Sample}
	/${params.gene_scatter}/custom_scatter_chrwise.py ${params.gene_scatter}/chrwise_list.txt ./${Sample}.final.cnr ./${Sample}.final.cns ${Sample}_chr_
	"""
}

process ifcnv_run {
	input:
		val Sample
	output:
		//tuple val (Sample), file ("*.html"), file ("*.tsv")
		val (Sample)
	script:
	"""
	${params.links} $PWD/Final_Output/ ${params.input}
	mkdir ifCNV
	${params.ifcnv} ./ ${params.bedfile}.bed ifCNV

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

process update_db {
	publishDir "$PWD/work/", mode: 'copy', pattern: 'freq.txt'
	input:
		val Sample
	output:
		file ("freq.txt")
	script:
	"""
	for i in `cat ${params.input}`
	do 
		if [ -f ${PWD}/Final_Output/\${i}/\${i}.somaticseq.vcf ]; then
			ln -s ${PWD}/Final_Output/\${i}/\${i}.somaticseq.vcf ./
		fi
	done
	files=\$(ls *.somaticseq.vcf)
	${params.updatedb} ${params.alpdb} \${files}
	${params.freq_py} ${params.alpdb} freq.txt
	"""
}

process annotSV {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_AnnotSV.tsv'
	input:
		tuple val (Sample), file (finalCns), file (finalCnr), file (geneScatter), file (finalScatter), file (finalDiagram)
	output:
		tuple val (Sample), file ("*_AnnotSV.tsv")
	script:
	"""
	/${params.annotsv} ${finalCns} ${Sample}
	/${params.substitute_null} ${Sample}_annotsv.tsv ${Sample}_AnnotSV.tsv
	"""
}

process igv_reports {
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno), file (cancervarMultianno)		
	output:
		val(Sample)
	script:
	"""
	perl ${params.annovarLatest_path}/table_annovar.pl ${somaticVcf} --out ${Sample}.annovar --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring . --otherinfo --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt -vcfinput 

	${params.igv_script} ${params.genome} ${Sample}.annovar.hg19_multianno.vcf $PWD/Final_Output/${Sample}/${Sample}.final.bam $PWD/Final_Output/${Sample}/${Sample}_igv.html
	"""		
}

process coverview_run {
	executor="local"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: "*.coverview_regions.csv"
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("*.coverview_regions.csv")
	script:
	"""
	${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile}.bed -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	"""
}

process cava {
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno), file(cancervarMultianno), file(combineVcf)
	output:
		tuple val(Sample), file ("*.cava.csv")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${somaticVcf} -o ${Sample}.somaticseq
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${combineVcf} -o ${Sample}.combined
	python3 ${params.cava_script_path} ${Sample}.somaticseq.txt ${Sample}.combined.txt ${Sample}.cava.csv
	"""
}

process format_somaticseq_combined {
	input:
		tuple val (Sample), file(somaticseqVcf), file (multianno), file (cancervarMultianno)
	output:
		tuple val (Sample), file("*.somaticseq.csv")
	script:
	"""
	python3 ${params.format_somaticseq_script} ${multianno} ${Sample}.somaticseq.csv
	"""
}

process format_concat_combine_somaticseq {
	input:
		tuple val (Sample), file (somaticseqCsv)
	output:
		tuple val (Sample), file ("${Sample}.final.concat.csv"), file ("${Sample}.artefacts.csv")
	script:
	"""
	sed -i '1d' ${somaticseqCsv}
	python3 ${params.format_remove_artefact_script} ${somaticseqCsv} ${params.artefactFile} ./${Sample}.final.concat.csv ./${Sample}.artefacts.csv
	sed -i '1iChr,Start,End,Ref,Alt,Variant_Callers,FILTER,SOMATIC_FLAG,VariantCaller_Count,REF_COUNT,ALT_COUNT,VAF,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,Gene_full_name.refGene,Function_description.refGene,Disease_description.refGene,cosmic84,PopFreqMax,1000G_ALL,ExAC_ALL,CG46,ESP6500siv2_ALL,InterVar_automated' ${Sample}.final.concat.csv
	sed -i '1iChr,Start,End,Ref,Alt,Variant_Callers,FILTER,SOMATIC_FLAG,VariantCaller_Count,REF_COUNT,ALT_COUNT,VAF,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,Gene_full_name.refGene,Function_description.refGene,Disease_description.refGene,cosmic84,PopFreqMax,1000G_ALL,ExAC_ALL,CG46,ESP6500siv2_ALL,InterVar_automated' ${Sample}.artefacts.csv
	"""
}

process merge_csv {
	input:
		tuple val (Sample), file (finalConcat), file (artefacts), file (cavaCsv), file (coverviewRegions) ,file (finalPindel), file(finalCns), file(finalCnr), file(geneScatter), file (finalScatter), file (finalDiagram), file (somaticVcf), file (somaticseqMultianno), file (cancervarMultianno)
	output:
		val Sample
	script:
	"""
	sed -i 's/\t/,/g' ${finalCnr}
	python3 ${params.pharma_marker_script} ${Sample} ./ ${params.pharma_input_xlxs} ./${Sample}_pharma.csv
	python3 ${params.merge_csvs_script} ${Sample} ./ ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ./ ${coverviewRegions} ${finalPindel} ${finalCnr} ./${Sample}_pharma.csv

	cp ${finalConcat} ${Sample}.final.concat_append.csv
	${params.vep_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.somaticseq.vcf ${PWD}/Final_Output/${Sample}/${Sample}
	${params.vep_extract_path} ${Sample}.final.concat_append.csv ${PWD}/Final_Output/${Sample}/${Sample}_vep_delheaders.txt > ${Sample}.vep
	${params.cancervar_extract} ${cancervarMultianno} ${Sample}.vep ${Sample}_cancervar.csv
	${params.pcgr_cpsr_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ${Sample}_cancervar.csv
	mv output_temp.xlsx ${PWD}/Final_Output/${Sample}/${Sample}.xlsx
	"""
}

process update_freq {
	input:
		val (Sample)
	output:
		val Sample
	script:
	"""
	ln -s ${PWD}/work/freq.txt ./
	for i in `cat ${params.input}`
	do
		if [ -f ${PWD}/Final_Output/\${i}/\${i}.xlsx ]; then
			${params.update_freq_excel} ${PWD}/Final_Output/\${i}/\${i}.xlsx freq.txt
		fi
	done
	"""  	
}

process Final_Output {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.Low_Coverage.png'
	input:
		tuple val (Sample), file(countsBed), file(pindelCountsBed), file(finalCns), file(finalCnr), file(geneScatter), file (finalScatter), file (finalDiagram)
	output:
		tuple val (Sample), file("*.Low_Coverage.png")
	script:
	"""
	python3 ${params.coveragePlot_script} ${Sample} ${countsBed} ./
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
	somaticSeq_run(generatefinalbam.out.join(mutect2_run.out.join(vardict_run.out.join(varscan_run.out.join(lofreq_run.out.join(strelka_run.out.join(freebayes_run.out.join(platypus_run.out))))))))
	CombineVariants(mutect2_run.out.join(vardict_run.out.join(varscan_run.out.join(lofreq_run.out.join(strelka_run.out.join(freebayes_run.out.join(platypus_run.out)))))))
	pindel(generatefinalbam.out)
	cnvkit_run(generatefinalbam.out)
	annotSV(cnvkit_run.out)
	ifcnv_run(generatefinalbam.out.collect())
	igv_reports(somaticSeq_run.out)
	update_db(somaticSeq_run.out.collect())
	coverview_run(generatefinalbam.out)
	cava(somaticSeq_run.out.join(CombineVariants.out))
	format_somaticseq_combined(somaticSeq_run.out)
	format_concat_combine_somaticseq(format_somaticseq_combined.out)
	format_pindel(pindel.out.join(coverage.out))
	merge_csv(format_concat_combine_somaticseq.out.join(cava.out.join(coverview_run.out.join(format_pindel.out.join(cnvkit_run.out.join(somaticSeq_run.out))))))
	update_freq(merge_csv.out.collect())
	Final_Output(coverage.out.join(cnvkit_run.out))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
