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


process filt3r {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_filt3r_out.csv'

	input:
		val (Sample)

	output:
		tuple val (Sample), file("*_filt3r_out.csv")

	script:
	"""
	${params.bedtools} bamtofastq -i ${params.sequences}/${Sample}_tumor.bam -fq ${Sample}_R1.fastq -fq2 ${Sample}_R2.fastq
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

process getitd {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_getitd'
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

process generatefinalbam{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.bam.bai'
	input:
		val Sample
	output:
		tuple val(Sample), file ("*.final.bam"), file ("*.final.bam.bai"), file ("*.old_final.bam"), file ("*.old_final.bam.bai")
		
	script:
	"""
	${params.bedtools} sort -i ${params.bedfile}.bed > sorted.bed

	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${params.sequences}/${Sample}_tumor.bam --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --targets sorted.bed --tmpdir ./ > abra.log

	${params.samtools} sort ${params.sequences}/${Sample}_tumor.bam > ${Sample}.old_final.bam
	${params.samtools} index ${Sample}.old_final.bam > ${Sample}.old_final.bam.bai
	${params.samtools} sort ${Sample}.abra.bam > ${Sample}.final.bam
	${params.samtools} index ${Sample}.final.bam > ${Sample}.final.bam.bai
	"""
}

process hsmetrics_run{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_hsmetrics*.txt'
	input:
		tuple val(Sample), file(finalBam), file (finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		path ("${Sample}_hsmetrics.txt"), emit : genewise
		path ("${Sample}_hsmetrics_exonwise.txt"), emit : exonwise
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT

	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${finalBam} O= ${Sample}_hsmetrics_exonwise.txt BAIT_INTERVALS= ${params.bedfile_exonwise}.interval_list TARGET_INTERVALS= ${params.bedfile_exonwise}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT

	#${params.hsmetrics_all} $PWD/Final_Output/hsmetrics.tsv ${Sample} ${Sample}_hsmetrics.txt
	"""
}

process hsmetrics_collect {
	publishDir "$PWD/Final_Output/", mode: 'copy'
	input:
		file (GeneWise) 
		file (ExonWise)
	output:
		file("hsmetrics_probewise.txt") 
		file("hsmetrics_exonwise.txt")
	script:
	"""
	echo -e "Sample name\tOn target\tOff target" > hsmetrics_probewise.txt
	for i in ${GeneWise}
	do
		samp_name=\$(basename -s .txt \${i})
		grep -v '#' \${i} | awk -v name=\${samp_name} 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print name,\$7,\$8}' >> hsmetrics_probewise.txt
	done

	echo -e "Sample name\tOn target\tOff target" > hsmetrics_exonwise.txt
	for i in ${ExonWise}
	do
		samp_name=\$(basename -s .txt \${i})
		grep -v '#' \${i} | awk -v name=\${samp_name} 'BEGIN{FS="\t"; OFS="\t"}NR==3{ print name,\$7,\$8}' >> hsmetrics_exonwise.txt
	done
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
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${finalBam} -O 50 -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf
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

process DeepSomatic{
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

process somaticSeqDragen_run {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.somaticseq.vcf'
	input:
		tuple val (Sample), file(finalBam), file(finalBamBai), file(oldfinalBam), file(oldfinalBamBai) , file(mutectVcf), file(vardictVcf), file(DeepSomaticVcf), file(lofreqVcf), file(strelkaVcf), file(freebayesVcf), file(platypusVcf)
	output:
		tuple val (Sample), file("*.somaticseq.vcf"), file("*.hg19_multianno.csv"), file("*.hg19_multianno.txt.cancervar.ensemble.pred")
	script:
	"""
	${params.vcf_sorter_path} ${freebayesVcf} ${Sample}.freebayes.sorted.vcf
	${params.vcf_sorter_path} ${platypusVcf} ${Sample}.platypus.sorted.vcf
	${params.vcf_sorter_path} ${params.sequences}/${Sample}.hard-filtered.vcf ${Sample}.dragen.sorted.vcf
	${params.vcf_sorter_path} ${DeepSomaticVcf} ${Sample}.deepsomatic.sorted.vcf

	python3 ${params.splitvcf_path} -infile ${Sample}.platypus.sorted.vcf -snv ${Sample}_platypus_cnvs.vcf -indel ${Sample}_platypus_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.freebayes.sorted.vcf -snv ${Sample}_freebayes_cnvs.vcf -indel ${Sample}_freebayes_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.dragen.sorted.vcf -snv ${Sample}_dragen_cnvs.vcf -indel ${Sample}_dragen_indels.vcf
	python3 ${params.splitvcf_path} -infile ${Sample}.deepsomatic.sorted.vcf -snv ${Sample}_deepsomatic_snvs.vcf -indel ${Sample}_deepsomatic_indels.vcf

	${params.vcf_sorter_path} ${Sample}_platypus_cnvs.vcf ${Sample}_platypus_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_platypus_indels.vcf ${Sample}_platypus_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_cnvs.vcf ${Sample}_freebayes_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_freebayes_indels.vcf ${Sample}_freebayes_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_dragen_cnvs.vcf ${Sample}_dragen_cnvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_dragen_indels.vcf ${Sample}_dragen_indels_sort.vcf
	${params.vcf_sorter_path} ${Sample}_deepsomatic_snvs.vcf ${Sample}_deepsomatic_snvs_sort.vcf
	${params.vcf_sorter_path} ${Sample}_deepsomatic_indels.vcf ${Sample}_deepsomatic_indels_sort.vcf

	somaticseq_parallel.py --output-directory ./${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --algorithm xgboost  --dbsnp-vcf  /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${mutectVcf} --vardict-vcf ${vardictVcf} --lofreq-vcf ${lofreqVcf} --strelka-vcf ${strelkaVcf} --sample-name ${Sample} --arbitrary-snvs ${Sample}_freebayes_cnvs_sort.vcf ${Sample}_platypus_cnvs_sort.vcf ${Sample}_dragen_cnvs_sort.vcf ${Sample}_deepsomatic_snvs_sort.vcf --arbitrary-indels ${Sample}_freebayes_indels_sort.vcf ${Sample}_platypus_indels_sort.vcf ${Sample}_dragen_indels_sort.vcf ${Sample}_deepsomatic_indels_sort.vcf
	
	${params.vcf_sorter_path} ./${Sample}.somaticseq/Consensus.sSNV.vcf ./${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ./${Sample}.somaticseq/somaticseq_snv.vcf > ./${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ./${Sample}.somaticseq/somaticseq_snv.vcf.gz
	
	${params.vcf_sorter_path} ./${Sample}.somaticseq/Consensus.sINDEL.vcf ./${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ./${Sample}.somaticseq/somaticseq_indel.vcf > ./${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ./${Sample}.somaticseq/somaticseq_indel.vcf.gz
	
	${params.bcftools_path} concat -a ./${Sample}.somaticseq/somaticseq_snv.vcf.gz ./${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ./${Sample}.somaticseq.vcf

	sed -i 's/##INFO=<ID=MDLK0123,Number=8,Type=Integer,Description="Calling decision of the 8 algorithms: MuTect, VarDict, LoFreq, Strelka, SnvCaller_0, SnvCaller_1, SnvCaller_2, SnvCaller_3">/##INFO=<ID=MDLKFPGS,Number=8,Type=String,Description="Calling decision of the 8 algorithms: MuTect, VarDict, LoFreq, Strelka, Freebayes, Platypus, Dragen, DeepSomatic">/g' ${Sample}.somaticseq.vcf

	sed -i 's/MDLK0123/MDLKFPGS/g' ${Sample}.somaticseq.vcf
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf  --outfile ${Sample}.somaticseq.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	${params.cancervar} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}

	"""
}

process mocha {
	publishDir "$PWD/Final_Output/${Sample}/MoChA/", mode: 'copy', pattern: '*.png'
	publishDir "$PWD/Final_Output/${Sample}/MoChA/", mode: 'copy', pattern: '*.pdf'
	publishDir "$PWD/Final_Output/${Sample}/MoChA/", mode: 'copy', pattern: '*.tsv'
	input:
		tuple val (Sample), file(somaticseqVcf), file (multianno), file (cancervarMultianno)
	output:
		tuple val(Sample), file ("*.png"), file ("*.pdf"), file ("*.tsv")
	script:
	"""
	mv ${somaticseqVcf} ${Sample}.somaticseq_old.vcf
	${params.bedtools} intersect -a ${Sample}.somaticseq_old.vcf -b ${params.mocha_bedfile} -header > ${Sample}.somaticseq.vcf
	${params.mocha} ${Sample} ./
	"""
}


process platypus_run{
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val(Sample), file ("*.platypus.vcf")
	script:
	"""
	python2.7 ${params.platypus_path} callVariants --bamFiles=${finalBams[0]} --refFile=${params.genome} --output=${Sample}.platypus.vcf --nCPU=15 --minFlank=0 --filterDuplicates=0 --minMapQual=50 --maxVariants=6 --minReads=6 --regions=${params.bedfile}_regions.txt
	"""
}

process coverage {
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai), file (oldfinalBam), file (oldfinalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.counts.bed"), file ("${Sample}_pindel.counts.bed"), file("${Sample}_pindel_ubtf.counts.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${finalBams[0]} > ${Sample}.bed
	${params.bedtools} coverage -counts -a ${params.bedfile_exonwise}.bed -b ${Sample}.bed > ${Sample}.counts.bed
	${params.bedtools} coverage -counts -a ${params.flt3_bedfile}.bed -b ${Sample}.bed > ${Sample}_pindel.counts.bed
	grep 'UBTF' ${params.bedfile}.bed > ubtf_pindel.bed
	${params.bedtools} coverage -counts -a ubtf_pindel.bed -b ${Sample}.bed > ${Sample}_pindel_ubtf.counts.bed
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

process pindel_UBTF {
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

process format_pindel {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_final.pindel.csv'
	input:
		tuple val (Sample), file (pindelMultianno), file (countsBed), file (pindelCountsBed), file (pindelUBTFCountsBed)
	output:
		tuple val (Sample), file("*_final.pindel.csv")
	script:
	"""
	python3 ${params.format_pindel_script} ${pindelCountsBed} ${pindelMultianno} ${Sample}_final.pindel.csv
	"""
}

process format_pindel_UBTF {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_final.pindel_ubtf.csv'
	input:
		tuple val (Sample), file (pindelMultianno), file (countsBed), file (pindelCountsBed), file (pindelUBTFCountsBed)
	output:
		tuple val (Sample), file("*_final.pindel_ubtf.csv")
	script:
	"""
	python3 ${params.format_pindel_script} ${pindelUBTFCountsBed} ${pindelMultianno} ${Sample}_final.pindel_ubtf.csv
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
	# Commenting the following command for CNV myeloid panel, uncomment for the usual AL panel
	#/${params.gene_scatter}/custom_scatter_v3.py ${params.gene_scatter}/chr_list_all.txt ./${Sample}.final.cnr ./${Sample}.final.cns ${Sample}

	/${params.gene_scatter}/custom_scatter_chrwise.py ${params.gene_scatter_list}/chrwise_list.txt ./${Sample}.final.cnr ./${Sample}.final.cns ${Sample}_chr_
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
	${params.coverview_path}/coverview -i ${finalBam} -b ${params.bedfile_exonwise}.bed -c ${params.coverview_path}/config/config.txt -o ${Sample}.coverview
	python3 ${params.coverview_script_path} ${Sample}.coverview_regions.txt ${Sample}.coverview_regions.csv
	"""
}

process combine_variants{
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
	${params.vcf_sorter_path} ${params.sequences}/${Sample}.hard-filtered.vcf ${Sample}.dragen.sorted.vcf

	${params.java_path}/java -jar ${params.GATK38_path} -T CombineVariants -R ${params.genome} --variant ${Sample}.mutect.sorted.vcf --variant ${Sample}.vardict.sorted.vcf --variant ${Sample}.DeepSomatic.sorted.vcf --variant ${Sample}.lofreq.sorted.vcf --variant ${Sample}.strelka.sorted.vcf --variant ${Sample}.freebayes.sorted.vcf --variant ${Sample}.platypus.sorted.vcf --variant ${Sample}.dragen.sorted.vcf -o ${Sample}.combined.vcf -genotypeMergeOptions UNIQUIFY
	"""
}


process cava {
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno), file(cancervarMultianno), file(combinedVcf)
	
	output:
		tuple val(Sample), file ("*.cava.csv")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${somaticVcf} -o ${Sample}.somaticseq
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${combinedVcf} -o ${Sample}.combined
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
	#python3 ${params.format_somaticseq_script} ${multianno} ${Sample}.somaticseq.csv
	python3 /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/somaticseqoutput-format_dragen.py ${multianno} ${Sample}.somaticseq.csv
	"""
}

process format_concat_combine_somaticseq {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.concat.csv'
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
		tuple val (Sample), file (finalConcat), file (artefacts), file (cavaCsv), file (coverviewRegions) ,file (finalPindel), file(finalCns), file(finalCnr), file(geneScatter), file (finalScatter), file (finalDiagram), file (somaticVcf), file (somaticseqMultianno), file (cancervarMultianno), file(avinput), file(filt3rCsv), file(finalPindelUBTF)
	output:
		val Sample
	script:
	"""
	sed -i 's/\t/,/g' ${finalCnr}
	python3 ${params.pharma_marker_script} ${Sample} ./ ${params.pharma_input_xlxs} ./${Sample}_pharma.csv
	python3 ${params.merge_csvs_script} ${Sample} ./ ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ./ ${coverviewRegions} ${finalPindel} ${finalCnr} ./${Sample}_pharma.csv ${finalPindelUBTF}

	cp ${finalConcat} ${Sample}.final.concat_append.csv
	${params.vep_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.somaticseq.vcf ${PWD}/Final_Output/${Sample}/${Sample}
	${params.vep_extract_path} ${Sample}.final.concat_append.csv ${PWD}/Final_Output/${Sample}/${Sample}_vep_delheaders.txt > ${Sample}.vep
	${params.cancervar_extract} ${cancervarMultianno} ${Sample}.vep ${Sample}_cancervar.csv
	${params.pcgr_cpsr_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ${Sample}_cancervar.csv ${filt3rCsv}

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
		tuple val (Sample), file(countsBed), file(pindelCountsBed), file(pindelUBTFCountsBed), file(finalCns), file(finalCnr), file(geneScatter), file (finalScatter), file (finalDiagram)
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
	generatefinalbam(samples_ch)
	getitd(generatefinalbam.out)
	hsmetrics_run(generatefinalbam.out)
	platypus_run(generatefinalbam.out)
	coverage(generatefinalbam.out)
	freebayes_run(generatefinalbam.out)
	mutect2_run(generatefinalbam.out)
	vardict_run(generatefinalbam.out)
	//varscan_run(generatefinalbam.out)
	DeepSomatic(generatefinalbam.out)
	lofreq_run(generatefinalbam.out)
	strelka_run(generatefinalbam.out)
	somaticSeqDragen_run(generatefinalbam.out.join(mutect2_run.out.join(vardict_run.out.join(DeepSomatic.out.join(lofreq_run.out.join(strelka_run.out.join(freebayes_run.out.join(platypus_run.out))))))))
	combine_variants(mutect2_run.out.join(vardict_run.out.join(DeepSomatic.out.join(lofreq_run.out.join(strelka_run.out.join(freebayes_run.out.join(platypus_run.out)))))))
	pindel(generatefinalbam.out)
	cnvkit_run(generatefinalbam.out)
	//annotSV(cnvkit_run.out)
	//ifcnv_run(generatefinalbam.out.collect())
	igv_reports(somaticSeqDragen_run.out)
	//update_db(somaticSeqDragen_run.out.collect())
	coverview_run(generatefinalbam.out)
	cava(somaticSeqDragen_run.out.join(combine_variants.out))
	format_somaticseq_combined(somaticSeqDragen_run.out)
	format_concat_combine_somaticseq(format_somaticseq_combined.out)
	format_pindel(pindel.out.join(coverage.out))
	merge_csv(format_concat_combine_somaticseq.out.join(cava.out.join(coverview_run.out.join(format_pindel.out.join(cnvkit_run.out.join(somaticSeqDragen_run.out))))))
	//update_freq(merge_csv.out.collect())
	Final_Output(coverage.out.join(cnvkit_run.out))
}

workflow MIPS_mocha {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }
	main:
	//filt3r(samples_ch)
	generatefinalbam(samples_ch)
	//getitd(generatefinalbam.out)
	//hsmetrics_run(generatefinalbam.out)

	//all_hsmetrics_gw = hsmetrics_run.out.genewise.collect()
	//all_hsmetrics_pw = hsmetrics_run.out.exonwise.collect()

	//hsmetrics_collect(all_hsmetrics_gw, all_hsmetrics_pw)

	//platypus_run(generatefinalbam.out)
	//coverage(generatefinalbam.out)
	//freebayes_run(generatefinalbam.out)
	//mutect2_run(generatefinalbam.out)
	//vardict_run(generatefinalbam.out)
	////varscan_run(generatefinalbam.out)
	//DeepSomatic(generatefinalbam.out)
	//lofreq_run(generatefinalbam.out)
	//strelka_run(generatefinalbam.out)
	//somaticSeqDragen_run(generatefinalbam.out.join(mutect2_run.out.join(vardict_run.out.join(DeepSomatic.out.join(lofreq_run.out.join(strelka_run.out.join(freebayes_run.out.join(platypus_run.out))))))))
	////mocha(somaticSeqDragen_run.out)
	//combine_variants(mutect2_run.out.join(vardict_run.out.join(DeepSomatic.out.join(lofreq_run.out.join(strelka_run.out.join(freebayes_run.out.join(platypus_run.out)))))))
	//pindel(generatefinalbam.out)
	//pindel_UBTF(generatefinalbam.out)
	cnvkit_run(generatefinalbam.out)
	//annotSV(cnvkit_run.out)
	//ifcnv_run(generatefinalbam.out.collect())
	//igv_reports(somaticSeqDragen_run.out)
	//update_db(somaticSeqDragen_run.out.collect())
	//coverview_run(generatefinalbam.out)
	//cava(somaticSeqDragen_run.out.join(combine_variants.out))
	//format_somaticseq_combined(somaticSeqDragen_run.out)
	//format_concat_combine_somaticseq(format_somaticseq_combined.out)
	//format_pindel(pindel.out.join(coverage.out))
	//format_pindel_UBTF(pindel_UBTF.out.join(coverage.out))
	//merge_csv(format_concat_combine_somaticseq.out.join(cava.out.join(coverview_run.out.join(format_pindel.out.join(cnvkit_run.out.join(somaticSeqDragen_run.out.join(filt3r.out.join(format_pindel_UBTF.out))))))))
	//update_freq(merge_csv.out.collect())
	//Final_Output(coverage.out.join(cnvkit_run.out))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
