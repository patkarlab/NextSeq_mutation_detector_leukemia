#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

// Adapter Trimming, alignment and GATK BQSR (based on https://github.com/GavinHaLab/fastq_to_bam_paired_snakemake)
include { FASTQTOBAM; ABRA_BAM } from './scripts/processes.nf'

// FLT3 ITD detection
include { FILT3R; GETITD } from './scripts/flt3_itd.nf'

// HSmetrics calculation
include { HSMETRICS; HSMETRICS_COLLECT } from './scripts/hsmetrics.nf'

// COVERAGE calculation
include { COVERAGE; COVERVIEW } from './scripts/coverage.nf'

// Variant calling
include { PLATYPUS; FREEBAYES; MUTECT2; VARDICT; DEEPSOMATIC; LOFREQ; STRELKA; PINDEL; PINDEL_UBTF} from './scripts/variant_call.nf'
// Variant integration 
include { SOMATICSEQ; COMBINE_VARIANTS } from './scripts/somaticseq.nf'

// CNV calling
include { CNVKIT; ANNOT_SV; IFCNV } from './scripts/cnv_call.nf'

// IGV reports
include { IGV_REPORTS } from './scripts/igv_reports.nf'

workflow MyoPool {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		// .map { sample -> def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz") 
		// 				def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz")
		// 				tuple(sample, r1, r2)}
		.map { sample -> def r1 = file("${params.sequences}/${sample}_R1.fastq.gz")
						def r2 = file("${params.sequences}/${sample}_R2.fastq.gz")
						tuple(sample, r1, r2)}
		.set { bam_ch }

	main:
	// Adapter Trimming, alignment and GATK BQSR
	final_bams_ch = FASTQTOBAM(bam_ch)
	ABRA_BAM (final_bams_ch)

	// FLT3 ITD detection
	FILT3R(ABRA_BAM.out)
	GETITD(ABRA_BAM.out)

	// HSmetrics calculation 
	HSMETRICS(ABRA_BAM.out)
	all_hsmetrics_gw = HSMETRICS.out.genewise.collect()
	all_hsmetrics_pw = HSMETRICS.out.exonwise.collect()
	HSMETRICS_COLLECT(all_hsmetrics_gw, all_hsmetrics_pw)

	// COVERAGE calculation
	COVERAGE(ABRA_BAM.out)
	COVERVIEW(ABRA_BAM.out)

	// Variant calling 
	PLATYPUS(ABRA_BAM.out)
	FREEBAYES(ABRA_BAM.out)
	MUTECT2(ABRA_BAM.out)
	VARDICT(ABRA_BAM.out)
	DEEPSOMATIC(ABRA_BAM.out)
	LOFREQ(ABRA_BAM.out)
	STRELKA(ABRA_BAM.out)
	PINDEL(ABRA_BAM.out)
	PINDEL_UBTF(ABRA_BAM.out)

	// Variant integration 
	SOMATICSEQ(ABRA_BAM.out.join(MUTECT2.out.join(VARDICT.out.join(DEEPSOMATIC.out.join(LOFREQ.out.join(STRELKA.out.join(FREEBAYES.out.join(PLATYPUS.out))))))))
	COMBINE_VARIANTS(MUTECT2.out.join(VARDICT.out.join(DEEPSOMATIC.out.join(LOFREQ.out.join(STRELKA.out.join(FREEBAYES.out.join(PLATYPUS.out)))))))

	// CNV calling
	CNVKIT(ABRA_BAM.out)
	ANNOT_SV(CNVKIT.out)
	IFCNV(ABRA_BAM.out.collect())

	// IGV reports
	IGV_REPORTS(SOMATICSEQ.out)

	// update_db(SOMATICSEQ.out.collect())
	// cava(SOMATICSEQ.out.join(COMBINE_VARIANTS.out))
	// format_somaticseq_combined(SOMATICSEQ.out)
	// format_concat_combine_somaticseq(format_somaticseq_combined.out)
	// format_pindel(pindel.out.join(COVERAGE.out))
	// format_pindel_UBTF(pindel_UBTF.out.join(COVERAGE.out))
	// merge_csv(format_concat_combine_somaticseq.out.join(cava.out.join(COVERVIEW.out.join(format_pindel.out.join(CNVKIT.out.join(SOMATICSEQ.out.join(FILT3R.out.join(format_pindel_UBTF.out))))))))
	// update_freq(merge_csv.out.collect())
	// Final_Output(COVERAGE.out.join(CNVKIT.out))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	println "Completed at: ${workflow.complete}"
	println "Total time taken: ${workflow.duration}"
}