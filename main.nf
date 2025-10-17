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
include { FASTQTOBAM; FASTQTOBAM_WGS; ABRA_BAM } from './modules/processes.nf'

// FLT3 ITD detection
include { FILT3R; GETITD } from './modules/flt3_itd.nf'

// HSmetrics calculation
include { HSMETRICS; HSMETRICS_COLLECT } from './modules/hsmetrics.nf'

// COVERAGE calculation
include { COVERAGE; COVERVIEW; COVERAGE_WGS; COVERAGE_WGS_COLLECT } from './modules/coverage.nf'

// Variant calling
include { PLATYPUS; FREEBAYES; MUTECT2; VARDICT; DEEPSOMATIC; LOFREQ; STRELKA; PINDEL; PINDEL_UBTF} from './modules/variant_call.nf'

// Variant integration 
include { SOMATICSEQ; COMBINE_VARIANTS } from './modules/somaticseq.nf'

// CNV calling
include { CNVKIT; ANNOT_SV; IFCNV } from './modules/cnv_call.nf'

// IGV reports
include { IGV_REPORTS } from './modules/igv_reports.nf'

// ichorCNA
include { ICHOR_CNA } from './modules/ichorCNA.nf'

// Format output
include {CAVA; FORMAT_SOMATICSEQ_COMBINED; FORMAT_CONCAT_SOMATICSEQ_COMBINED; FORMAT_PINDEL; FORMAT_PINDEL_UBTF; MERGE_CSV; FINAL_OUTPUT; UPDATE_FREQ; UPDATE_DB} from './modules/format_output.nf'

// DND SCV
include { DNDSCV } from './modules/dnd_scv.nf'

workflow MyoPool {
	leukemia = Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample_full = row[0].trim()
			def sample_base = sample_full.tokenize('-')[0]

			def r1 = file("${params.sequences}/${sample_full}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample_full}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample_full}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample_full}*_R2.fastq.gz", checkIfExists: false)
			}

			tuple(sample_full, sample_base, r1, r2)
		}
		.branch {
			myopool: it[0].contains("MyOPool")
			wgs:     it[0].contains("WGS")
		}

	myopool_ch = leukemia.myopool.map { full, base, r1, r2 -> tuple(base, r1, r2) }
	wgs_ch     = leukemia.wgs.map     { full, base, r1, r2 -> tuple(base, r1, r2) }

	myopool_ch.view { "MYOPOOL: $it" }
	wgs_ch.view     { "WGS: $it" }

	main:
	// Adapter Trimming, alignment and GATK BQSR - MYOPOOL
	myo_bam_ch = FASTQTOBAM(myopool_ch)
	ABRA_BAM (myo_bam_ch)

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

	// Format Output
	CAVA(SOMATICSEQ.out.join(COMBINE_VARIANTS.out))	
	FORMAT_SOMATICSEQ_COMBINED(SOMATICSEQ.out)
	FORMAT_CONCAT_SOMATICSEQ_COMBINED(FORMAT_SOMATICSEQ_COMBINED.out)
	FORMAT_PINDEL(PINDEL.out.join(COVERAGE.out))
	FORMAT_PINDEL_UBTF(PINDEL_UBTF.out.join(COVERAGE.out))
	MERGE_CSV(FORMAT_CONCAT_SOMATICSEQ_COMBINED.out.join(CAVA.out.join(COVERVIEW.out.join(FORMAT_PINDEL.out.join(CNVKIT.out.join(SOMATICSEQ.out.join(FILT3R.out.join(FORMAT_PINDEL_UBTF.out))))))))	
	FINAL_OUTPUT(COVERAGE.out.join(CNVKIT.out))
	UPDATE_FREQ(MERGE_CSV.out.collect())
	UPDATE_DB(SOMATICSEQ.out.collect())
	DNDSCV(UPDATE_FREQ.out.collect())



	// Adapter Trimming, alignment and GATK BQSR - WGS
	wgs_bam_ch = FASTQTOBAM_WGS(wgs_ch)

	// Coverage calculation for WGS
	COVERAGE_WGS(wgs_bam_ch)
	all_coverages_wgs = COVERAGE_WGS.out.collect()
	COVERAGE_WGS_COLLECT(all_coverages_wgs)

	// ichorCNA
	ICHOR_CNA(wgs_bam_ch)

}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	log.info ( "Completed at: ${workflow.complete}")
	log.info ( "Total time taken: ${workflow.duration}")
}
