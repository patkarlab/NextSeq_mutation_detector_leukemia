#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { generatefinalbam; hsmetrics_run; hsmetrics_collect } from './main_bamin.nf'
include { platypus_run; coverage; freebayes_run; mutect2_run; vardict_run; varscan_run; DeepSomatic; lofreq_run; strelka_run; somaticSeqDragen_run } from './main_bamin.nf'

workflow CLL_TP53 {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }
		// samples_ch.view()
	main:
	// filt3r(samples_ch)
	generatefinalbam(samples_ch)
	// getitd(generatefinalbam.out)
	hsmetrics_run(generatefinalbam.out)
	// 
	all_hsmetrics_gw = hsmetrics_run.out.genewise.collect()
	all_hsmetrics_pw = hsmetrics_run.out.exonwise.collect()
	// 
	hsmetrics_collect(all_hsmetrics_gw, all_hsmetrics_pw)
	
	platypus_run(generatefinalbam.out)
	coverage(generatefinalbam.out)
	freebayes_run(generatefinalbam.out)
	mutect2_run(generatefinalbam.out)
	vardict_run(generatefinalbam.out)
	varscan_run(generatefinalbam.out)
	DeepSomatic(generatefinalbam.out)
	lofreq_run(generatefinalbam.out)
	strelka_run(generatefinalbam.out)
	somaticSeqDragen_run(generatefinalbam.out.join(mutect2_run.out.join(vardict_run.out.join(DeepSomatic.out.join(lofreq_run.out.join(strelka_run.out.join(freebayes_run.out.join(platypus_run.out))))))))
	combine_variants(mutect2_run.out.join(vardict_run.out.join(DeepSomatic.out.join(lofreq_run.out.join(strelka_run.out.join(freebayes_run.out.join(platypus_run.out)))))))
	pindel(generatefinalbam.out)
	igv_reports(somaticSeqDragen_run.out)
	update_db(somaticSeqDragen_run.out.collect())
	coverview_run(generatefinalbam.out)
	cava(somaticSeqDragen_run.out.join(combine_variants.out))
	format_somaticseq_combined(somaticSeqDragen_run.out)
	format_concat_combine_somaticseq(format_somaticseq_combined.out)
	format_pindel(pindel.out.join(coverage.out))
	// format_pindel_UBTF(pindel_UBTF.out.join(coverage.out))
	merge_csv(format_concat_combine_somaticseq.out.join(cava.out.join(coverview_run.out.join(format_pindel.out.join(cnvkit_run.out.join(somaticSeqDragen_run.out.join(filt3r.out.join(format_pindel_UBTF.out))))))))
	update_freq(merge_csv.out.collect())
	// Final_Output(coverage.out.join(cnvkit_run.out))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}