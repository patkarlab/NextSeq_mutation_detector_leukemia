#!/usr/bin/bash
##BSUB -J smMIPS_pipeline
##BSUB -n 25
##BSUB -q normal
##-m cn2" to submit jobs to cn2
## or " -m cn3"

##########
#for ENTRY : BEDFILES#
##for LEUKEMIA/MIPS: 06112021_Leukemia_Panel_sorted
##for MIPS (IDT-MRD): /home/pipelines/mutation_detector_nextflow/bedfile/04243058_MRD_Panel_V1_final_sorted 
##for CNVpanel+ALP:ALP_CNV_backbone_sorted.bed

source activate new_base
nextflow -c /home/pipelines/NextSeq_mutation_detector_leukemia/nextflow.config run test_nf_scripts/main_v3_abra2_pindel_cancervar.nf -entry MIPS --bedfile /home/pipelines/NextSeq_mutation_detector_leukemia/Varanasi_bedfile/Probes_merged_ok_Premas_HomiBhaba_Ver2_TE-92262737_hg19_220428032712 --sequences /home/pipelines/NextSeq_mutation_detector_leukemia/sequences/ --input /home/pipelines/NextSeq_mutation_detector_leukemia/samplesheet.csv -resume -bg

#conda deactivate 
