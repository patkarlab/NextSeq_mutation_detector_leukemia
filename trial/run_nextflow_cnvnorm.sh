#!/usr/bin/bash
##BSUB -J smMIPS_pipeline
##BSUB -n 25
##BSUB -q normal
##-m cn2" to submit jobs to cn2
## or " -m cn3"

##########
#for ENTRY : BEDFILES#
##for LEUKEMIA/MIPS: /home/pipelines/mutation_detector_nextflow/bedfile/06112021_Leukemia_Panel_sorted
##for MIPS (IDT-MRD): /home/pipelines/mutation_detector_nextflow/bedfile/04243058_MRD_Panel_V1_final_sorted 
##for CNVpanel+ALP: /home/pipelines/mutation_detector_nextflow/bedfile/ALP_CNV_backbone_sorted
##for only CNV panel: /home/pipelines/mutation_detector_nextflow/bedfile/xgen-human-cnv-backbone-hyb-panel-probes

source activate new_base

nextflow -c /home/pipelines/NextSeq_mutation_detector_leukemia/nextflow.config run old_nf_scripts/main_v3_abra2_pindel_cancervar_cnvnorm.nf -entry CNVpanel --bedfile /home/pipelines/mutation_detector_nextflow/bedfile/xgen-human-cnv-backbone-hyb-panel-probes --sequences /home/pipelines/NextSeq_mutation_detector_leukemia/sequences/ --input /home/pipelines/NextSeq_mutation_detector_leukemia/samplesheet.csv -resume -bg
