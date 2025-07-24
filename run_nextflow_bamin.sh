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
##for CNVpanel+ALP:/home/pipelines/mutation_detector_nextflow/bedfile/ALP_CNV_backbone_sorted
##for CNVpanel:/home/pipelines/mutation_detector_nextflow/bedfile/xgen-human-cnv-backbone-hyb-panel-probes
##for Lungpanel:/home/pipelines/mutation_detector_nextflow/bedfile/lung_panel_egfr_kras_tp53_sortd
##for Twistmyeloid:/home/pipelines/mutation_detector_nextflow/bedfile/Leukemia_Bed_file_MYFU_grch37_sorted
##for Twistlymphoid:/home/pipelines/mutation_detector_nextflow/bedfile/Leukemia_Bedfile_ALL_grch37hglft_genome_ucsc
##for combined_panel:/home/pipelines/mutation_detector_nextflow/bedfile/Leukemia_Panel_Myeloid_2023_Feb_hg37_sortd
##for multiple_myeloma:/home/pipelines/mutation_detector_nextflow/bedfile/myeloma_combined_sortd

echo "WARNING : change the bedfile and the cnv reference"
# for cnvkit reference 
# 06112021_Leukemia_Panel_sorted.bed : "/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_ref_GeneNames/Reference_labelled.cnn" 
# Leukemia_Panel_Myeloid_2023_Feb_hg37_sortd.bed : "/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_combpanel/Reference_combpanel.cnn"
# Leukemia_Panel_Myeloid_2023_Feb_hg37_sortd.bed (dragen bam reference): "/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_dragen/Reference_combpanel.cnn"

# For ALP
#source activate new_base
#nextflow -c /home/pipelines/NextSeq_mutation_detector_leukemia/nextflow.config run test_nf_scripts/main_bamin.nf -entry MIPS \
#--bedfile /home/pipelines/mutation_detector_nextflow/bedfile/Leukemia_Panel_Myeloid_2023_Feb_hg37_sortd \
#--cnvkitRef /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_dragen/Reference_combpanel.cnn \
#--sequences /home/pipelines/NextSeq_mutation_detector_leukemia/sequences/ \
#--input /home/pipelines/NextSeq_mutation_detector_leukemia/samplesheet.csv \
#-resume -bg
#conda deactivate 

#source activate new_base
#nextflow -c /home/pipelines/NextSeq_mutation_detector_leukemia/nextflow.config run test_nf_scripts/main_v5_acmg.nf -entry CNVpanel \
#--bedfile /home/pipelines/mutation_detector_nextflow/bedfile/myeloma_combined_tp53_nras_kras_sortd \
#--cnvkitRef /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_mmpanel/Reference_combpanel.cnn \
#--sequences /home/pipelines/NextSeq_mutation_detector_leukemia/sequences/ \
#--input /home/pipelines/NextSeq_mutation_detector_leukemia/samplesheet.csv \
#-resume -bg

# For CNV myeloid bed panel
#source activate new_base
#nextflow -c /home/pipelines/NextSeq_mutation_detector_leukemia/nextflow.config run test_nf_scripts/main_bamin.nf -entry MIPS_mocha \
#--bedfile /home/pipelines/mutation_detector_nextflow/bedfile/CNV_Small_hg19_newmyeloid_ubtf_sortd \
#--bedfile_exonwise  /home/pipelines/mutation_detector_nextflow/bedfile/MYOPOOL_231224_Rebalanced_sortd \
#--cnvkitRef /home/pipelines/MMpanel/scripts/cnvkit_cnvmyeloid_dragen/Reference_combpanel.cnn \
#--gene_scatter_list /home/pipelines/MMpanel/scripts/cnvkit_cnvmyeloid \
#--gene_scatter /home/pipelines/MMpanel/scripts/gene_scatter \
#--sequences /home/pipelines/NextSeq_mutation_detector_leukemia/sequences/ \
#--input /home/pipelines/NextSeq_mutation_detector_leukemia/samplesheet.csv \
#-resume -bg
#conda deactivate

# For MyOPool
#source activate new_base
#nextflow -c /home/pipelines/NextSeq_mutation_detector_leukemia/nextflow.config run test_nf_scripts/main_bamin.nf -entry MIPS_mocha \
#--bedfile /home/pipelines/NextSeq_mutation_detector_leukemia/bedfiles/MYOPOOL_240125_UBTF_sortd \
#--bedfile_exonwise /home/pipelines/mutation_detector_nextflow/bedfile/MYOPOOL_231224_Rebalanced_sortd \
#--cnvkitRef /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_myopool_lt_2x_ver2/Reference_combpanel.cnn \
#--gene_scatter_list /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_MyOPool_exonwise/ \
#--gene_scatter /home/pipelines/MMpanel/scripts/gene_scatter \
#-resume -bg 
#conda deactivate

# For MyOPool
source activate new_base
nextflow -c /home/pipelines/NextSeq_mutation_detector_leukemia/nextflow.config run scripts/alp.nf -entry MyoPool \
--bedfile /home/pipelines/NextSeq_mutation_detector_leukemia/bedfiles/MYOPOOL_240125_UBTF_sortd \
--bedfile_exonwise /home/pipelines/mutation_detector_nextflow/bedfile/MYOPOOL_231224_Rebalanced_sortd \
--cnvkitRef /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_myopool_lt_2x_ver2/Reference_combpanel.cnn \
--gene_scatter_list /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_MyOPool_exonwise/ \
--gene_scatter /home/pipelines/MMpanel/scripts/gene_scatter \
-resume -bg
conda deactivate
