#!/usr/bin/bash 


input_vcf=$1
pcgr_output_dir=$2
sampleid=$3
cpsr_output_dir=$4
ensembl_ids=$5

source activate SIGVEN
conda activate /home/vishram/pcgr

pcgr --input_vcf $input_vcf --pcgr_dir /home/vishram/ --output_dir $pcgr_output_dir --genome_assembly grch37 --sample_id $sampleid --no_docker --call_conf_tag MVDLKFP --force_overwrite --assay TARGETED --target_size_mb 0.8 --vep_regulatory --vep_gencode_all --estimate_signatures --estimate_msi_status --estimate_tmb --tumor_only

cpsr --input_vcf $input_vcf --pcgr_dir /home/vishram/ --output_dir $cpsr_output_dir --genome_assembly grch37 --sample_id $sampleid --no_docker --force_overwrite --custom_list $ensembl_ids --classify_all

conda deactivate
conda deactivate
