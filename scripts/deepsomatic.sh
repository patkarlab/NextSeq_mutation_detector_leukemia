#!/usr/bin/bash

BIN_VERSION="1.6.1"
BAM_DIR=$1
YOUR_OUTPUT_DIR=$2
YOUR_OUTPUT_VCF=$3
YOUR_BAM=$4
CONTROL_BAM=$5
genome=$6
bedfile=$7

genome_name=$(basename ${genome})
genome_location=$(echo ${genome} | awk 'BEGIN{OFS=FS="/"} {$NF=""; print $0}')
echo "${genome_name} ${genome_location}"

control_bam=$(basename ${CONTROL_BAM})
control_bam_location=$(echo ${CONTROL_BAM} | awk 'BEGIN{OFS=FS="/"} {$NF=""; print $0}')

bedfile_name=$(basename ${bedfile})
bedfile_path=$(echo ${bedfile} | awk 'BEGIN{OFS=FS="/"} {$NF=""; print $0}')

#docker run \
#-v ${INPUT_DIR}:${INPUT_DIR} \
#-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
#google/deepsomatic:"${BIN_VERSION}" \
#run_deepsomatic \
#--model_type=WGS \
#--ref=${INPUT_DIR}/GRCh38_no_alts_chr1.fasta \
#--reads_normal=${INPUT_DIR}/S1395_WGS_ilm_normal.bwa.dedup.chr1.quickstart.bam \
#--reads_tumor=${INPUT_DIR}/S1395_WGS_ilm_tumor.bwa.dedup.chr1.quickstart.bam \
#--output_vcf=${OUTPUT_DIR}/HCC1395_deepsomatic_quickstart.vcf.gz \
#--output_gvcf=${OUTPUT_DIR}/HCC1395_deepsomatic_quickstart.g.vcf.gz \
#--sample_name_tumor="tumor" \
#--sample_name_normal="normal" \
#--num_shards=1 \
#--logging_dir=${OUTPUT_DIR}/logs \
#--intermediate_results_dir ${OUTPUT_DIR}/intermediate_results_dir \
#--regions=chr1:10,000,000-10,100,000

docker run -v ${BAM_DIR}:"/input" -v ${YOUR_OUTPUT_DIR}:"/output" -v ${genome_location}:"/genome_dir" -v ${control_bam_location}:"/control_dir" -v ${bedfile_path}:"/bedfile_dir" google/deepsomatic:"${BIN_VERSION}" run_deepsomatic --model_type=WGS --ref=/genome_dir/${genome_name} --reads_tumor=/input/${YOUR_BAM} --reads_normal=/control_dir/${control_bam} --output_vcf=/output/${YOUR_OUTPUT_VCF} --num_shards=16 --regions=/bedfile_dir/${bedfile_name}

grep '^#' output/${YOUR_OUTPUT_VCF} > temp.vcf
grep -v '^#' output/${YOUR_OUTPUT_VCF} | awk '{ if ($7=="PASS") print $0}' >> temp.vcf
mv temp.vcf ${YOUR_OUTPUT_VCF}
