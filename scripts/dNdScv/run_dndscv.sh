#!/usr/bin/env bash

input_excel=$1
sample_id=$2

source deactivate
# Making the input
/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/dNdScv/dndscv_input.py ${input_excel} ${sample_id}

export R_LIBS="/home/vishram/R/x86_64-pc-linux-gnu-library/3.6/:$R_LIBS"
Rscript /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/dNdScv/dndscv.R -f ${sample_id}_dndscv.tsv -g ${sample_id}_genelist.tsv -o ${sample_id}_geneout.tsv -v ${sample_id}_varout.tsv

/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/dNdScv/dnds_filter.py ${sample_id}_varout.tsv ${sample_id}_geneout.tsv ${input_excel}