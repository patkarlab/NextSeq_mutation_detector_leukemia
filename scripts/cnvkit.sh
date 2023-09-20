#/usr/bin/bash 
# This script will run cnvkit in a conda environment

final_bam=$1
cnvkit_ref=$2
outdir=$3

source activate cnvkit

cnvkit.py batch ${final_bam} -r ${cnvkit_ref} -m hybrid --drop-low-coverage --output-dir ${outdir} --diagram --scatter

# cnvkit.py scatter -s /home/pipelines/NextSeq_mutation_detector_leukemia/39654-MM/cnvkit/39654-MM.final.cn{r,s} -o myplot.svg
conda deactivate
