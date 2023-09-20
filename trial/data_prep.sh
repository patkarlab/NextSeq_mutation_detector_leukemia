#! /usr/bin/bash

sample=$1
outdir=$2

mocha_path="/home/programs/MoChA"
export PATH="${mocha_path}:$PATH"
export BCFTOOLS_PLUGINS="${mocha_path}/bcftools-1.16/plugins/"

#${mocha_path}/Allelic_Depth_calc.py ${sample}.somaticseq.vcf ${sample}_AD.somaticseq.vcf # Adding AD values
#${mocha_path}/bcftools +mochatools ${sample}_AD.somaticseq.vcf --no-version -Ob -o ${sample}_AD_GC.somaticseq.vcf -- -t GC -f /home/reference_genomes/hg19_broad/hg19_all.fasta	# Adding GC values

#${mocha_path}/bcftools +fill-tags ${sample}_AD_GC.somaticseq.vcf -Ob -o persamp_DP.vcf -- -t 'FORMAT/DP:1=int(smpl_sum(FORMAT/AD))'
#${mocha_path}/bcftools +fill-tags persamp_DP.vcf -Ob -o out.vcf -- -t 'DP:1=int(sum(FORMAT/DP))'

#sed -i 's/chr//g' out.vcf

${mocha_path}/binary_vcf.sh ${sample} ${outdir} ${mocha_path}
