#!/usr/bin/bash

sample=$1
bamfile=${sample}_tumor.bam

/home/arpit/miniconda3/bin/readCounter --window 10000 --quality 20 \
--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
$bamfile > ${sample}_10kb.tumor.wig

sed -i 's/chr//g' ${sample}_10kb.tumor.wig
sed -i 's/om/chrom/g' ${sample}_10kb.tumor.wig
