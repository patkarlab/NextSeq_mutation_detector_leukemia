#!/usr/bin/bash

sample=$1
bamfile=$2
outdir=$3

/home/arpit/miniconda3/bin/readCounter --window 1000000 --quality 20 \
--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
${bamfile} > tumor.wig

sed -i 's/chr//g' tumor.wig
sed -i 's/om/chrom/g' tumor.wig
mkdir ${outdir}

/usr/bin/Rscript /home/diagnostics/pipelines/ichorCNA-0.4.0/scripts/runIchorCNA.R --id ${sample} \
  --libdir /home/diagnostics/pipelines/ichorCNA-0.4.0/ \
  --WIG tumor.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/gc_hg19_1000kb.wig \
  --mapWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/map_hg19_1000kb.wig \
  --centromere /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --normalPanel /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ${outdir}

