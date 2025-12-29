#!/usr/bin/bash

source deactivate

Sample=$1
offtarget_bam=$2
bedfile=$3
bedfile_exonwise=$4

/home/arpit/miniconda3/bin/readCounter --window 1000000 --quality 20 \
--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
${offtarget_bam} > ${Sample}.tumor.wig

sed -i 's/chr//g' ${Sample}.tumor.wig
sed -i 's/om/chrom/g' ${Sample}.tumor.wig

mkdir ${Sample}_normalized

Rscript /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/ichorCNA_offtarget/normalize_offtarget.R \
  --id ${Sample} \
  --libdir /home/diagnostics/pipelines/ichorCNA-0.4.0/ \
  --offTargetFuncs /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/ichorCNA_offtarget/utils.R \
  --TUMWIG ${Sample}.tumor.wig \
  --NORMWIG /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/ichorCNA_offtarget/BNC2-CNVValMyOPoolV2.tumor.wig \
  --baitBedTum ${bedfile} \
  --gcWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/gc_hg19_1000kb.wig \
  --mapWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/map_hg19_1000kb.wig \
  --centromere /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --genomeBuild hg19 \
  --outDir ${Sample}_normalized


mkdir ${Sample}_offtarget

Rscript /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/ichorCNA_offtarget/runIchorCNA_offTarget.R \
  --libdir /home/diagnostics/pipelines/ichorCNA-0.4.0/ \
  --id ${Sample} \
  --logRFile ${Sample}_normalized/${Sample}_offTarget_cor.txt \
  --statsFile ${Sample}_normalized/${Sample}_readStats.txt \
  --gcWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/gc_hg19_1000kb.wig \
  --mapWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/map_hg19_1000kb.wig \
  --repTimeWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/RepTiming_hg19_1000kb.wig \
  --normalPanel /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
  --ploidy "c(2,3)" \
  --normal "c(0.5,0.6,0.7,0.8,0.9)" \
  --maxCN 5 \
  --includeHOMD False \
  --chrs "c(1:22, 'X')" \
  --chrTrain "c(1:22)" \
  --genomeBuild hg19 \
  --estimateNormal True \
  --estimatePloidy True \
  --estimateScPrevalence True \
  --scStates "c(1,3)" \
  --likModel t \
  --centromere /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --exons.bed ${bedfile_exonwise} \
  --txnE 0.9999 \
  --txnStrength 10000 \
  --minMapScore 0.9 \
  --fracReadsInChrYForMale 0.002 \
  --maxFracGenomeSubclone 0.5 \
  --maxFracCNASubclone 0.5 \
  --normal2IgnoreSC 0.95 \
  --scPenalty 5 \
  --plotFileType pdf \
  --plotYLim "c(-2,4)" \
  --outDir ${Sample}_offtarget
