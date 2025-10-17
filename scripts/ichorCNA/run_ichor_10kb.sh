sample=$1

mkdir ${sample}_ichorCNA

Rscript /home/diagnostics/pipelines/ichorCNA-0.4.0/scripts/runIchorCNA_v2.R --id $sample \
  --WIG /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/${sample}_10kb.tumor.wig --ploidy "c(2)" \
  --normal "c(0.99)" --maxCN 4 \
  --gcWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/gc_hg19_10kb.wig \
  --mapWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/map_hg19_10kb.wig \
  --centromere /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --normalPanel /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/PON_MYOPOOL_280525_10kb_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence False \
  --exons.bed /home/pipelines/mutation_detector_nextflow/bedfile/MYOPOOL_231224_Rebalanced_sortd.bed \
  --scStates "c()" --txnE 0.9999 --txnStrength 10000 --outDir ./${sample}_ichorCNA
