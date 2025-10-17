sample=$1

mkdir ${sample}_ichorCNA_nobed

Rscript /home/diagnostics/pipelines/ichorCNA-0.4.0/scripts/runIchorCNA_v2.R --id $sample \
  --WIG /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/${sample}.tumor.wig --ploidy "c(2)" \
  --normal "c(0.99)" --maxCN 4 \
  --gcWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/gc_hg19_1000kb.wig \
  --mapWig /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/map_hg19_1000kb.wig \
  --centromere /home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --normalPanel /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/PON_MYOPOOL_280525_nobed_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence False \
  --scStates "c()" --txnE 0.9999 --txnStrength 10000 --outDir ./${sample}_ichorCNA_nobed 
