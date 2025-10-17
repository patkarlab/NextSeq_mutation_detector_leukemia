Rscript createPanelOfNormals_v2.R \
    -f /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/wig_files_10kb.txt \
    --gcWig=/home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/gc_hg19_10kb.wig \
    --mapWig=/home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/map_hg19_10kb.wig \
    --centromere=/home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
    -e /home/pipelines/mutation_detector_nextflow/bedfile/MYOPOOL_231224_Rebalanced_sortd.bed \
    -o PON_MYOPOOL_280525_10kb
