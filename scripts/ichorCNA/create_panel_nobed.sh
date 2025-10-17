Rscript createPanelOfNormals_v2.R \
    -f /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/ichorCNA/wig_files_nobed.txt \
    --gcWig=/home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/gc_hg19_1000kb.wig \
    --mapWig=/home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/map_hg19_1000kb.wig \
    --centromere=/home/diagnostics/pipelines/ichorCNA-0.4.0/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
    -o PON_MYOPOOL_280525_nobed
