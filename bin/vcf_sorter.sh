#!/usr/bin/bash 

infile=$1
outfile=$2

grep "^#" $infile > $outfile
grep -v "^#" $infile | sort -k1,1V -k2,2g >> $outfile

#sed -i 's/##INFO=<ID=MVDLK01,Number=7,Type=Integer,Description="Calling decision of the 7 algorithms: MuTect, VarScan2, VarDict, LoFreq, Strelka, SnvCaller_0, SnvCaller_1">/##INFO=<ID=MVDLKFP,Number=7,Type=String,Description="Calling decision of the 7 algorithms: MuTect, VarScan2, VarDict, LoFreq, Strelka, Freebayes, Platypus">/g' 22VNGS183-ALP_new.somaticseq/22VNGS183-ALP.somaticseq.vcf

#sed -i 's/MVDLK01/MVDLKFP/g' 22VNGS183-ALP_new.somaticseq/22VNGS183-ALP.somaticseq.vcf
