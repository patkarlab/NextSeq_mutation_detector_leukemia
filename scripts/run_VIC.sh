#!/usr/bin/env bash
Sample=$1
vcf=$2

source activate new_base

sed -i 's/chrX/chr23/g' ${avinput}
sed -i 's/chr//g' ${avinput}

java -jar /home/programs/VIC/VIC/target/VIC-1.0.1.jar -b hg19 -i ${vcf} -o ${Sample} -input_type VCF -db /home/programs/VIC/VIC/vicdb -table_annovar /home/programs/annovar_latest/annovar/table_annovar.pl -convert2annovar /home/programs/annovar_latest/annovar/convert2annovar.pl -annotate_variation /home/programs/annovar_latest/annovar/annotate_variation.pl -d /home/programs/annovar_latest/annovar/humandb -otherinfo true
