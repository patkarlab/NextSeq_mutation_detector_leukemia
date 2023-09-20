#! /usr/bin/env bash

annovar_input_file=$1
sample=$2
cancervar_infile=${sample}_cancervar_input.dat
cancervar_outprefix=${sample}myanno

source activate new_base
/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cancervar_input.py ${annovar_input_file} ${cancervar_infile}

cp /home/programs/CancerVar/config.ini config.ini
sed -i -r "s/inputfile = .*/inputfile = ${cancervar_infile}/g" config.ini
sed -i -r "s/outfile = .*/outfile = ${cancervar_outprefix}/g" config.ini
python3 /home/programs/CancerVar/CancerVar.py -c config.ini

python3 /home/programs/CancerVar/OPAI/scripts/feature_preprocess.py -a ${cancervar_outprefix}.hg19_multianno.txt.grl_p -c ${cancervar_outprefix}.hg19_multianno.txt.cancervar -m ensemble -n 5 -d /home/programs/CancerVar/OPAI/saves/nonmissing_db.npy -o ${cancervar_outprefix}.hg19_multianno.txt.cancervar.ensemble.csv

if [ -s ${cancervar_outprefix}.hg19_multianno.txt.cancervar.ensemble.csv ]; then
	# File is not empty
	python3 /home/programs/CancerVar/OPAI/scripts/opai_predictor.py -i ${cancervar_outprefix}.hg19_multianno.txt.cancervar.ensemble.csv -m ensemble -c /home/programs/CancerVar/OPAI/saves/ensemble.pt -d cpu -v ${cancervar_outprefix}.hg19_multianno.txt.cancervar -o ${cancervar_outprefix}.hg19_multianno.txt.cancervar.ensemble.pred

else
	# File is empty
	echo "#Chr    Start   End     Ref     Alt     Ref.Gene        Func.refGene    ExonicFunc.refGene      Gene.ensGene    avsnp147        AAChange.ensGene        AAChange.refGene        clinvar: Clinvar   CancerVar: CancerVar and Evidence      Freq_ExAC_ALL   Freq_esp6500siv2_all    Freq_1000g2015aug_all   Freq_gnomAD_genome_ALL  CADD_raw        CADD_phred      SIFT_score        GERP++_RS       phastCons20way_mammalian        dbscSNV_ADA_SCORE       dbscSNV_RF_SCORE        Interpro_domain AAChange.knownGene      MetaSVM_score   Freq_gnomAD_genome_POPs   OMIM    Phenotype_MIM   OrphaNumber     Orpha   Pathway Therap_list     Diag_list       Prog_list       Polyphen2_HDIV_score    FATHMM_score    MetaLR_score    MutationAssessor_score    cosmic91        icgc28  Otherinfo       ensemble_score" > ${cancervar_outprefix}.hg19_multianno.txt.cancervar.ensemble.pred
fi
#python3 /home/programs/CancerVar/OPAI/scripts/opai_predictor.py -i ${cancervar_outprefix}.hg19_multianno.txt.cancervar.ensemble.csv -m ensemble -c /home/programs/CancerVar/OPAI/saves/ensemble.pt -d cpu -v ${cancervar_outprefix}.hg19_multianno.txt.cancervar -o ${cancervar_outprefix}.hg19_multianno.txt.cancervar.ensemble.pred
