#! /usr/bin/bash 

samplesheet=$1

source activate new_base
for samples in `cat ${samplesheet}`
do
	echo $samples
	#mkdir -p Final_Output/$samples
	#mkdir -p /home/pipelines/NextSeq_mutation_detector_leukemia/${samples}/cnvkit

	#/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit.sh Final_Output_old/${samples}/${samples}.final.bam /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_ref_GeneNames/Reference_labelled.cnn /home/pipelines/NextSeq_mutation_detector_leukemia/${samples}/cnvkit/

	#/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/gene_scatter/custom_scatter_v3.py /home/pipelines/NextSeq_mutation_detector_leukemia/scripts/gene_scatter/chr_list_all.txt /home/pipelines/NextSeq_mutation_detector_leukemia/${samples}/cnvkit/${samples}.final.cnr /home/pipelines/NextSeq_mutation_detector_leukemia/${samples}/cnvkit/${samples}.final.cns ${samples}

	#cp ${samples}gene_scatter.pdf /home/pipelines/NextSeq_mutation_detector_leukemia/${samples}/cnvkit/
	#cp ${samples}gene_scatter.pdf /home/pipelines/NextSeq_mutation_detector_leukemia/Final_Output/${samples}/

	#cp /home/pipelines/NextSeq_mutation_detector_leukemia/${samples}/cnvkit/${samples}.final-scatter.png /home/pipelines/NextSeq_mutation_detector_leukemia/${samples}/cnvkit/${samples}.final-diagram.pdf /home/pipelines/NextSeq_mutation_detector_leukemia/Final_Output/${samples}/

	rm -r ${samples}
	rm ${samples}gene_scatter.pdf
done
