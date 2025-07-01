# MYOPOOL 

This is a nextflow pipeline for analysing target DNA sequencing data from the custom panel -MYOPOOL 

For running this pipeline, following programs need to be installed and their complete paths need to be added in the params section of the `nextflow.config`.

- adaptors = Fasta file of adapter sequences for trimming
- genome = Genomic fasta file
- samtools = samtools executable path
- bedtools = bedtools executable path
- flt3_bedfile = path to the flt3 bedfile
- ubtf_bedfile =  path to the ubtf bedfile
- site1 = known_polymorphic_sites 1 (Mills_and_1000G_gold_standard.indels)
- site2 = known_polymorphic_sites 2 (dbsnp_138)
- site3 = known_polymorphic_sites 3 (1000G_phase1.snps.high_confidence)
- picard_path = path to the picard.jar file
- GATK38_path = path to the GenomeAnalysisTK-3.8 jar file
- freebayes_path = freebayes executable path 
- platypus_path = path to Platypus.py 
- vardict_path = VarDict executable path
- bcftools_path = bcftools executable path
- strelka_path = path to the strelka bin folder
- lofreq_path = lofreq executable path
- coverview_path = path to the CoverView-1.4.3 folder
- cava_path = path to the CAVA directory
- somaticseq_path = path to the somaticseq_parallel.py
- annovarLatest_path = path to the ANNOVAR folder
- pindel = path to the pindel folder
- get_itd_path = path to the get_itd folder
- java_path = directory containing the java executable
- abra2_path = directory containing the abra jar file
- filt3r_ref = reference fasta for Filt3r
- cnvkit_path = custom CNVkit script path
- vep_script_path = custom VEP script path
- ifcnv = custom ifcnv script path
- deepsomatic = custom deepsomatic script path
- fastp = fastp executable path

## Usage:

1. Keep the `bam`, `bam.bai` and `hard-filtered.vcf` files obtained from Dragen into the `sequences/` folder.

2. Change the `samplesheet.csv`. It should have a list of IDs of the samples. 

3.  Run the following script.

```
./run_nextflow_bamin.sh > script.log
```
This script contains the nextflow command used to execute the workflow.


```
source activate new_base

nextflow -c $PWD/nextflow.config run $PWD/test_nf_scripts/main_bamin.nf -entry MIPS_mocha \
--bedfile <probewise bedfile> \
--bedfile_exonwise <exonwise bedfile> \
--cnvkitRef <cnvkit Reference> \
--gene_scatter_list <chromosome wise list of the exons> \
--gene_scatter <path to the folder containing custom_scatter_chrwise.py> \
-resume -bg 

conda deactivate
```