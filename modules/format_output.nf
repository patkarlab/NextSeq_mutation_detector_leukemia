#!/usr/bin/env nextflow

process CAVA {
	tag "${Sample}"
	input:
		tuple val(Sample), file (somaticVcf), file (somaticseqMultianno), file(cancervarMultianno), file(combinedVcf)
	output:
		tuple val(Sample), file ("*.cava.csv")
	script:
	"""
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${somaticVcf} -o ${Sample}.somaticseq
	${params.cava_path}/cava -c ${params.cava_path}/config_v2.txt -t 10 -i ${combinedVcf} -o ${Sample}.combined
	cava.py ${Sample}.somaticseq.txt ${Sample}.combined.txt ${Sample}.cava.csv
	"""
}

process FORMAT_SOMATICSEQ_COMBINED {
	tag "${Sample}"
	input:
		tuple val (Sample), file(somaticseqVcf), file (multianno), file (cancervarMultianno)
	output:
		tuple val (Sample), file("*.somaticseq.csv")
	script:
	"""
	somaticseqoutput-format.py ${multianno} ${Sample}.somaticseq.csv
	"""
}

process FORMAT_CONCAT_SOMATICSEQ_COMBINED {
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.concat.csv'
	input:
		tuple val (Sample), file (somaticseqCsv)
	output:
		tuple val (Sample), file ("${Sample}.final.concat.csv"), file ("${Sample}.artefacts.csv")
	script:
	"""
	sed -i '1d' ${somaticseqCsv}

	python3 ${params.format_remove_artefact_script} ${somaticseqCsv} ${params.artefactFile} ./${Sample}.final.concat.csv ./${Sample}.artefacts.csv
	sed -i '1iChr,Start,End,Ref,Alt,Variant_Callers,FILTER,SOMATIC_FLAG,VariantCaller_Count,REF_COUNT,ALT_COUNT,VAF,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,Gene_full_name.refGene,Function_description.refGene,Disease_description.refGene,cosmic84,PopFreqMax,1000G_ALL,ExAC_ALL,CG46,ESP6500siv2_ALL,InterVar_automated' ${Sample}.final.concat.csv
	sed -i '1iChr,Start,End,Ref,Alt,Variant_Callers,FILTER,SOMATIC_FLAG,VariantCaller_Count,REF_COUNT,ALT_COUNT,VAF,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,Gene_full_name.refGene,Function_description.refGene,Disease_description.refGene,cosmic84,PopFreqMax,1000G_ALL,ExAC_ALL,CG46,ESP6500siv2_ALL,InterVar_automated' ${Sample}.artefacts.csv
	"""
}

process FORMAT_PINDEL {
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_final.pindel.csv'
	input:
		tuple val (Sample), file (pindelMultianno), file (countsBed), file (pindelCountsBed), file (pindelUBTFCountsBed)
	output:
		tuple val (Sample), file("*_final.pindel.csv")
	script:
	"""
	python3 ${params.format_pindel_script} ${pindelCountsBed} ${pindelMultianno} ${Sample}_final.pindel.csv
	"""
}

process FORMAT_PINDEL_UBTF {
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_final.pindel_ubtf.csv'
	input:
		tuple val (Sample), file (pindelMultianno), file (countsBed), file (pindelCountsBed), file (pindelUBTFCountsBed)
	output:
		tuple val (Sample), file("*_final.pindel_ubtf.csv")
	script:
	"""
	python3 ${params.format_pindel_script} ${pindelUBTFCountsBed} ${pindelMultianno} ${Sample}_final.pindel_ubtf.csv
	"""
}

process MERGE_CSV {
	input:
		tuple val (Sample), file (finalConcat), file (artefacts), file (cavaCsv), file (coverviewRegions) ,file (finalPindel), file(finalCns), file(finalCnr), file(geneScatter), file (finalScatter), file (finalDiagram), file (somaticVcf), file (somaticseqMultianno), file (cancervarMultianno), file(filt3rCsv), file(finalPindelUBTF), file(ext_vcf)
	output:
		val Sample
	script:
	"""
	sed -i 's/\t/,/g' ${finalCnr}
	python3 ${params.pharma_marker_script} ${Sample} ./ ${params.pharma_input_xlxs} ./${Sample}_pharma.csv
	flt3_ext_format.py -v ${ext_vcf} -o ${Sample}_ext.tsv
	python3 ${params.merge_csvs_script} ${Sample} ./ ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ./ ${coverviewRegions} ${finalPindel} ${finalCnr} ./${Sample}_pharma.csv ${finalPindelUBTF} ${Sample}_ext.tsv

	cp ${finalConcat} ${Sample}.final.concat_append.csv
	${params.vep_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.somaticseq.vcf ${PWD}/Final_Output/${Sample}/${Sample}
	${params.vep_extract_path} ${Sample}.final.concat_append.csv ${PWD}/Final_Output/${Sample}/${Sample}_vep_delheaders.txt > ${Sample}.vep
	${params.cancervar_extract} ${cancervarMultianno} ${Sample}.vep ${Sample}_cancervar.csv
	${params.pcgr_cpsr_script_path} ${PWD}/Final_Output/${Sample}/${Sample}.xlsx ${Sample}_cancervar.csv ${filt3rCsv}

	mv output_temp.xlsx ${PWD}/Final_Output/${Sample}/${Sample}.xlsx
	
	"""
}

process FINAL_OUTPUT {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.Low_Coverage.png'
	input:
		tuple val (Sample), file(countsBed), file(pindelCountsBed), file(pindelUBTFCountsBed), file(finalCns), file(finalCnr), file(geneScatter), file (finalScatter), file (finalDiagram)
	output:
		tuple val (Sample), file("*.Low_Coverage.png")
	script:
	"""
	python3 ${params.coveragePlot_script} ${Sample} ${countsBed} ./
	"""
}

process UPDATE_FREQ {
	input:
		val (Sample)
	output:
		path "*.xlsx"
	script:
	"""
	ln -s ${PWD}/work/freq.txt ./

	for i in `cat ${params.input} | grep -i 'myo' | sed 's/-[[:alpha:]]*//g'`
	do
		if [ -f ${PWD}/Final_Output/\${i}/\${i}.xlsx ]; then
			${params.update_freq_excel} ${PWD}/Final_Output/\${i}/\${i}.xlsx freq.txt

			ln -s ${PWD}/Final_Output/\${i}/\${i}.xlsx ./\${i}.xlsx
		fi
	done
	"""
}

process UPDATE_DB {
	publishDir "$PWD/work/", mode: 'copy', pattern: 'freq.txt'
	input:
		val Sample
	output:
		file ("freq.txt")
	script:
	"""
	for i in `cat ${params.input} | grep -i 'myo' | sed 's/-[[:alpha:]]*//g'`
	do 
		if [ -f ${PWD}/Final_Output/\${i}/\${i}.somaticseq.vcf ]; then
			ln -s ${PWD}/Final_Output/\${i}/\${i}.somaticseq.vcf ./
		fi
	done
	files=\$(ls *.somaticseq.vcf)
	${params.updatedb} ${params.alpdb} \${files}
	${params.freq_py} ${params.alpdb} freq.txt
	"""
}