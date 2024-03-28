#!/bin/bash -euo pipefail

# The command loops through sarek results/VariantCalling directory and writes to
# the respective subdirectories.

# Run with sarek-2.7.2/singularity-images/nfcore-sarek-2.7.2.img 

input_path=$1
threads=$2

ref_dir=igenomes_ref

if [[ $1 != *results/ ]]; then
	echo "Error: You must provide a sarek-2.7.2 results/ directory path." 1>&2
	exit 1
fi
			
for vcf in ${1}VariantCalling/*_vs_*/*/*.vcf.gz
do
	outfile=$(echo ${vcf} | sed 's/.vcf.gz/.norm.vcf.gz/g')

	# we do not normalize twice
	if [[ ${vcf} == *.norm.* ]]; then
		echo "normalized file already exists."
		continue

	# we do not annotate Manta and unfiltered Mutect2 results.
	elif [[ ${vcf} == *Manta* ]] || [[ ${vcf} == *unfiltered* ]]; then
		continue

	# tumor-normal-paired Strelka vcfs are already normalized.
	elif [[ ${vcf} == *Strelka*vs* ]]; then
		cp ${vcf} ${outfile}
		echo -e "tumor-normal-paired Strelka vcf is already normalized.\n"

	else
		echo -e "normalizing ${vcf}\n"

		if [[ ${vcf} == *Mutect2_filtered* ]]; then
			# Changing the "," delimiter to a "|" delimiter in "AS_FilterStatus" 
			# and "AS_SB_TABLE because readVcf() has issues with "," delimiters 
			# within fields.
			# AS_SB_TABLE="REF forward|REF reverse|ALT forward|ALT reverse read counts" 
			mkdir tmp
			cd tmp
			zcat ${vcf} | \
				grep \# > HEAD.txt
			zcat ${vcf} | \
				grep -v \# > BODY.txt
			awk '{print $8}' BODY.txt > INFO.txt
			awk -F ";" '{print $1 ";" $2}' INFO.txt | \
				sed 's/,/|/g' > AS_format.txt
			cut --complement -d ";" -f 1,2 INFO.txt > no_AS_INFO.txt
			paste -d ";" AS_format.txt no_AS_INFO.txt > AS_format_INFO.txt
			cut -f 1-7 BODY.txt > BODY_l.txt
			cut --complement -f 1-8 BODY.txt > BODY_r.txt
			paste -d "\t" BODY_l.txt AS_format_INFO.txt BODY_r.txt > AS_format_BODY.txt
			cat HEAD.txt AS_format_BODY.txt > AS_format_VCF.vcf

			bcftools norm \
				--threads $2 \
				--check-ref e --fasta-ref=${ref_dir}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
				--multiallelics=-both \
				--output=${outfile} \
				--output-type=z AS_format_VCF.vcf

			cd ..
			rm -rf tmp
		else
			bcftools norm \
				--threads $2 \
				--check-ref e --fasta-ref=${ref_dir}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
				--multiallelics=-both \
				--output=${outfile} \
				--output-type=z ${vcf}
		fi
		echo -e "done normalizing ${vcf}\n"
	fi
done

