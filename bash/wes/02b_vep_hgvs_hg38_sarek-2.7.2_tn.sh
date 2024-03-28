#!/bin/bash -euo pipefail

# The command loops through sarek results/VariantCalling directory and writes to
# results/Annotation/${sample}/VEP and results/reports/${sample}/VEP directories.

# Run with sarek-2.7.2/singularity-images/nfcore-sarekvep-2.7.2.GRCh38.img

input_path=$1
threads=$2

ref_dir=igenomes_ref

if [[ $1 != *results/ ]]; then
	echo "Error: You must provide a sarek-2.7.2 results/ directory path." 1>&2
	exit 1
fi

for vcf in ${1}VariantCalling/*_vs_*/*/*.norm.vcf.gz
do
	# we do not annotate Manta and unfiltered Mutect2 results.
	if [[ ${vcf} == *Manta* ]] || [[ ${vcf} == *unfiltered* ]]; then
		continue
	else
		echo "annotating ${vcf}"

		sample=$(echo ${vcf} | sed -e 's/.*VariantCalling\///g' -e 's/\/.*//g')

		outname=$(basename ${vcf} | sed 's/.norm.vcf.gz/_VEP_hgvs.norm.ann.vcf/g')
		outfile=${1}Annotation/${sample}/VEP/${outname}

		statsname=$(basename ${vcf} | sed 's/.vcf.gz/_VEP_hgvs.summary.html/g')
		statsfile=${1}Reports/${sample}/VEP/${statsname}

		outdir="$(dirname ${outfile})"
		statsdir="$(dirname ${statsfile})"

		if [ -f "${outfile}.gz" ]; then
			echo "annotated file already exists."
			continue
		fi

		mkdir -p ${outdir}
		mkdir -p ${statsdir}

		# https://github.com/nf-core/sarek/blob/e8f56e5bbb6f8a34793a6b5a2945981eb43090aa/main.nf#L3762
		vep \
			-i ${vcf} \
			-o ${outfile} \
			--assembly GRCh38 \
			--species homo_sapiens \
			--fasta ${ref_dir}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
			--offline \
			--cache \
			--cache_version 99 \
			--dir_cache /.vep \
			--everything \
			--filter_common \
			--fork $2 \
			--format vcf \
			--per_gene \
			--stats_file ${statsfile} \
			--total_length \
			--vcf

		echo "bgzip ${outfile}"

		# https://github.com/nf-core/sarek/blob/e8f56e5bbb6f8a34793a6b5a2945981eb43090aa/main.nf#L3804
		bgzip --threads $2 -c ${outfile} > ${outfile}.gz
		tabix ${outfile}.gz
		rm ${outfile}

		echo -e "done annotating ${vcf}\n"
	fi
done

