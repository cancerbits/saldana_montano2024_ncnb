#!/bin/bash

# parse config parameters:
source parse_yaml.sh
eval $(parse_yaml ../config.yaml CONF_)

PARSE_PIPE_CFG="$1"

if [ $# -eq 0 ] ; then
  echo "usage: ./run_parse_pipe.sh PIPE_CONFIG.yaml"
  exit 1
fi

eval $(parse_yaml $PARSE_PIPE_CFG CONF_)


OUT_ROOT=${CONF_out_root_host}/parse/parse_pipe
DATA_ROOT=${CONF_data_root_host}
RESOURCE_ROOT=${CONF_resource_root_host}
TMP_ROOT=${CONF_tmp_root_host}/ncnb_parse/
PARSE_PIPE=${CONF_parse_pipe_dir}
GENOME_INPUT=$(CONF_genome)
GENOME_NAME=$(basename ${CONF_genome})
SAMTOOLS=${CONF_samtools}
MAX_THREADS=${CONF_max_threads}
PARSE_RUN_ID=${CONF_run_id}
SUBLIB_IDS=${CONF_sublib_ids}
SAMPLE_SHEET=${CONF_project_root_host}/metadata/${CONF_sample_sheet}

mkdir -p ${OUT_ROOT}/
mkdir -p ${TMP_ROOT}/



########## Build reference genome ##########

# use the same reference as for 10x and build Parse-specific indices

if [ ! -d ${RESOURCE_ROOT}/parse_data/${GENOME_NAME} ] ; then
	echo "Build Parse genome reference.."
	${PARSE_PIPE}/make_ref.sh \
		${RESOURCE_ROOT}/parse_data/${GENOME_NAME} \
		${RESOURCE_ROOT}/${GENOME_INPUT}/genes/genes.gtf \
		${RESOURCE_ROOT}/${GENOME_INPUT}/fasta/genome.fa
fi


########## Prepare FASTQ files ##########

# convert uBAM to FASTQs and concatenate multiple sequencing lanes into one

echo "Prepare FASTQs for all sub-libraries (destination $TMP_ROOT )..."

for sublib in $SUBLIB_IDS ; do

	echo "- ${sublib}"

	if [ ! -f ${TMP_ROOT}/${sublib}_R1.fq.gz ] ; then
		
		for bam in ${DATA_ROOT}/${PARSE_RUN_ID}/*/*${sublib}*.bam ; do
			base_bam=$(basename $bam)
			echo "-- convert ${base_bam}"
			fq0=${TMP_ROOT}/${sublib}_tmp_${base_bam}_R0.fq.gz
			fq1=${TMP_ROOT}/${sublib}_tmp_${base_bam}_R1.fq.gz
			fq2=${TMP_ROOT}/${sublib}_tmp_${base_bam}_R2.fq.gz		
			$SAMTOOLS fastq -@ $MAX_THREADS -i --i1 /dev/null --index-format "i*"  -1 $fq1 -2 $fq2 $bam
		done
		
		echo "-- concatanete fastqs"
		cat ${TMP_ROOT}/${sublib}_tmp_*_R1.fq.gz > ${TMP_ROOT}/${sublib}_R1.fq.gz
		cat ${TMP_ROOT}/${sublib}_tmp_*_R2.fq.gz > ${TMP_ROOT}/${sublib}_R2.fq.gz
		
		echo "-- remove intermediates"
	fi
	
done



########## Write sample sheet ##########

# convert the sample annotaiton table to the format required by the pipeline (e.g., tabs instead of commata)

echo "Sample sheet: ${TMP_ROOT}/samples.txt"
if [ ! -f ${TMP_ROOT}/samples.txt ] ; then
	echo "- prepare from  ${SAMPLE_SHEET}"
	awk -F, '{print $3 " " $2}' ${SAMPLE_SHEET} | grep -v sample_name > ${TMP_ROOT}/samples.txt
fi



########## Execute pipeline ##########

# run per sub-library

echo "Run Parse pipeline.."

if [ ! -d ${TMP_ROOT}/${GENOME_NAME} ] ; then
	echo "- make copy of genome index in temp dir ${TMP_ROOT}/${GENOME_NAME}"
	cp -r ${RESOURCE_ROOT}/parse_data/${GENOME_NAME} ${TMP_ROOT}/
fi

for sublib in $SUBLIB_IDS ; do
	if [ ! -f ${TMP_ROOT}/${sublib}/split-pipe_v1_0_6p.log ] ; then
		mkdir -p ${TMP_ROOT}/${sublib}/
		cd ${TMP_ROOT}/${sublib}/

		echo "- run pipeline for ${sublib}.."	
		
		echo "	samples: ${TMP_ROOT}/samples.txt"
		echo "	fq1: ${TMP_ROOT}/${sublib}_R1.fq.gz"
		echo "	out: ${TMP_ROOT}/${sublib}/"
		echo "	genome: ${TMP_ROOT}/${GENOME_NAME}"
		
		echo "${PARSE_PIPE}/run_pipe.sh ${TMP_ROOT}/samples.txt ${TMP_ROOT}/${sublib}_R1.fq.gz ${TMP_ROOT}/${sublib}/ ${RESOURCE_ROOT}/parse_data/${GENOME_NAME}"
		
		${PARSE_PIPE}/run_pipe.sh \
			${TMP_ROOT}/samples.txt \
			${TMP_ROOT}/${sublib}_R1.fq.gz \
			${TMP_ROOT}/${sublib}/ \
			${TMP_ROOT}/${GENOME_NAME}			
	fi
done



rm -f ${TMP_ROOT}/sublib.lis
for sublib in $SUBLIB_IDS ; do
	echo "${TMP_ROOT}/${sublib}" >> ${TMP_ROOT}/sublib.lis
	rsync -a -c -v -h --progress ${TMP_ROOT}/${sublib} ${CONF_out_root_host}/results/parse_pipe/
done
			
${PARSE_PIPE}/combine_libs.sh ${TMP_ROOT}/${GENOME_NAME} ${TMP_ROOT}/sublib.lis ${TMP_ROOT} ${TMP_ROOT}/Parse_combined
	