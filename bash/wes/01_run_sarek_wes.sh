#!/bin/bash

# Usage
# from base dir run: bash bash/wes/01_run_sarek_wes.sh

ref_dir=igenomes_ref
NXF_SINGULARITY_CACHEDIR=sarek-2.7.2/singularity-images/

export NXF_SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR
export NXF_OPTS='-Xms1g -Xmx10g'

target_bed=metadata/Ensembl_e99_Twist_RefSeq_exome_targets_hg38.pad150.merged.sorted.bed

mkdir ./sarek_run
cp bash/wes/input.tsv sarek_run/
cd sarek_run

nextflow info -d > nextflow.info
nextflow run sarek-2.7.2/workflow/ \
	-with-trace trace_complete.txt \
	-with-timeline timeline_complete.html \
	-with-dag dag_complete.pdf \
	--tools 'mutect2,strelka,manta' \
	--input 'input.tsv' \
	--genome GRCh38 \
	--igenomes_base ${ref_dir} \
	--target_bed ${target_bed} \
	--max_time '240.h' \
	--max_memory '12.GB' \
	--max_cpus 4 \
	--cpus 4 \
	--single_cpu_mem '12.GB' \
	-profile singularity > complete.log 2>&1
