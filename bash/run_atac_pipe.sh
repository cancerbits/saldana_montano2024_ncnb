#!/bin/bash
#
# This script was used to perform the primary data processing pipelines for the ATAC-seq datasets generated in this study.
# We used the PEPATAC and looper pipeline (http://looper.databio.org/), which requires third-party software and auxiliary resources
# for execution (e.g., genome indices). We recommend starting downstream analysis in R/python instead with the ready-processed read 
# counts per peak/gene files provided in GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219153). 
#
# Dependencies: 
# * local install of looper
# * prebuilt docker container cancerbits/pipes
# * `databio/pepatac` and `cancerbits/pipes` in $CODEBASE
# * genome indices in $RESOURCES/genomes/
# * unalignmed BAM files in $DATA

source parse_yaml.sh
eval $(parse_yaml ../config.yaml CONF_)

# run looper pipelines:
export PEPENV="${CONF_pipes_dir}/looper_config.yaml"
nohup looper run --sp atacseq ${CONF_project_root_host}/metadata/atac_pipe.yaml > "nohup_$(hostname)_${USER}_atacseq.log" 2>&1 &

# run additional TSS score calculation (this step was skipped in the initial pipeline run because of a missing TSS annotation file):
REF=$(basename ${CONF_genome_chromatin})
OUT_ROOT=${CONF_out_root_host}/results/pipeline/
for S in $(cut -f1 -d, ${CONF_project_root_host}/metadata/samples_atac.csv) ; do
	INFILE=$OUT_ROOT/${S}/aligned_${REF}/${S}_sort_dedup.bam
	if [ -f $INFILE ] ; then
		echo $S
		OUTFILE=$OUT_ROOT/${S}/QC_${REF}/${S}_TssEnrichment.txt
		if [ ! -f $OUTFILE ] ; then
			$CODEBASE/pepatac/tools/pyTssEnrichment.py -a $INFILE -b ${RESOURCES}/genomes/${REF}/${REF}_TSS.tsv -p ends -c 8 -e 2000 -u -v -s 4 -o $OUTFILE
		fi
	fi
done
