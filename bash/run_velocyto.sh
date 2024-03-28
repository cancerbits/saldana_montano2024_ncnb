#! /bin/bash

# execute this script from within the main project directory, i.e. ./bash/run_velocyto.sh

# run the R script to generate the filtered barcode files
Rscript R/velocyto_helper.R

# activate environment
source /opt/conda/etc/profile.d/conda.sh
conda activate /opt/conda/envs/velocyto

# parse config parameters:
source bash/parse_yaml.sh
eval $(parse_yaml config_bob.yaml CONF_)

IN=metadata/velocyto_libraries.csv

tmpdir=${CONF_tmp_root}/${CONF_project_name}_velocyto
mkdir -p ${tmpdir}
echo "tmpdir: [${tmpdir}]"

while IFS="," read -r path_prefix library_name orig_ident;
do
  echo "--------------------"
  echo "path_prefix: [${path_prefix}]"
  echo "library_name: [${library_name}]"
  echo "orig_ident: [${orig_ident}]"
  
  outdir=${CONF_out_root}/scrna/velocyto
  mkdir -p $outdir
  echo "outdir: [${outdir}]"
  
  indir=${CONF_out_root}/${path_prefix}${library_name}/outs
  inbams=$(find ${indir} -name *.bam)
  bigbam=${tmpdir}/${library_name}_possorted_alignments.bam
  
  echo need to merge these bam files: ${inbams}
  
  if [[ -f ${bigbam} ]]
  then
    echo file ${bigbam} already exists, do not run samtools merge
  else
    echo run samtools merge for library ${library_name}
    samtools merge --threads ${CONF_n_threads_max} -o ${bigbam} ${inbams}
    echo run samtools index for library ${library_name}
    samtools index --threads ${CONF_n_threads_max} ${bigbam}
  fi
  
  sortedbam=${tmpdir}/cellsorted_${library_name}_possorted_alignments.bam
  if [[ -f ${sortedbam} ]]
  then
    echo file ${sortedbam} already exists, do not run samtools sort
  else
    echo run samtools sort for library ${library_name}
    samtools sort -@ ${CONF_n_threads_max} -t CB -O BAM -o $sortedbam $bigbam
  fi
  
  # define velocyto input besides the bam file
  barcodefile=${tmpdir}/${library_name}_filtered_barcodes.tsv
  repeatmaskerfile=${CONF_resource_root}/path/to/GRCh38/GRCh38_RepeatMasker_track_rmsk.gtf
  gtffile=${CONF_resource_root}/${CONF_gtf}

  if [[ -f ${outdir}/${library_name}.loom ]]
  then
    echo file already ${outdir}/${library_name}.loom exists, skipping
  else
    echo starting velocyto on library ${library_name}
    velocyto run --verbose --verbose \
                 --sampleid ${library_name} \
                 --bcfile ${barcodefile} \
                 --outputfolder ${outdir} \
                 --mask ${repeatmaskerfile} \
                 --samtools-threads ${CONF_n_threads_max} \
                 ${bigbam} ${gtffile}
  fi

done < <(tail -n +2 ${IN})
