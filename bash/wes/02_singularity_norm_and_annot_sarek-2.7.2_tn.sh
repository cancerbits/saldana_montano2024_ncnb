#!/bin/bash

# Usage:
# from base dir run: bash bash/wes/02_singularity_norm_and_annot_sarek-2.7.2_tn.sh sarek_run/results/ <n_cores> > sarek_run/vep_hgvs_tn.log 2>&1

# This script expects the nf-core/sarek release 2.7.2 singularity images.
# You can download the sarek-2.7.2 workflow including the image with the download_sarek-2.7.2.sh script provided in this repository.

# This script expects an igenomes reference base.
# You can download the igenomes reference base with the download_igenomes.sh script provided in this repository.

ref_dir=igenomes_ref

singularity exec \
	--bind $1:$1 \
	--bind $ref_dir:$ref_dir \
	sarek-2.7.2/singularity-images/nfcore-sarek-2.7.2.img \
	bash bash/wes/02a_bcftools_norm_tn.sh $1 $2

singularity exec \
	--bind $1:$1 \
	--bind $ref_dir:$ref_dir \
	sarek-2.7.2/singularity-images/nfcore-sarekvep-2.7.2.GRCh38.img \
	bash bash/wes/02b_vep_hgvs_hg38_sarek-2.7.2_tn.sh $1 $2

