#!/bin/bash

# to download references from igenomes awscli is needed
# conda create --name aws
# conda activate aws
# conda install -c conda-forge awscli

# conda activate aws
# from base dir run: bash bash/wes/00_download_igenomes.sh

aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ ./igenomes_ref
