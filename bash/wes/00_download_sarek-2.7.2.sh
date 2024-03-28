#!/bin/bash

# set up a nf-core conda env in your terminal
# conda create --name nf-core
# conda activate nf-core
# conda install nf-core

# conda activate nf-core
# from base dir run: bash bash/wes/00_download_sarek-2.7.2.sh > 00_download_sarek-2.7.2.log 2>&1
nf-core download sarek -r 2.7.2 -s -o ./sarek-2.7.2 
