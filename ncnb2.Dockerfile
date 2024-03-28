# pull base image

FROM rocker/tidyverse:4.1.3

# who maintains this image
LABEL maintainer Luis Montano "luis.montano@ccri.at"
LABEL version 4.1.3-v1

# how to build this image
#docker build -t cancerbits/dockr:ncnb2 -f ncnb2.Dockerfile .

# The default CRAN mirror is set to a snapshot from 2022-04-21
# Bioconductor is at version 3.14

# fix it
ENV CRAN=https://packagemanager.posit.co/cran/2022-04-21
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2022-04-21'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site

# change some permissions
RUN chmod -R a+rw ${R_HOME}/site-library # so that everyone can dynamically install more libraries within container
RUN chmod -R a+rw ${R_HOME}/library

# add custom options for rstudio sessions
# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf
# make sure authentication is not needed so often
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf
RUN echo "auth-stay-signed-in-days=365" >> /etc/rstudio/rserver.conf



#install OS-level dependencies
RUN DEBIAN_FRONTEND=noninteractive \
  apt clean && \
  apt autoclean && \
  apt-get update && \
  apt-get install --assume-yes \
  htop \
  nano \
  libigraph-dev \
  libcairo2-dev \
  libxt-dev \
  libcurl4-openssl-dev \
  libcurl4 \
  libxml2-dev \
  openssl \
  libssl-dev \
  wget \
  curl \
  bzip2 \
  libbz2-dev \
  libpng-dev \
  libhdf5-dev \
  pigz \
  libudunits2-dev \
  libgdal-dev \
  libgeos-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libz-dev \
  libarmadillo-dev \
  libglpk-dev \
  jags \
  libgsl-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  cmake \
  libxt6 \
  libglpk40 \
  libgsl0-dev \
  imagemagick \
  libmagick++-dev \
  libmagickcore-dev \
  default-jdk \
  liblzma-dev \
  libgeos-dev \
  python3 \
  python3-pip \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# install a couple of R packages that we commonly use
#RUN apt-get -y update && apt-get -y install \
#  libglpk40 \
#  libxt6 \
#  && apt-get clean \
#  && rm -rf /tmp/* /var/tmp/*
RUN R -e "BiocManager::install(c('markdown', 'Seurat', 'SeuratObject', 'sparseMatrixStats', 'edgeR', 'apeglm', 'DESeq2', 'fgsea', 'hypeR', 'patchwork', 'dplyr', 'tidyr'))"
RUN R -e "BiocManager::install(c('Rhtslib', 'biovisBase', 'Rsamtools'))"
RUN R -e "BiocManager::install(c('GenomicAlignments', 'BSgenome', 'VariantAnnotation'))"
RUN R -e "BiocManager::install(c('rtracklayer', 'GenomicFeatures', 'OrganismDbi'))"
RUN R -e "BiocManager::install(c('ensembldb', 'ggbio'))"
RUN R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')"
RUN R -e "BiocManager::install('ggvenn')"
RUN R -e "BiocManager::install('ggrepel')"



# for Quant-seq read quantification:
WORKDIR /home/tools/
RUN wget -O mmquant.zip https://bitbucket.org/mzytnicki/multi-mapping-counter/get/5e2ee253a55e.zip && \
	unzip mmquant.zip && \
	mv mzy* mmquant
	
WORKDIR /home/tools/
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
  tar xf htslib-1.18.tar.bz2 && \
  cd htslib-1.18 && \
  ./configure && \
  make && \
  make install

WORKDIR /home/tools/
RUN git clone https://github.com/single-cell-genetics/cellsnp-lite.git  && \
  cd cellsnp-lite  && \
  autoreconf -iv  && \
  ./configure  && \
  make && \
  make install 
  
 
# basic R stuff:
RUN R -e "BiocManager::install(c('XML', 'RCurl', 'R.utils', 'digest', 'optparse', 'Cairo', 'RColorBrewer', 'data.table', 'pheatmap', 'viridis', 'cowplot', 'matrixStats', 'reshape2', 'stringdist', 'ggplotify', 'fpc', 'magick', 'ComplexHeatmap', 'ggdendro', 'simpleCache', 'ggrastr', 'fastcluster', 'eulerr'))"

#  for NGS data:
RUN R -e "BiocManager::install(c('rhdf5', 'BiocParallel', 'Biostrings', 'Rsubread', 'ShortRead'))"

# for dealing with genome regions:
RUN R -e "BiocManager::install(c('GenomicRanges', 'liftOver', 'AnnotationDbi', 'plyranges'))"

# for differential analysis:
RUN R -e "BiocManager::install(c('limma'))"

# for DNA sequence motifs:
RUN R -e "BiocManager::install(c('motifmatchr', 'BSgenome.Hsapiens.UCSC.hg38'))"

# for heatmaps:
RUN R -e "BiocManager::install('seriation')"
RUN R -e "BiocManager::install('changepoint')"

# for GRNs:
#RUN pip3 install arboreto
RUN R -e "BiocManager::install('WGCNA')"

# install our packages:
RUN R -e "remotes::install_github(repo = 'cancerbits/canceRbits', ref = '554b34c1f9266ec07c940f722b476bd41d44727c')"
RUN R -e "remotes::install_github(repo = 'cancerbits/DElegate', ref = '9f6b99c0d5949331757b4ab6872dd8e9dc7854f2')"

RUN R -e "BiocManager::install(c('DropletUtils', 'scDblFinder', 'infercnv', 'gghighlight', 'harmony', 'broman'))"

# install sceasy (for scVI):
RUN R -e "remotes::install_github('cellgeni/sceasy', ref = '0cfc0e39da4b3ce4abf19d2171daa0e4d2acdd03')"

# install further packages for trajectories:
RUN R -e "BiocManager::install(c('sp', 'glmGamPoi','slingshot', 'tradeSeq', 'sccore'))"



# install conda (mambaforge):
RUN wget -O conda.sh https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-Linux-x86_64.sh && \
  /bin/sh conda.sh -b -p /opt/conda && \
  rm -f miniconda.sh
RUN chmod -R a+rw /opt/conda
ENV PATH=/opt/conda/bin:$PATH

# create a conda environment for scVI:
RUN R -e "reticulate::conda_create(envname='/opt/conda/envs/scvi', python_version='3.11.5', pip=TRUE)" && \
  R -e "reticulate::conda_install(envname='/opt/conda/envs/scvi', packages=c('scvi_tools==0.20.3', 'scanpy==1.9.3'), pip=TRUE)"

# create a conda environment for arboreto/GRNboost2:
RUN R -e "reticulate::conda_create(envname='/opt/conda/envs/arboreto', python_version='3.8.17', pip=TRUE)" && \
  R -e "reticulate::conda_install(envname='/opt/conda/envs/arboreto', packages=c('arboreto==0.1.6'), pip=TRUE)"


WORKDIR /home/tools/
RUN wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
    tar xf bcftools-1.13.tar.bz2 && \
    cd bcftools-1.13 && \
    ./configure && \
    make && \
    make install


WORKDIR /home/tools/
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar xf samtools-1.18.tar.bz2 && \
    cd samtools-1.18 && \
    ./configure --without-curses && \
    make && \
    make install

# create a conda environment for arboreto/GRNboost2:
RUN export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True &&  R -e "reticulate::conda_create(envname='/opt/conda/envs/vireo', python_version='3.8', pip=TRUE)" && \
  R -e "reticulate::conda_install(envname='/opt/conda/envs/vireo', packages=c('git+https://github.com/single-cell-genetics/vireo@701274d6c7a43c4a6f3cb285e30a34d4a400ecff', 'git+https://github.com/single-cell-genetics/MQuad@6c1ae59996c679687f96ecf715bd7afd1dfe049e', 'scikit-learn', 'numpy==1.19.5'), pip=TRUE)"

# install further packages for trajectories:
RUN R -e "BiocManager::install(c('vcfR'))"

# create a conda environment for velocyto:
RUN conda create --prefix /opt/conda/envs/velocyto python=3.11.5 pip numpy scipy cython numba matplotlib scikit-learn h5py click
RUN conda install -p /opt/conda/envs/velocyto -c bioconda samtools
RUN . /opt/conda/etc/profile.d/conda.sh && \
  conda activate /opt/conda/envs/velocyto && \
  pip install pysam && \
  pip install velocyto
  
# install R packages for velocyto analysis
RUN R -e "BiocManager::install('pcaMethods')"
RUN R -e "remotes::install_github('velocyto-team/velocyto.R', ref = '83e6ed9')"
