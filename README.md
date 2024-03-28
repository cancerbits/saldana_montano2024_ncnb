# Supplementary code repository for: A human neural crest model reveals the developmental impact of neuroblastoma-associated chromosomal aberrations

Authors: Ingrid M. Saldana-Guerrero<sup>1,2,3,\*</sup>, Luis F. Montano-Gutierrez<sup>4,\*</sup>, Katy Boswell<sup>1,2,$</sup>, Christoph Hafemeister<sup>4,$</sup>, Evon Poon<sup>5,$</sup>, Lisa E. Shaw<sup>6,$</sup>, Dylan Stavish<sup>1,2,$</sup>, Rebecca A. Lea<sup>1</sup>, Sara Wernig-Zorc<sup>4</sup>, Eva Bozsaky<sup>4</sup>, Irfete S. Fetahu<sup>4</sup>, Peter Zoescher<sup>4</sup>, Ulrike Pötschger<sup>4</sup>, Marie Bernkopf<sup>4,7</sup>, Andrea Wenninger-Weinzierl<sup>4</sup>, Caterina Sturtzel<sup>4</sup>, Celine Souilhol<sup>1,2,8</sup>, Sophia Tarelli<sup>1,2</sup>, Mohamed R. Shoeb<sup>4</sup>, Polyxeni Bozatzi<sup>4</sup>, Magdalena Rados<sup>4</sup>, Maria Guarini<sup>9</sup>, Michelle C. Buri<sup>4</sup>, Wolfgang Weninger<sup>6</sup>, Eva M. Putz<sup>4</sup>, Miller Huang<sup>10,11</sup>, Ruth Ladenstein<sup>4</sup>, Peter W. Andrews<sup>1</sup>, Ivana Barbaric<sup>1,2</sup>, George D. Cresswell<sup>4</sup>, Helen E. Bryant<sup>3</sup>, Martin Distel<sup>4</sup>, Louis Chesler<sup>5,+</sup>, Sabine Taschner-Mandl<sup>4,+</sup>, Matthias Farlik<sup>6,+</sup>, Anestis Tsakiridis<sup>1,2,#</sup>, Florian Halbritter<sup>4,#</sup>

Affiliations:

1 Centre for Stem Cell Biology, School of Biosciences, The University of Sheffield, Sheffield, UK
2 Neuroscience Institute, The University of Sheffield, Sheffield, UK	
3 Sheffield Institute for Nucleic Acids (SInFoNiA), School of Medicine and Population Health, The University of Sheffield, Sheffield, UK
4 St. Anna Children’s Cancer Research Institute (CCRI), Vienna, Austria
5  Division of Clinical Studies, The Institute of Cancer Research (ICR) & Royal Marsden NHS Trust, London, UK
6 Medical University of Vienna, Department of Dermatology, Vienna, Austria
7 Labdia Labordiagnostik GmbH, Vienna, Austria
8 Biomolecular Sciences Research Centre, Department of Biosciences and Chemistry, Sheffield Hallam University, Sheffield, UK
9 CeMM Research Center for Molecular Medicine of the Austrian Academy of Science, Vienna, Austria 
10 Children’s Hospital Los Angeles, Cancer and Blood Disease Institutes, and The Saban Research Institute, Los Angeles, CA, USA
11 Keck School of Medicine, University of Southern California, Los Angeles, CA, USA

<sup>*,$,+</sup>, These authors contributed equally to this work.
<sup>#</sup> Equally contributing senior authors.

## Abstract:

Early childhood tumours arise from transformed embryonic cells, which often carry large copy number alterations (CNA). However, it remains unclear how CNAs contribute to embryonic tumourigenesis due to a lack of suitable models. Here we employ female human embryonic stem cell (hESC) differentiation and single-cell transcriptome and epigenome analysis to assess the effects of chromosome 17q/1q gains, which are prevalent in the embryonal tumour neuroblastoma (NB). We show that CNAs impair the specification of trunk neural crest (NC) cells and their sympathoadrenal derivatives, the putative cells-of-origin of NB. This effect is exacerbated upon overexpression of MYCN, whose amplification co-occurs with CNAs in NB. Moreover, CNAs potentiate the pro-tumourigenic effects of MYCN and mutant NC cells resemble NB cells in tumours. These changes correlate with a stepwise aberration of developmental transcription factor networks. Together, our results sketch a mechanistic framework for the CNA-driven initiation of embryonal tumours. 

## Repository structure:

* `ncnb2.Dockerfile` defines the environment used to carry out all single-cell and ATAC-seq analyses. 
* `config.yaml` is used to set paths 
* `R/` holds R function definitions and misc utility scripts
* `Rmd/` holds R markdown documents for the individual steps of the project (THIS IS THE MAIN SOURCE OF SOURCE CODE FOR THIS PROJECT)
* `bash/` holds shell scripts to build and run the docker image, and to parse the config file
* `metadata/` holds custom geneset definitions

## Data exploration:

We aimed to provide multiple means for researchers to easily explore and reuse the data generated in this study. All raw and processed data are available from GEO (see Links section below) -- this includes counts and similar files as produced by the initial processing steps (CellRanger etc.) that can be readily downloaded and loaded into R or python.

Additionally, you can access our data using interactive visualization and data exploration tools like:

### R2

You can directly browse the scRNA-seq data generated in this study via the [R2](https://r2.amc.nl/) platform at this link (you can choose between the WT-only dataset as seen in Fig. 2 and the complete, integrated dataset as seen in Fig. 4):
 
[http://r2platform.com/halbritter24/](http://r2platform.com/halbritter24/)

### cellxgene

If you already have the [cellxgene](https://cellxgene.cziscience.com/) software set up on your computer, you can launch an interface to explore the scRNA-seq data in the paper using the following commands (using data hosted on GEO):

WT-only dataset (as seen in Fig. 2; switch to the view "wtumap.fullscvi822" to see the same UMAP layout as in the figure):

```{bash}
cellxgene launch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221853/suppl/GSE221853%5Fncnb%5Fseurat%5Fwt.h5ad"
```

Complete integrated dataset (as seen in Fig. 4; switch to the view "umap.fullscvi8" to see the same UMAP layout as in the figure; note, this dataset is quite large and might take a while to load):

```{bash}
cellxgene launch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221853/suppl/GSE221853%5Fncnb%5Fseurat%5Ffull.h5ad"
```

## Analysis workflow:

The majority of the analysis for this project can be carried out from within RStudio using the R-markdown documents provided in the `Rmd` folder. All software dependencies have been installed in `ncnb2.Dockerfile`, so prospective users will only need to build the Docker (or Podman) image and execute it to start a live Rstudio session. Please configure local paths in the `config.yaml` file. 

Specific instructions / entry points for different parts of the analysis below.

### scRNA-seq Part I. sequential analysis (execution recommended in order):

This code analyses data split into three datasets: 10X G1-G13(without G2), 10X G14-G27 and  One Parse dataset. 

The scripts bash/run_cellranger.sh and run_cellranger_newdata.sh are used to process the 10X Genomics datasets and produce the count matrices. 

The single-cell analysis code is inside the Rmd/scrna folder, with files numbered roughly according to the order of execution. 

*`Scripts 00...` deal with preparing the config files and commands for the execution of cell ranger (ready outputs in metadata and bash folders) 

*`Scripts 01...` deal with import, quality control, demultiplexing and mapping to embryonic adrenal gland references

*`Scripts 02...` deal with assembly of all the sub-datasets into one Seurat object and applying scvi integration to them.

*`Scripts 03 and 05` relate to the global analysis of the wild type dataset and full dataset, and figures requiring only each dataset. These scripts are the basis and produce the main clusters and umaps usd throughout the figures. 

*`Scripts 052, 04 and 06` involve analysis comparing mutant cells, patients' cells and other genotypes to cells in the wild type reference. 

### scRNA-seq Part II. Subsequent analyses that that interconnect the WT, mutant and patient datasets in one way or another, and thus depend on the above scripts. They need not be executed in a precise order. 

*`Script 07` performs an analysis of the chromosomal breakpoints in the study's cell lines as compared to an Austrian cohort of Neuroblastoma patients, and finds relevant genes within these break points. This sectionsection depends on output from the CNA analysis. 

*`Script 08` performs Pseudotime trajectory analysis of the WT dataset using Slingshot.

*`Script 09` reanalyses bulk RNA-seq data from patient datasets and searches for signatures found in the in vitro differentiation dataset. 

*`Script 10` performs RNA velocity analysis with Velocyto. 

*`Script 11` performs mapping, to the WT differentiation reference, of the extended split pool (Parse Biosciences) dataset containing multiple other genotypes.

### ATAC-seq:

The script `bash/run_atac_pipe.sh` executes the primary data analysis pipeline (e.g., alignment and peak calling) using [PEPATAC](https://pepatac.databio.org/) and config parameters from `config.yaml` and `metadata/samples_atac.csv`. This requires some setup of third-party software. Alternatively, you can obtain peak coordinates and read counts from GEO and skip these initial processing steps.

To reproduce the remaining analysis and figures in the paper, use our provided Docker/Podman containers (use `./run_podman_rstudio.sh`) to start up an Rstudio session and follow through the notebooks in `Rmd/atac/` in numeric order. N.B., these scripts also rely on some outputs from the scRNA-seq analysis, so make sure to run this first.

###  Copy number Alteration (CNA) and phylogenetic analysis

NB: Analysis was performed in a separate docker (R version 4.2, Bioconductor 3.16). 
This can be pulled from `docker pull ccribioinf/dockrstudio:4.2.0-v1`.
Session information for each analysis can be found in a `rds` file in each directory (`R/cna_analysis` and `R/phylogeny`).

#### Estimation of copy number changes in whole exome sequencing bulk analysis by Sequenza

To run the analysis you must first use Sequenza on the command line to bin reads in the BAM files.

Step 1: Create the conda environment using `sequenza_conda.yaml`.

Step 2: Run `bash/01_CNA_run_sequenza_seqz_binning.sh` located in the `bash` folder.

Step 3: Within R/Rstudio run `R/cna_analysis/02_CNA_sequenza_analysis.R` to do copy number analysis and create the cohort heatmap.

Step 4: Run `R/cna_analysis/03_CNA_breakpoint_output.R` to get a summary of the aberrations.

### Phylogenetic analysis on single nucleotide variants

Step 1: Run `R/phylogeny/01_PHY_phylogenetic_analysis.R` script on the variant tsv.

### Whole Exome Sequencing Single Nucleotide Variant analysis:

This code performs the analysis of Single Nucleotide variants in the human stem cell lines used in this study, and reports neuroblastoma-relevant variants. The code depends on output generated by the CNA analysis. 

To reproduce this analysis, follow the following steps:
 
Prepare conda environments:
1. Follow the comments in ```bash/wes/00_download_igenomes.sh```
2. Follow the comments in ```bash/wes/00_download_sarek-2.7.2.sh```
 
From the repository base directory, execute
1. ```bash bash/wes/00_download_igenomes.sh```  
This will create an igenomes_ref directory in the base directory
2. ```bash bash/wes/00_download_sarek-2.7.2.sh```  
This will create a sarek-2.7.2 directory in the base directory
3. ```bash bash/wes/01_run_sarek_wes.sh```  
This will create a sarek_run directory in the base directory
4. ```bash bash/wes/02_singularity_norm_and_annot_sarek-2.7.2_tn.sh```  
This will create an Annotation directory in the sarek_run directory
5. ```bash bash/wes/03_collect_hsmetrics.sh```  
This will create an hsmetrics directory in the base directory
 
After the above steps are complete, the R code can be executed in a ccribioinf/dockrstudio:4.0.3-v1 container:
* ```Rmd/wes/ncnb_wes-tn_multisample_sarek-2.7.2_v3.Rmd``` creates per-sample variant reports and "ncnb_smlVars_all_tn_pass_tiers.tsv" (needed for ```R/phylogeny/01_PHY_phylogenetic_analysis.R```)
* ```Rmd/wes/ncnb_coverage_plots.Rmd``` creates coverage barplots and "mean_coverages_all_samples.tsv"
* ```Rmd/wes/ncnb_cnv_filtering.Rmd``` expects "NB_cell_line_all_aberrations_coordinates_wgenes.csv" (from ```Rmd/scrna/07_ChrBreakpoints.Rmd```) and creates "NB_cell_line_all_aberrations_coordinates_whitelistgenes_nbaf100.csv"

 
Prepare conda environments:
1. Follow the comments in ```bash/wes/00_download_igenomes.sh```
2. Follow the comments in ```bash/wes/00_download_sarek-2.7.2.sh```
 
From the repository base directory, execute
1. ```bash bash/wes/00_download_igenomes.sh```  
This will create an igenomes_ref directory in the base directory
2. ```bash bash/wes/00_download_sarek-2.7.2.sh```  
This will create a sarek-2.7.2 directory in the base directory
3. ```bash bash/wes/01_run_sarek_wes.sh```  
This will create a sarek_run directory in the base directory
4. ```bash bash/wes/02_singularity_norm_and_annot_sarek-2.7.2_tn.sh```  
This will create an Annotation directory in the sarek_run directory
5. ```bash bash/wes/03_collect_hsmetrics.sh```  
This will create an hsmetrics directory in the base directory
 
After the above steps are complete, the R code can be executed in a ccribioinf/dockrstudio:4.0.3-v1 container:
* ```Rmd/wes/ncnb_wes-tn_multisample_sarek-2.7.2_v3.Rmd``` creates per-sample variant reports and "ncnb_smlVars_all_tn_pass_tiers.tsv" (needed for ```R/phylogeny/01_PHY_phylogenetic_analysis.R```)
* ```Rmd/wes/ncnb_coverage_plots.Rmd``` creates coverage barplots and "mean_coverages_all_samples.tsv"
* ```Rmd/wes/ncnb_cnv_filtering.Rmd``` expects "NB_cell_line_all_aberrations_coordinates_wgenes.csv" (from ```Rmd/scrna/07_ChrBreakpoints.Rmd```) and creates "NB_cell_line_all_aberrations_coordinates_whitelistgenes_nbaf100.csv"


## Links:

* Pre-print version of the paper: [doi:10.1101/2022.11.21.515753](https://doi.org/10.1101/2022.11.21.515753)
* Gene Expression Omnibus (GEO) entry: [GSE219153](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219153)
* Data scope in the R2 database: [http://r2platform.com/halbritter24/](http://r2platform.com/halbritter24/)
