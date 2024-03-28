
library(data.table)
library(ggplot2)


resultsDir <- function(...) {
  paste0(config$out_root, "/atac/", ...)
}

pipelineDir <- function(...) {
  paste0(config$out_root, "/pipe/", ...)
}

getBamPath <- function(pipe_name, genome_ver = basename(config$genome_chromatin)) {
  return(pipelineDir("out/", pipe_name, "/aligned_", genome_ver, "/", pipe_name, "_sort_dedup.bam"))
}

getPeaksPath <- function(pipe_name, genome_ver = basename(config$genome_chromatin)) {
  return(pipelineDir("out/", pipe_name, "/peak_calling_", genome_ver, "/", pipe_name, "_peaks.narrowPeak"))
}

dtToGr <- function(dt, chrCol="chrom", startCol="start", endCol="end", metaCols=c()) {
  library("GenomicRanges")
  
  argList <- list()
  for(n in metaCols) {
    argList[[n]] <- dt[,get(n)]
  }
  
  argList$ranges <- IRanges(dt[,as.numeric(get(startCol))],dt[,as.numeric(get(endCol))])
  argList$seqnames <- dt[,get(chrCol)]
  
  do.call(GRanges, args=argList)
}
grToDt <- function(gr, chrCol="chrom", startCol="start", endCol="end") {		
  dt <- data.table(chrom=as.character(seqnames(gr)), start=start(gr), end=end(gr))
  setnames(dt, c(chrCol, startCol, endCol))
  dt
}

liftOver <- function(gr, from="hg19", to=config$genome_build) {
  
  if(from!=to) {
    # download liftOver chain file:
    chnF <- resultsDir("downloaded_files/", from, "To", capFirst(to), ".over.chain")
    if(!file.exists(chnF)) {
      dir.create(resultsDir("downloaded_files"), showWarnings=FALSE)
      download.file(paste0("https://hgdownload-test.gi.ucsc.edu/goldenPath/",from,"/liftOver/", basename(chnF), ".gz"), destfile=paste0(chnF,".gz"))
      R.utils::gunzip(paste0(chnF,".gz"))
    }
    
    chn <- rtracklayer::import.chain(chnF)	
    gr <- rtracklayer::liftOver(gr, chn)
  }
  
  return(gr)
}
simpleMotifName <- function(x) {
  gsub("^(.+) .+$", "\\1",x)
}

# adapted from https://github.com/satijalab/seurat/blob/a1294c4d363780548dbf9cc4a4abb3a6078a6d64/R/utilities.R
# to make it such that it can be run on a matrix without dependecies on other Seurat features
# (also made it a lot faster)
calculateModuleScore <- function(
    assay.data,
    features,
    nbin = 24,
    ctrl = 100,
    name = 'Cluster',
    seed = 1,
    nthread = 8
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  features.old <- features
  
  cluster.length <- length(x = features)
  
  data.avg <- Matrix::rowMeans(x = assay.data)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  
  # merged three for-loops into one parallelized apply function for speed:
  feature.scores.allinone <- BiocParallel::bplapply(1:length(features),  function(i) {
    features.use <- features[[i]]
    # pick out `ctrl` random features at a similar expression level as controls for each feature in the query set:
    ctrl.use <- unique(unlist(lapply(features.use, function(f) {
      names(x = sample(
        x = data.cut[which(x = data.cut == data.cut[f])],
        size = ctrl,
        replace = FALSE
      ))
    })))
    
    ctrl.scores <- Matrix::colMeans(x = assay.data[ctrl.use, , drop = FALSE])
    features.scores <- Matrix::colMeans(x = assay.data[features.use, , drop = FALSE])
    
    return(features.scores - ctrl.scores)
  }, BPPARAM=BiocParallel::MulticoreParam(nthread))
  
  names(feature.scores.allinone) <- names(features)
  return(simplify2array(feature.scores.allinone))
}

options("RCACHE.DIR"=paste0(config$out_root, "/rcache/"))

CACHE_GENE_ANNOT <- "gene_annot"
CACHE_ATAC_PEAKS <- "atac_peaks"
CACHE_ATAC_PEAKS_ANNOTATED <- "atac_peaks_annot"
CACHE_ATAC_DDS <- "atac_dds"
CACHE_ATAC_DDS_RES <- "atac_dds_res"
CACHE_ATAC_META <- "atac_meta"
CACHE_ATAC_REGION_SETS <- "atac_roi"
CACHE_ATAC_COMPARISONS <- "atac_cmps"
CACHE_ATAC_MOTIF_LABELS <- "atac_mtf_lbl"
CACHE_ATAC_ENRICHDB_PREFIX <- "atac_enrich"
CACHE_ATAC_GENE_ASSIGNMENTS <- "atac_gene_assign"
CACHE_ATAC_ENRICHMENT_RESULTS <- "atac_enrich_res"
CACHE_ATAC_REGION_ORDER <- "atac_region_order"
CACHE_SCRNA_MARKERS <- "scrna_markers"
CACHE_SCRNA_DATA <- "scrna_data"
CACHE_SCRNA_DATA_AGGREGATED <- "scrna_data_agg"
CACHE_SCRNA_EXPR <- "scrna_expr"
CACHE_SCRNA_UMAP <- "scrna_umap"
CACHE_SCRNA_MUT_GENES <- "scrna_mut_genes"
CACHE_SCRNA_MUT_SCORE <- "scrna_mut_score"

COLOR_PALETTES$all_modules <- COLOR_PALETTES$all_mods_id

COLOR_PALETTES$day_fac <- COLOR_PALETTES$day <- COLOR_PALETTES$day[paste0("D",c(0, 3, 9, 14, 19))]
COLOR_PALETTES$condition <- COLOR_PALETTES$condition[c("WT", "17q", "17q1q", "17q1qMYCN")]

CONDITION_ORDER <- names(COLOR_PALETTES$condition)
DAY_ORDER <- names(COLOR_PALETTES$day)

dir.create(resultsDir(), showWarnings=FALSE)