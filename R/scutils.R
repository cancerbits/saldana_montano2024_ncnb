#!/usr/bin/env Rscript

###############################################
#loading libraries
###############################################

cat("loading libraries...\n")

require(scales)
library(devtools)
library(Seurat)
library(dplyr)
library(tidyr)
library(data.table)
library(simpleCache)
library(ggplot2)
library(patchwork)
library(gghighlight)
library(broman)
library(canceRbits)
options(scipen=999)
options(max.print=500)
options(RCACHE.DIR = params$cachedir)



cat("adding important variables and constants...\n")
##########################################
#important variables and constants
##########################################
recreateAll=F
cfg <- list(
        PCA_DIMS = 1:30,
        MIN_FEATURES = 500,
        MIN_CELLS = 20,
        PERC_MITO = 10,
        MAX_DOUBLET_SCORE = 3
)


stagelevels=c("D0", "D3", "D9", "D14", "D19", "Negative")
conditionlevels=c("WT", "17q", "17qx","17qMYCN", "Negative")

###########################################################################################################################
#  PLOTTING PARAMETERS
###########################################################################################################################
cat("adding plotting parameters...\n")
#size unit of png in pixels
w=250
#size unit of pdg in inches
pw=5



##########################################################################################################################
#general functions
##########################################################################################################################
###########################################################################################################################
cat("adding general functions...\n")
#permutate an array
shuffle=function(x) sample(x, length(x), replace=F)

#prepare gene DE table by removing NAs, filtering by cutoff p-value
filterorderDE<-function(x, thresh=0.05, lfthresh=1, netchange=FALSE, top=NA){
  
  if(netchange==TRUE){
  x=x[is.na(x[,pvaluecol])==FALSE & x[,pvaluecol ]<thresh & abs(x[,lfccol])>lfthresh,]
  ord=order(abs(x[,lfccol]), decreasing=T)
  }else{
  x=x[is.na(x[,pvaluecol])==FALSE & x[,pvaluecol ]<thresh & x[,lfccol]>lfthresh,]
  ord=order(x[, lfccol], decreasing=TRUE)
  }
  if(!is.na(top)){
x[ord[1:top], ]
}else{
x[ord, ]
}
}


reddf=function(sc, reduc="pca"){
cbind(metadata(sc), Embeddings(sc, reduction=reduc))
}

harmonydf=function(sc){
cbind(metadata(sc), Embeddings(sc, reduction="harmony"))
}

umapdf=function(sc, reductions=c("umap"), change.names=NULL){
  cbind(metadata(sc), lapply(1:length(reductions), function(x) Embeddings(sc, reduction=reductions[x])) %>% Reduce(cbind, .)   )

  
  }

pcadf=function(sc){
cbind(metadata(sc), Embeddings(sc, reduction="pca"))
}

getdate= function() format(Sys.time(), "%Y%b%d")


setlevels=function(so, varr, levls){
Idents(so)=varr
so@active.ident <- factor(x = so@active.ident, levels = levls)
so}

cat("checkpoint1...\n")

col2names= function(dff, coll){
rownames(dff)= dff[, coll]
dff
}

condrbind=function(x,y){

if(is.null(x)){
return(y)}else{
				if(is.null(y)){return(x)}
			  else{
			  return(rbind(x,y))
			  }}

}

condcbind=function(x,y){

if(is.null(x)){
return(y)}else{
if(is.null(y)){return(x)}
  else{
  return(cbind(x,y))
  }}

}



#recursive function to merge a bunch of clusters at the same time
# function takes several vectors c(a,b,c...), c(d,e,f)... so that clusters in each vector will be merged together into a single category. 
#output: a final vector defining new labels for the merged clusters. 

triagecells.multi= function(so, ...){
grps=list(...)# assemble groups into a list

allgroups=lapply(grps, function(x) combineclusters2(so, x))
mergername=function(x, refnames){
if(x==0){
return(refnames)
}else{
refnames[allgroups[[x]]]=LETTERS[x]
return(mergername(x-1, refnames))
}
}
out=as.character(so@meta.data[, "seurat_clusters"])
final=mergername(length(allgroups), out)
final= factor(final)
final=factor(final, labels=1:length(levels(final)))
final
}



###############
#integrate names into the dataframe
###############
namestocol=function(df){
df[, "id"]=rownames(df)
df
}

###############
#shortcut paste with "_" in between
###############

paste_=function(...) paste(..., sep="_")

###############
#open png or pdf with a timestamp (not close)
###############
timestamp= function(nm) paste0(Sys.Date(),"_", nm)


tpng=function(nm,path="./",  ...) png(paste0(path, timestamp(nm), ".png"),...)

tpdf=function(nm, path="./", ...) pdf(paste0(path, timestamp(nm), ".pdf"),...)

#############################################################
#                   
#############################################################

replacepoint= function(num) paste(strsplit(as.character(num), split='\\.')[[1]], collapse="p")


#######################################################################################
#nearest neighbor operations.
#######################################################################################
#get a particular variable from nearest neighbors of a cell. 
#######################################################################################
getnnvar= function(so, cll, vr) metadata(so)[so@graphs$RNA_nn[cll, ]==1, vr]
getsnnvar= function(so, cll, vr) metadata(so)[so@graphs$RNA_snn[cll, ]==1, vr]

getnn= function(so, cll) metadata(so)[so@graphs$RNA_snn[cll, ]==1, ]



#######################################################################################
#get data slots from a seurat object
#######################################################################################

getcounts<- function(so) so@assays$RNA@counts



#######################################################################################
# Cluster merging operations
########################################################################################
####look for cells that belong to any of the cluster ids in the list and true them. specifically for scM
combineclusters= function(lst, outvar="seurat_clusters") Reduce('|', lapply(lst, function(el) metadata(scM)[,outvar]==el))
### more general to be applied to any seurat object
combineclusters2= function(so, lst, outvar="seurat_clusters") Reduce('|', lapply(lst, function(el) metadata(so)[,outvar]==el))
####################################
#merge clusters in list lst, in factor column
####################################
renamelabel=function(column, ..., renameto){
grps=list(...)
lapply(grps, function(lst){
groupp=combineclusters(lst)
sapply(1:length(scM@meta.data[,column]), function(x) ifelse(groupp[x],renameto, as.character(scM@meta.data[x, column])), USE.NAMES=F)
})
}

cat("checkpoint2...\n")
triagecells.multi2=function(so, groupvar="seurat_clusters", ...){
  #... are the clusters in the group one after the other, ordered arbitrarily and grouped using c() when they are meant to be merged
  # warning: if some cluster groups overlap in any extent, cell indices will appear in two or more occasions and will be sequentially replaced!

  grps=list(...)
  grps2=grps %>% lapply(., function(x) paste0(addprefix(x, "^"), "$"))

original=so %>% metadata %>% pull(get(groupvar))
final.arr= rep(NA, length(original))
indlist=lapply(grps2, function(gr) grepl( gr %>% paste0(., collapse="|"), original) %>% which )
qlist=lapply(grps2, function(gr) gr %>% paste0(., collapse="|"))
vector.stamp=function( vector, ind){
  if (ind<= length(grps)){ 
 vector[indlist[[ind]]]=ind
return(vector.stamp( vector, ind+1))
 }else{
   vector[is.na(vector)]= paste_("original", original[is.na(vector)])
  return(vector) 
 }
}

out=vector.stamp(final.arr, 1)
}


###########################################
#replace each value in an array for another one.
###########################################

generalreplace=function(value, reff, target){

target[which(reff==value)]

}

##############################################
#making the row names into a gene column
##############################################

names2col=function(dff, coln="gene") {
dff[, coln]= rownames(dff) 
dff
}

##########################################################
#mohamed's function to accomodate the slingshot data into a dataframe
##########################################################

tidy_ss_output <- function(ss_output){
library(tidyr)
library(dplyr)
library(tibble)
  #cell embeddings
  cell_embedd <- ss_output@reducedDim
  #coordinates of the curve
  curve_coord <- lapply(names(ss_output@curves), function(x){
    `[[`(ss_output@curves[[x]],1) %>% 
      as.data.frame() %>% 
      rename_with(.,~paste0(x,"_",.))
    }) %>%
    bind_cols()
  #order of cell over infered curves
  cell_order <- lapply(names(ss_output@curves), function(x){
    `[[`(ss_output@curves[[x]],3) %>% 
      as.data.frame() %>% 
      rename_with(.,~paste0(x,"_lambda"))
    }) %>%
    bind_cols()
  #bind outputs in a single data frame
  cbind(cell_embedd,curve_coord,cell_order) %>% 
    as.data.frame() %>% 
    rownames_to_column("cell_id")
}


#############################################################################################################################################
#      PLOTTING FUNCTIONS
#############################################################################################################################################


#### modify hex colors up or down. 


dehash=function(color) substr(color, 2,7)
color2dec=function(color) hex2dec(substr(color, 2,7))
dec2color=function(dec) paste0('#', dec2hex(dec))
color2=function(dec) paste0('#', dec2hex(dec))

colormore=function(color, pts){

newdec=color2dec(color)+pts

if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
if(newdec<hex2dec("000000")){newdec="000000"}

return(dec2color(newdec))
}

colorless=function(color, pts){ #virtually the same as above but automatically subtracts the number so you don't have to worry about signs. 

newdec=color2dec(color)-pts

if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
if(newdec<hex2dec("000000")){newdec="000000"}
return(dec2color(newdec))
}


colormix=function(color1, color2){ #mixing two colors!
getsubcolors= function(co) lapply(c(1,3,5), function(x) substr(dehash(co), x, x+1))
  
  
newdec=color2dec(color)+colordec(color2)

if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
if(newdec<hex2dec("000000")){newdec="000000"}
return(dec2color(newdec))
}

colorsubtract=function(color1, color2){ #subtracting two colors!

newdec=abs(color2dec(color)-colordec(color2))

if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
if(newdec<hex2dec("000000")){newdec="000000"}
return(dec2color(newdec))
}



colorbleach=function(color, pts){ ###ADDS PTS TO EACH OF THE r g b SECTIONS

  
newdec=color2dec(color)+pts

if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
if(newdec<hex2dec("000000")){newdec="000000"}

return(dec2color(newdec))
}

#########################################################
#ACCESSORY FUNCTIONS
#########################################################

metadata= function(x) x@meta.data


seuratprocessing=function(o){
o=NormalizeData(o)
o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
o=ScaleData(o)
o=RunPCA(o, dims=cfig$PCA_DIMS)
o=ProjectDim(o)
#o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
o=FindNeighbors(o)
o=FindClusters(o)
o
}

SCTprocessing=function(o){
  o=NormalizeData(o)
  o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  o=ScaleData(o)
  o=RunPCA(o, dims=cfig$PCA_DIMS)
  o=ProjectDim(o)
  #o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
  o=FindNeighbors(o)
  o=FindClusters(o)
  o
}

seuratprocessing_umap=function(o){
o=NormalizeData(o)
o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
o=ScaleData(o)
o=RunPCA(o, dims=cfig$PCA_DIMS)
o=ProjectDim(o)
o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
o=FindNeighbors(o)
o=FindClusters(o)
o
}
getpname=function(x) unique(x@meta.data[, "orig.ident"])


seuratprocessing_markers1=function(so, verbose=F, umapdims=30, mincpm=3, minfc=1, fdr=0.05, selecttop=5, vfeatures=2000){
  so=SCTransform(so, verbose=verbose, variable.features.n=vfeatures)
so=RunPCA(so, verbose=verbose)
so=FindNeighbors(so, verbose=verbose)
so=FindClusters(so, verbose=verbose)
#so=RunUMAP(so, dims=1:di, verbose=verbose, return.model=T)
cat("checkpoint3...\n")

cat(paste0("finding markers using the Riffle permutation test...\n"))  

m=cb_pos_markers(counts=so@assays$RNA@counts, grouping=so %>% metadata %>% pull("seurat_clusters") )

   clusters=m %>% pull(group) %>% unique %>% as.numeric %>% sort
 ## arrange and filter best markers
    
  if(is.null(selecttop)){
  lt=lapply(clusters, function(x) m %>% filter(group==x,FDR<=fdr, logFC>=minfc, logCPM>=mincpm) %>% arrange(-logFC))
  topstr=""  
}else{
  lt=lapply(clusters, function(x) m %>% filter(group==x,FDR<=fdr, logFC>=minfc, logCPM>=mincpm) %>% arrange(-logFC) %>% head(n=selecttop))
}
    mp=lt %>% Reduce(rbind,.)

cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
so$RNA@misc$markers <- m
so$RNA@misc$top_markers <- mp
cat("Dataset", nm, ": calculating residuals for missing genes in scale data...\n")
so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
 
  gc()
  so
}


seuratmarkers.delegate=function(so, method="deseq", group_column="seurat_clusters", replicate_column=NULL, padj=0.05, minfc=0.5,selecttop=NULL,minrate=0.5, fullreload=T, fullrecreate=F, getresidual=F, min.cells.per.group=3){
  fcat(paste0("finding markers using the DElegate package...\n") ) 
#simpleCache(paste_("markers_DElegate_seuratproject",Project(so), "on_groups_of", group_column, "splitbyreps", replicate_column, "method", method), {
  
#function (object, meta_data = NULL, group_column = NULL, replicate_column = NULL, 
    #method = "edger", min_rate = 0.05, min_fc = 1, lfc_shrinkage = NULL, 
    #verbosity = 1)   
if(min.cells.per.group){
 
  countstab=so %>% metadata %>% group_by(!!sym(group_column)) %>% summarise(counts=n()) %>% as.data.frame
  
  goodcats=countstab[ countstab$counts>=min.cells.per.group, group_column]
  fcat("enough cells found for", paste(goodcats, collapse=","))
  
  so = so[, (metadata(so) %>% filter(!!sym(group_column) %in% goodcats) %>% rownames)]
   
}
  
  
 m= DElegate::FindAllMarkers2(object=so, group_column=group_column, replicate_column=replicate_column, method=method, min_rate=minrate)
 m
#m=cb_pos_markers(scwt@assays$RNA@counts, grouping=sma %>% metadata %>% pull() )
#}, assignToVar="m", reload=fullreload, recreate=fullrecreate)

#print(m)

if(group_column %in% (allcolors %>% names)){
clusters=allcolors[[group_column]][allcolors[[group_column]] %>% names %in% (so %>% metadata %>% pull(group_column) %>% unique)] %>% names
}else{
clusters=so %>% metadata %>% pull(group_column) %>% unique 
   
}

 ## arrange and filter best markers.
## currently filtering such that the minimum rate for in group cells is at least minrate. 
    pa=padj
    mi=minfc
  if(is.null(selecttop)){
  lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc))
  topstr=""  
}else{
  lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc) %>% head(n=selecttop))
}
    mp=lt %>% Reduce(rbind,.)

cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
so[["RNA"]]@misc$markers <- m
so[["RNA"]]@misc$top_markers <- mp
if(getresidual){
so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
}
so
}



################################################################################
# calculate markers using DE, each versus each other
################################################################################


seuratmarkers.delegateDE=function(so, method="deseq", group_column="seurat_clusters", compare="all_vs_all", replicate_column=NULL, padj=0.05, minfc=0.5,selecttop=NULL,minrate=0.5, fullreload=T, fullrecreate=F, getresidual=F, min.cells.per.group=3, lfc_shrinkage=NULL){
  #the other option for "compare" is each vs rest"
  fcat(paste0("finding markers using the DElegate package...\n") ) 
simpleCache(paste_("markers_DElegate_seuratproject",Project(so), "on_groups_of", group_column, "splitbyreps", replicate_column, "method", method), {
  
#function (object, meta_data = NULL, group_column = NULL, replicate_column = NULL, 
    #method = "edger", min_rate = 0.05, min_fc = 1, lfc_shrinkage = NULL, 
    #verbosity = 1)   
if(min.cells.per.group){
 
  countstab=so %>% metadata %>% group_by(!!sym(group_column)) %>% summarise(counts=n()) %>% as.data.frame
  
  goodcats=countstab[ countstab$counts>=min.cells.per.group, group_column]
  fcat("enough cells found for", paste(goodcats, collapse=","))
  
  so = so[, (metadata(so) %>% filter(!!sym(group_column) %in% goodcats) %>% rownames)]
   
}
  
  
 m= DElegate::findDE(object=so, group_column=group_column, compare=compare, replicate_column=replicate_column, method=method,lfc_shrinkage=NULL)
 m
#m=cb_pos_markers(scwt@assays$RNA@counts, grouping=sma %>% metadata %>% pull() )
}, assignToVar="m", reload=fullreload, recreate=fullrecreate)

#print(m)

if(group_column %in% (allcolors %>% names)){
clusters=allcolors[[group_column]][allcolors[[group_column]] %>% names %in% (so %>% metadata %>% pull(group_column) %>% unique)] %>% names
}else{
clusters=so %>% metadata %>% pull(group_column) %>% unique 
   
}

 ## arrange and filter best markers.
## currently filtering such that the minimum rate for in group cells is at least minrate. 
    pa=padj
    mi=minfc
  if(is.null(selecttop)){
  lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc))
  topstr=""  
}else{
  lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc) %>% head(n=selecttop))
}
    mp=lt %>% Reduce(rbind,.)

cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
so[["RNA"]]@misc$markers <- m
so[["RNA"]]@misc$top_markers <- mp
if(getresidual){
so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
}
so
}


############################################
#function to  get the sample from a barcode.
############################################
getsamplegroupfrombc=function(bc){
if (bc=="Doublet"){
return("Doublet")
}else{
if (bc=="Negative"){
return("Negative")
}else{

grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(sample_group, stage,treatment, experiment_group)][1])[1,c(1:4)]), collapse="_")
grp

}
}

}
####### get only stage from barcode
getstagefrombc=function(bc){
if (bc=="Doublet"){
return("Doublet")
}else{
if (bc=="Negative"){
return("Negative")
}else{

grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(stage)][1])[1,1]), collapse="_")
grp

}
}
}


####### get only stage from barcode
getconditionfrombc=function(bc){
if (bc=="Doublet"){
return("Doublet")
}else{
if (bc=="Negative"){
return("Negative")
}else{

grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(sample_group)][1])[1,1]), collapse="_")
grp

}
}
}

####### get only induction from barcode
gettreatmentfrombc=function(bc){
if (bc=="Doublet"){
return("Doublet")
}else{
if (bc=="Negative"){
return("Negative")
}else{

grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(treatment)][1])[1,1]), collapse="_")
grp

}
}
}


cat("checkpoint4...\n")
####replicate
getreplicatefrombc=function(bc){
if (bc=="Doublet"){
return("Doublet")
}else{
if (bc=="Negative"){
return("Negative")
}else{

grp=paste(unname(as.matrix(dA[multiseq_sequence== bc][,.(replicate)][1])[1,1]), collapse="_")
grp

}
}
}


recreateAll <- F

#dA <- fread(metapath)
#dA <- merge(dA[bsf_name!="",], fread(barcodepath)[,.(multiseq_id, multiseq_sequence)], by="multiseq_id", all.x=T)
#setkey(dA, sample_name)

#################make roc curve
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


############Quality Ccontrol QC plots function
logvarname= function(x) paste0("log10_")


qcplots=function(so, 
                 plotred=FALSE, 
                 plotsize=10, 
                 plotname="qcplots",
                 plot.title=NULL,
                 demux= NULL, 
                 topx=10000,
                 topf=20000, 
                 ww=1, 
                 hh=1,
                 colorby="seurat_clusters",
                 format="png",
                 nott=scales::comma,
                 path="./", 
                 cutoffs=NULL, 
                 cellnumbers=NULL, 
                 seriate.method=NULL,
                 plot.empty=NULL, 
                 groupvar="group1",
                 lfcvar="log_fc"){

Idents(so)=colorby  

   markers= so@assays$RNA@misc$top_markers
   
   if(is.null(markers)){
    fcat("There are no markers! please add markers under @assays$RNA@misc$top_markers" )
     
     
     
   cells=NULL
     }else{
if(is.null(seriate.method)){            
cells <- WhichCells(so, downsample = 100)
}else{
cells= seriatecells(so, clusvar=colorby, meth=seriate.method, groupvarname=groupvar, lfcvarname = lfcvar)
}
     }   
  
  
if(is.null(plot.title)){
  plot.title=Project(so)
  }

  
  require(scales)
  tsz=7
  
  
  so$log10_nCount_RNA <- log10(so$nCount_RNA)
  so$log10_nFeature_RNA <- log10(so$nFeature_RNA)

qc1=ggplot(so@meta.data, aes(x=nCount_RNA))+
  geom_histogram(binwidth=20, alpha=0.5)+
  labs(x = "Read depth",y = "Frequency")+
  guides(fill=guide_legend())+
  coord_cartesian(xlim=c(0, topf))+
  scale_y_continuous(labels = nott)+
  scale_x_continuous(labels = nott)+
  theme_classic()+
  labs(title="Reads per barcode")+theme_classic()+theme(axis.text=element_text(size=tsz), axis.title=element_text(size=tsz));

#distribution of counts per barcode
qc2=ggplot(so@meta.data, aes(x=nFeature_RNA))+
  labs(x = "Num. Genes",y = "Frequency")+
  geom_histogram(binwidth=20, alpha=0.5)+
  labs(title="Genes per barcode")+
  scale_y_continuous(labels = nott)+
  scale_x_continuous(labels = nott)+
theme_classic()+
coord_cartesian(xlim=c(0, topx))+theme_classic()+theme(axis.text=element_text(size=tsz), axis.title=element_text(size=tsz));


if(!(colorby %in% names(allcolors)) || allcolors[[colorby]] %>% length != so %>% metadata %>% pull(!!sym(colorby)) %>% unique %>% length){
  
  cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
  allcolors[[colorby]]=  randomcolors(cats %>% length)
  names(allcolors[[colorby]])= cats
}

qc3 <- FeatureScatter(so, group.by=colorby, feature1 = "nFeature_RNA", feature2 = "percent.mt", raster=T)+
  scale_color_manual(values=allcolors[[colorby]])+
  labs(x = "Num. Genes",y = "% Mitochondrial")+
  ggtitle("")+theme_classic()+
  scale_x_continuous(labels = scientific)+
  theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz),axis.title=element_text(size=tsz))+NoLegend()
qc4 <- FeatureScatter(so, group.by=colorby,feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=T)+
  scale_color_manual(values=allcolors[[colorby]])+
  labs(x = "Read depth",y = "Num.Genes")+
  ggtitle("")+
  scale_y_continuous(labels = nott)+
  scale_x_continuous(labels = nott)+
  theme_classic()+
  theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz), axis.title=element_text(size=tsz),legend.position="top",legend.title=element_blank())
qc5 <- FeatureScatter(so, group.by=colorby,feature1 = "nCount_RNA", feature2 = "percent.mt", raster=T)+
  scale_color_manual(values=allcolors[[colorby]])+
  labs(x = "Read depth",y ="% Mitochondrial")+
  ggtitle("")+theme_classic()+
  scale_x_continuous(labels = scientific)+
  theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz), axis.title=element_text(size=tsz))+NoLegend()  


if(plotred!=FALSE){

  dms=10 
     umap1=DimPlot(so, group.by = colorby, reduction=plotred, label=T, label.size=6, label.color="black", cols= allcolors[[colorby]], raster=T)+
     ggtitle(paste0(plot.title, ": ", plotred, "\n", ncol(so), " cells"))+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 

   umapmt=FeaturePlot(so, feature="percent.mt", reduction=plotred, raster=T)+
     ggtitle("% Mitochondrial")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
   
      umapreads=FeaturePlot(so,  feature="nCount_RNA",reduction=plotred, raster=T)+
     ggtitle("Reads/BC")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
      
            umapgenes=FeaturePlot(so,feature="nFeature_RNA",  reduction=plotred, raster=T)+
     ggtitle("Genes/BC")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
            
            vlnreads= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=nCount_RNA))+geom_violin()+scale_y_continuous(trans='log10')+scale_fill_manual(values=allcolors[[colorby]])+theme_classic()+NoLegend()
            vlngenes= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=nFeature_RNA))+geom_violin()+scale_fill_manual(values=allcolors[[colorby]])+theme_classic()+NoLegend()
            vlnmt= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=percent.mt))+scale_fill_manual(values=allcolors[[colorby]])+geom_violin()+theme_classic()+NoLegend()
            
            
            if( !is.null(cutoffs)){
              lcolor="red"
             vlnreads= vlnreads+geom_hline(yintercept= cutoffs$counts, color=lcolor)
             vlngenes= vlngenes+geom_hline(yintercept= cutoffs$features, color=lcolor)
             vlnmt= vlnmt+geom_hline(yintercept= cutoffs$mito, color=lcolor)
            }
            
             if( !is.null(cellnumbers)){
            ncolor="black"
               sumdata= so %>% metadata %>% group_by(seurat_clusters) %>% summarise(counts=n(), nCount_RNA=median(nCount_RNA), nFeature_RNA=median(nFeature_RNA), percent.mt=median(percent.mt))
                  
             vlnreads= vlnreads+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
             vlngenes= vlngenes+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
             vlnmt= vlnmt+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
            }
            
      if(!is.null(demux) && !is.null(so$RNA@misc$top_markers)){
        
      demuxplot=ggplot(so %>% metadata, aes(x=seurat_clusters))+geom_bar(aes(fill=!!sym(demux)), position="fill")+theme_classic()
      
            
            
             layout="AAAAMMMMM
             AAAALLLLL
                    AAAALLLLL
                    AAAALLLLL
                    AAAALLLLL
                    EEBBIIIII
                    EEBBIIIII
                    FFCCJJJJJ
                    FFCCJJJJJ
                    GGDDKKKKK
                    GGDDKKKKK
                    "
      

    
fcat("Preparing heatmap...")

#if(!is.null(plot.empty)){
 # emptydropsplot=ggplot(so %>% metadata[cells, ] %>% mutate(ord=1:length(cells)), aes(x=ord, y=1, fill=isemptydroplet))+geom_bar(color=NA)
#}else{
#emptydropsplot=ggplot(so %>% metadata[cells, ] %>% mutate(ord=1:length(cells)), aes(x=ord, y=1, fill=isemptydroplet))+geom_bar(color=NA)
#}             
             
 cat("checkpoint5...\n")                       

ph <- DoHeatmap(so, features = so$RNA@misc$top_markers$feature, group.colors=allcolors[[colorby]], slot = "scale.data", cells = cells) + NoLegend()


    
    
    qc_assembly=wrap_plots(A = umap1, B = qc3, C = qc4,
                  D = qc5, E = umapmt, F=umapreads, G=umapgenes, I=vlnmt, J=vlnreads, K=vlngenes, L=ph,M=demuxplot,
                   design = layout)
  }else{

    
  qc_assembly=((umap1)/(qc3+qc4+qc5))|(umapmt+umapreads+umapgenes)|(vlnmt+vlnreads+vlngenes)#/(qc6+qc7+qc8)
}
###################
#plotting with tsne  
###################
tpng(plotname, path=path, width=w*plotsize*ww, height=w*plotsize*hh, res=600)
print(qc_assembly)
dev.off()

}else{

#ppc1=DimPlot(object = subset(so, subset = percent.mt<10), cols= c("grey", "blue"), reduction = "pca");
###################
#plotting without tsne  
###################
  
layout <- "
####NN
AABBNN
AABBNN
CCDDEE
CCDDEE
"
  qc_assembly=wrap_plots(A = qc1, qc2 = qc2, C = qc3, D=qc4+NoLegend(), E=qc5, N=qc4,  design = layout)
  

png.pdf.prop=12
if(format=="png"){
tpng(plotname, path=path, width=w*plotsize*ww/png.pdf.prop, height=w*plotsize*hh/png.pdf.prop, res=600)
print(qc_assembly)
dev.off()
}
if(format=="pdf"){
  tpdf(plotname, path=path, width=pw*plotsize*ww/png.pdf.prop, height=pw*plotsize*hh/png.pdf.prop)
  print(qc_assembly)
  dev.off()
}

}
list(plot=qc_assembly, cells=cells)
}


os <- function() {
  require(dplyr)
  gc();
  message("Objects in MB:")
  objects = ls(envir=.GlobalEnv);
  classes = sapply(objects, function(object) { class(get(object))[1] });
  sizes = sapply(objects, function (object) { object.size(get(object)) } )
  a = data.frame(MB=sizes/1e6, class=classes)
  ord = order(sizes, decreasing=TRUE)
  a2 = a[ord,];
  a2 = a2[! a2$class == "function",]
  print(head(a2, 30))
  message("Sum: ", signif(sum(sizes/1e6),4), " MB (", signif(sum(sizes/1e9),4), " GB)")
 # fcat("total: ", a2 %>% summarise(total=sum(MB)/1000) %>% pull(total) %>% as.numeric , "GB" 
  a2
}
#os()



###########################################################################
# FUNCTION TO KEEP THE BARE BONES MINIMUM SEURAT DATA, REMOVING ALL UNNNECESSARY ASSAYS ETC.
###########################################################################

slimseurat=function(so){
  DietSeurat( so,counts = TRUE,data = TRUE,scale.data = TRUE,features = NULL,
              assays = c("RNA","refAssay"))#,dimreducs = c("pca", "umap"))
}


###data=FALSE curently not supported
slimerseurat=function(so){
  DietSeurat( so,counts = TRUE,scale.data = FALSE,features = NULL,
              assays = "RNA", dimreducs=NULL)#,dimreducs = c("pca", "umap"))
}

cat("checkpoint6...\n")
############################################################################
#  function to get a plotting table with some genes and some cells together with the metadata
############################################################################
join_meta_exp=function(so, genes, cells=NULL, normrow=FALSE, lognormrow=FALSE, metacols=NULL, assay="RNA"){
  
  if(is.null(cells)){cells= colnames(so)}
  
  m=as.data.frame(t(as.matrix(so[genes,cells ]@assays[[assay]]@data)))
  if(normrow){
    m=apply(m, 1, normvec)
  }else{
    if(lognormrow){
      m=apply(m, 1, lognormvec)
    }
    
  }
  
  if (is.null(metacols)){metacols=colnames(metadata(so))}  
  cbind( metadata(so)[cells, metacols], m)
  
}

#############################################################################
#function to attach metadata and gene reads.
#############################################################################


genesdf=function(so, genes) as.matrix((so@assays$RNA@data))[genes,  ]    
  

####### function to change the legend size inside guide
legendsize=function(x) guide_legend(override.aes = list(size = x) )

##############################################################################
#resizing overall text in the plot
resizetext=function(x) theme_classic(base_size = x)

rotatex= function(x) theme(axis.text.x = element_text(angle = x))


################################################################################
#Heatmap accessory functions
################################################################################

###function to find the names of the genes in each cluster
genesinclusters=function(phe, hvarmat,  top=NULL){
  
  lapply(row_order(phe), function(x) {
    
    if(!is.null(top)){
      rownames(hvarmat)[x][1:top]
    }else{
      rownames(hvarmat)[x]
    }
  }
  )
}

cat("checkpoint7...\n")

topnclustmembers=function(phe, top=3){
  
  Reduce(c, lapply(row_order(phe), function(x) {
    x[1:top]
    
  }
  ))
}
###############################################################

getheatmapcolors=function(ph, var, orientation='top_annotation'){
base::eval(glue::glue("pha@{orientation}@anno_list${var}@color_mapping@colors"))
colls}

#####################################################################

namenums=function(x, prefix=NULL, zero=F){
  if(is.null(prefix)){
    
    if(zero){n=1}else{n=0}
    names(x)= (1:length(x))-n; return(x)}else{
      names(x)= paste0(prefix, 1:length(x))
    }
  
  x
}

########################################################################

showpalette=function(p, pname="generic"){
  sca=10
  rr=600
  w=250
  tpng(paste_("palette", pname), res=rr, he=sca*w, wi=sca*w*2)

  pl=ggplot(1:length(p) %>% as.data.frame %>% mutate(h=1, colors=p, name=names(p)), aes(x=factor(.), y=factor(h), label=colors))+
    geom_col(aes(fill=factor(.)))+
    scale_fill_manual(values=namenums(p))+
    theme_classic()+
    coord_cartesian(ylim=c(0,1))+
    geom_text(y=0.5, color="black", angle=90)+
    #NoAxes()+
    NoLegend()+
    ggtitle(pname)
  
  if(is.null(names(p))){
 print(pl )
  }else{
 print(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))  
  }
  dev.off()
  
  return(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))
}
#showpalette(p)

##################################################################################
#GGUMAP wrapper  and derivatives
##################################################################################

ggumap=function(so, umap.df=NULL, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_", legend.size=3, color.labels.by="defaultcolor"){
  require(ggrepel)
  require(ggrastr)
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  
  ##############################################################################
  # sourcing data
  ##############################################################################
  if(!is.null(umap.df)){
     cats=umap.df %>% pull(!!sym(colorby)) %>% unique
   dat= umap.df %>% mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% mutate(none="")
  }else{
      cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
  dat=umapdf(so, reductions=reductions) %>% mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% mutate(none="")
  }
  ##############################################################################
  # preparing colorlist
  ##############################################################################
    
  # colorlist doesn't exist.
  if (is.null(colorlist[[colorby]])){
  colorlist[[colorby]]= randomcolors(length(cats)) %>% givename(., cats)
  fcat("colors used:")
  
  colorlist[[colorby]] %>% dput
  }
  
  ### actual categories are smaller than total number of colors
     if (length(colorlist[[colorby]] %>% names) > length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
       fcat("adjusting colorlist to present labels")
       
  colorlist[[colorby]]= allcolors[[colorby]][cats]
     }
    ### there are more categories than colors present
     if (length(colorlist[[colorby]] %>% names) < length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
       fcat("adjusting colorlist to present labels")
       missinglabels= setdiff(cats, colorlist[[colorby]] %>% names)
       
  colorlist[[colorby]]= c(allcolors[[colorby]][cats], randomcolors(missinglabels) %>% givename(., missinglabels))
     }
   colorlist[["defaultcolor"]]=c(one="black")
  colorlist[[color.labels.by]]=c(colorlist[[color.labels.by]] )
  colorlist[[colorby]]=c(colorlist[[colorby]] )
  ##############################################################################
  # constructing plot
  ##############################################################################
 
  dat2=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=jitter(mean(UMAP_1)), UMAP_2=jitter(mean(UMAP_2)), defaultcolor="one")
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point_rast(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]], guide = guide_legend(override.aes = list(size = legend.size, shape=15) ))+
    theme_classic()+
    #NoLegend()+
    #geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
    geom_text(data= dat2,  aes(color=!!sym(color.labels.by)), size=sz*ssca)
  
  gpuc}


cat("checkpoint8...\n")
ggumap3=function(so, umap.df=NULL,
                 colorby="seurat_clusters",
                 labelby="seurat_clusters",
                 glassworkby="seurat_clusters",
                 glasswork.params=NULL,
                 sz=0.02,
                 ssca=300,
                 colorlist=allcolors,
                 reductions=c("umap"),
                 reduction.key="UMAP_", 
                 legend.size=3){
  require(ggrepel)
require(concaveman)

  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  
  ##############################################################################
  # sourcing data
  ##############################################################################
  if(!is.null(umap.df)){
     cats=umap.df %>% pull(!!sym(colorby)) %>% unique
   dat= umap.df %>% mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% mutate(none="")
  }else{
      cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
  dat=umapdf(so, reductions=reductions) %>% mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% mutate(none="")
  }
  ##############################################################################
  # preparing colorlist
  ##############################################################################
    
  # colorlist doesn't exist.
  if (is.null(colorlist[[colorby]])){
  colorlist[[colorby]]= randomcolors(length(cats)) %>% givename(., cats)
    }
  
  ### actual categories are smaller than total number of colors
     if (length(colorlist[[colorby]] %>% names) > length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
       fcat("adjusting colorlist to present labels")
       
  colorlist[[colorby]]= allcolors[[colorby]][cats]
     }
    ### there are more categories than colors present
     if (length(colorlist[[colorby]] %>% names) < length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
       fcat("adjusting colorlist to present labels")
       missinglabels= setdiff(cats, colorlist[[colorby]] %>% names)
       
  colorlist[[colorby]]= c(allcolors[[colorby]][cats], randomcolors(missinglabels) %>% givename(., missinglabels))
     }
  
  colorlist[["defaultcolor"]]=c(one="black")
  
  ##############################################################################
  # constructing plot
  ##############################################################################

  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point_rast(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]], guide = guide_legend(override.aes = list(size = legend.size) ))+
    theme_classic()+
    #NoLegend()+
    #geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
    geom_text(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2), defaultcolor="one"), aes(color=color.labels.by), size=sz*ssca)
  
  if(glassworkby=="none" |is.na(glassworkby) |is.null(glassworkby)| glassworkby=="null" ){
    gpuc
  }else{
    
    if(!is.null(glasswork.params)){
      msg("making glasswork hulls")
    glasswork.params=list(
      size=1,
      fillvar="null",
      linecolor="#000000",
      concavity=2, 
      fillcols="#FFFFFF00",
      clus=dat %>% pull(!!sym(glassworkby)) %>% unique %>% as.character
      )
    }
      
    polys=clusterhull3(NULL, umap.df=dat, clusvar=glassworkby, clus=glasswork.params$clus, fillvar="null", linecolor=glasswork.params$linecolor, size=glasswork.params$size, reduction.key=reduction.key, fillcols=glasswork.params$fillcols, concavity=glasswork.params$concavity)
    gpuc+polys+NoLegend()
  
}}




















ggumap2=function(so, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  require(ggrepel)
  cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
  if (length(colorlist[[colorby]])!= length(cats)){
  colorlist[[colorby]]= randomcolors(length(cats))
  }
  #udimnames=c(paste0(reduction.key, 1), paste0(reduction.key, 2))
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% mutate(none="")
  gpuc=ggplot(dat, aes(x=!!sym(udimnames[1]), y=!!sym(udimnames[2]), label=!!sym(labelby)))+
    geom_point_rast(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()+
    #NoLegend()+
    #geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
  
    geom_text(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
  
  gpuc}


####return a list containng the umap and the labels separate. 
ggumapl=function(so, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  library(ggrepel)
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% rename(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% mutate(none="")
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()
    #NoLegend()+
    #geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
    
    list(gpuc, 
    geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
    )
  }

###beta
ggumapc=function(so, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  library(ggrepel)
  if (sapply(so %>% metadata %>% pull(!!sym(colorby)), is.double) %>% any){
    colorf=scale_color_manual(values=colorlist[[colorby]])
  }else{colorf=scale_color_gradient(values=colorlist[[colorby]])}
    
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% rename(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2))
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()+
    NoLegend()+
    geom_text_repel(data=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=mean(UMAP_1), UMAP_2=mean(UMAP_2)), size=sz*ssca)
  
  gpuc}

cat("checkpoint9...\n")

ggumaph=function(so, colorby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  library(ggrepel)
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% rename(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% mutate(none="", ccondition=factor(ccondition, levels=names(colorlist[["ccondition"]])))
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    gghighlight()+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()
    #NoLegend()
  
  gpuc}

ggumaps=function(so, colorby="seurat_clusters", sz=0.02, ssca=300, colorlist=allcolors, reductions=c("umap"), reduction.key="UMAP_"){
  library(ggrepel)
  dim1=paste0(reduction.key, 1)
  dim2=paste0(reduction.key, 2)
  dat=umapdf(so, reductions=reductions) %>% rename(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% mutate(none="", ccondition=factor(stage, levels=names(colorlist[["stage"]])))
  gpuc=ggplot(dat, aes(x=UMAP_1, y=UMAP_2))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    gghighlight()+
    scale_color_manual(values=colorlist[[colorby]])+
    theme_classic()
  
  gpuc}
#################################

addprefix= (function(cc, prefix="C")  paste0(prefix, cc)) %>% Vectorize(., USE.NAMES=F)

givename=
  function(x, nam) {
    names(x)= nam
    x
  } 



givecolnames=function(x, ind=NULL, nms){
  if (is.null(ind)){
   ind= 1:length(nms) 
  }
 colnames(x)[ind]= nms
 x
  
}
 
giverownames=function(x, nms){
 rownames(x)= nms
 x
  
}
 
#####

findgeneindex=function(so, x)   (so %>% rownames) %>% grepl(paste0("^", x,"$"), .) %>% which


########################


    na2zero= function(x) generalreplace(is.na(x), c( FALSE, TRUE), c(x, 0))



maptoref= function(query, refdataset, mapvar, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NULL, mappinglabel=paste0(querylabel,"_to_",reflabel)){
  ###if query or ref do not have an SCT refAssay, compute them
  nameflag=0
  if ((intersect(query %>% colnames, refdataset %>% colnames) %>% length)>0){
    colnames(query)= paste0("query_", colnames(query))
    
    colnames(refdataset)=paste0("ref_", colnames(refdataset))
    nameflag=1
  }             
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)="refAssay"
  DefaultAssay(refdataset)="refAssay"
  refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F)
  refdataset <- RunUMAP(refdataset, reduction = "pca", dims = pdms, , assay="refAssay", return.model=TRUE, verbose=F)
  ref <- AzimuthReference(
    object = refdataset,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "refAssay",
    metadata = mapvar,
    dims = pdms,
    k.param = 31,
    reference.version = "1.0.0"
  )
  cat("Finding transfer anchors...\n\n\n")
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "refAssay",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = intersect(rownames(x = ref), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = 20,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = ref,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = 20,
    store.weights = TRUE
  )
  
  
  cat("Integrating embeddings...\n\n\n")
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = ref,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    dims.to.integrate=pdms
  )
  
  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(ref[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = ref[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = ref[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  cat("checkpoint10...\n")
  if(!is.na(customcol)){
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}



maptoref2= function(query, refdataset, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  nameflag=0
  if ((intersect(query %>% colnames, refdataset %>% colnames) %>% length)>0){
    colnames(query)= paste0("query_", colnames(query))
    
    colnames(refdataset)=paste0("ref_", colnames(refdataset))
    nameflag=1
  }             
  
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
  
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before running."))
    }else{
    
  DefaultAssay(refdataset)=ref.assay  
  refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
  refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
  }
  }

  cat("Finding transfer anchors...\n\n\n")
  

  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", #
    reference.neighbors = NULL,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
  
  
  cat("Integrating embeddings...\n\n\n")
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = refdataset,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    dims.to.integrate=pdms,
    new.reduction.name="integrated_dr"
  )
  
  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(refdataset[["mappingpca"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = refdataset[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = refdataset[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.na(customcol)){
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}

cat("checkpoint11...\n")

maptoref.pt1= function(query, 
                       refdataset, 
                       ref.assay="refAssay",
                       query.assay="refAssay",
                       mapvar, 
                       normalization.method="SCT",
                       reference.pca="pca", 
                       reference.umap="umap",
                       reference.reduction=reference.pca,
                       reference.neighbors=NULL,
                       project.query=FALSE,
                       n.trees=20,
                       pdms=1:50,
                       return.merge.metadata=T, 
                       return.query.metadata=F,
                       return.lean=T,
                       querylabel="query", 
                       reflabel="reference", 
                       customcol=NULL, 
                       mappinglabel=paste0(querylabel,"_to_",reflabel),
                       process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  
  nameflag=0
  if ((intersect(colnames(query), colnames(refdataset)) %>% length )>0){
    msg("labelling cell ids according to dataset of origin...")
    query=RenameCells(query, new.name= paste0("query_", colnames(query)))
    
    refdataset=RenameCells(refdataset, new.name= paste0("ref_", colnames(refdataset)))
    nameflag=1
  }             
  fcat("finished labeling cell ids")
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
  
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before reprocessing."))
    }else{
    
  DefaultAssay(refdataset)=ref.assay  
  refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
  refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
  }
  }else{
    fcat("setting", reference.pca, "and", reference.umap, "as mapping reductions")
    refdataset@reductions$mappingpca= refdataset@reductions[[reference.pca]]
    refdataset@reductions$mappingumap= refdataset@reductions[[reference.umap]]
  }

   # if(!is.null(reference.neighbors){
  #  nnname=paste0(umapred, "_snn")
  #  fcat("Attempting to use", nnname, "as reference neighbors...")
   # reference.neighbors=nnname
  #}
  
  if(is.null(reference.reduction)){
    fcat("Setting", reference.pca, "as the reference reduction")
  reference.reduction=reference.pca
    }
  
  
  cat("Finding transfer anchors...\n\n\n")
  

  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", # method of projection
    reference.neighbors = reference.neighbors,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
list(query=query, anchors=anchors, ref=refdataset)

      
      
    }
     
    
  


maptoref.pt2= function(query, refdataset, anchors, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap",reference.reduction=reference.pca, reference.neighbors=NULL, project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NULL, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){

  
    query=Seurat::MapQuery(
  anchor=anchors,
  query=query,
  reference=refdataset,
  refdata = mapvar,
  new.reduction.name = "mapfun.pca",
  reference.reduction = reference.reduction,
  reference.dims = pdms,
  query.dims = pdms,
  reduction.model = reference.umap,
  transferdata.args = list(),
  integrateembeddings.args = list(),
  projectumap.args = list(
    reduction.name=reference.umap,
    reduction.key= prepare.rk(reference.umap)
    ),
  verbose = TRUE
)   


    
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )    
    
    
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refumap1=prepare.rk(reference.umap) %>% paste0(., "1")
  refumap2=prepare.rk(reference.umap) %>% paste0(., "2")
  
  #refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, reference.umap)[, refumap1], col.name=refumap1)
  query=AddMetaData(query, reddf(query, reference.umap)[, refumap2], col.name=refumap2)
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
   if(!is.null(customcol)){
    getcols=function(x) x[, c("origin", refumap1, refumap2,nm1, nm2, nm3, customcol)]
   }else{
     getcols=function(x) x[, c("origin", refumap1, refumap2,nm1, nm2, nm3)  ]
   }
   if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      #if (nameflag==1){
       # colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      #}
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }   
    
    
  
  
}

cat("checkpoint13...\n")

maptoref.workflow=function(query, refdataset, ...){
#function(query, refdataset, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F)
  mappedcells=maptoref.pt1(query, refdataset, ...)
mapmetadata=maptoref.pt2(query=mappedcells$query, anchors=mappedcells$anchors, refdataset=mappedcells$ref, ...)
  mapmetadata
  
}


labeltransfer= function(query, 
                       refdataset, 
                       ref.assay="refAssay",
                       query.assay="refAssay",
                       mapvar, 
                       normalization.method="SCT",
                       reference.pca="pca", 
                       reference.umap="umap",
                       reference.reduction=reference.pca,
                       reference.neighbors=NULL,
                       project.query=FALSE,
                       n.trees=20,
                       pdms=1:50,
                       return.lean=T,
                       querylabel="query", 
                       reflabel="reference", 
                       customcol=NULL, 
                       mappinglabel=paste0(querylabel,"_to_",reflabel),
                       process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  
  nameflag=0
  if ((intersect(colnames(query), colnames(refdataset)) %>% length )>0){
    msg("labelling cell ids according to dataset of origin...")
    query=RenameCells(query, new.name= paste0("query_", colnames(query)))
    
    refdataset=RenameCells(refdataset, new.name= paste0("ref_", colnames(refdataset)))
    nameflag=1
  }             
  fcat("finished labeling cell ids")
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
  
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before reprocessing."))
    }else{
    
  DefaultAssay(refdataset)=ref.assay  
  refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
  refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
  }
  }else{
    fcat("setting", reference.pca, "and", reference.umap, "as mapping reductions")
    refdataset@reductions$mappingpca= refdataset@reductions[[reference.pca]]
    refdataset@reductions$mappingumap= refdataset@reductions[[reference.umap]]
  }

   # if(!is.null(reference.neighbors){
  #  nnname=paste0(umapred, "_snn")
  #  fcat("Attempting to use", nnname, "as reference neighbors...")
   # reference.neighbors=nnname
  #}
  
  if(is.null(reference.reduction)){
    fcat("Setting", reference.pca, "as the reference reduction")
  reference.reduction=reference.pca
    }
  
  
  cat("Finding transfer anchors...\n\n\n")
  

  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", # method of projection
    reference.neighbors = reference.neighbors,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
#list(query=query, anchors=anchors, ref=refdataset)

  
    cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )    
    

  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  
  #refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
   if(!is.null(customcol)){
    getcols=function(x) x[, c("origin", nm1, nm2, nm3, customcol)]
   }else{
     getcols=function(x) x[, c("origin",nm1, nm2, nm3)  ]
   }
   

      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
  
  }   
    
  










maptoref.pt3= function(query, refdataset, anchors, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){

  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(refdataset[[reference.pca]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  

  
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = refdataset[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = refdataset[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.null(customcol)){
    getcols=function(x) x[, c("origin",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}

cat("checkpoint14...\n")
#######################################################################################################
#Function to generate hulls from umap clusters that light up according to some metric
#######################################################################################################
clusterhullfill=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar=NULL, fillcols=NULL, size=0.05){
  # if squares...
  #maxx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% max
  #maxy=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>%pull(UMAP_2) %>% max
  #minx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% min
  #miny=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_2) %>% min
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf
    
  }
  
  if(is.null(clus)){
    clus=umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }

  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    clhull=cl[chull(cl), ]
    
    if(is.null(pols)){
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      pols<-append(pols, geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }
    
    if(cc==length(clus)){
      return(pols)
    }else{
      return(catpols(cc+1, pols))}
    
    
  }
  
  catpols(1)
}



################################################################################
#
################################################################################

clusterhullfill2=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar=NULL, fillcols=NULL, size=0.05, concavity=2){
  # if squares...
  #maxx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% max
  #maxy=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>%pull(UMAP_2) %>% max
  #minx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% min
  #miny=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_2) %>% min
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf
    
  }
  
  if(is.null(clus)){
    clus=umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }

  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ]
    clhull=concaveman(cl, concavity=concavity)
    
    if(is.null(pols)){
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      pols<-append(pols, geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }
    
    if(cc==length(clus)){
      return(pols)
    }else{
      return(catpols(cc+1, pols))}
    
    
  }
  
  catpols(1)
}

#############################################################################
#Code for seurat utilities which are useful
#############################################################################
NoLegend=function (...)
{
    no.legend.theme <- theme(legend.position = "none", validate = TRUE,
        ...)
    return(no.legend.theme)
}


NoAxes=function (..., keep.text = FALSE, keep.ticks = FALSE)
{
    blank <- element_blank()
    no.axes.theme <- theme(axis.line.x = blank, axis.line.y = blank,
        validate = TRUE, ...)
    if (!keep.text) {
        no.axes.theme <- no.axes.theme + theme(axis.text.x = blank,
            axis.text.y = blank, axis.title.x = blank, axis.title.y = blank,
            validate = TRUE, ...)
    }
    if (!keep.ticks) {
        no.axes.theme <- no.axes.theme + theme(axis.ticks.x = blank,
            axis.ticks.y = blank, validate = TRUE, ...)
    }
    return(no.axes.theme)
}



#######################importing manually functions to run hull sample from mvgps

hull_sample=function (X, num_grid_pts = 500, grid_type = "regular", 
          trim_hull = FALSE, trim_quantile = NULL) 
{
  X_rslt <- X_check(X)
  assign("X", X_rslt$X)
  assign("m", X_rslt$m)
  grid_type <- match.arg(grid_type, choices = c("regular", 
                                                "random", "hexagonal"))
  if (trim_hull == TRUE) {
    if (is.null(trim_quantile)) 
      stop("trim_hull set to TRUE but trim_quantile not specified.", 
           call. = FALSE)
    if (trim_quantile < 0.5 | trim_quantile > 1) 
      stop("trim_quantile must be between [0.5, 1]", 
           call. = FALSE)
    trim_upper <- apply(X, 2, quantile, trim_quantile)
    trim_lower <- apply(X, 2, quantile, 1 - trim_quantile)
    X_trim <- sapply(seq_len(m), function(x) {
      ifelse(X[, x] > trim_upper[x], NA, ifelse(X[, x] < 
                                                  trim_lower[x], NA, X[, x]))
    })
    colnames(X_trim) <- colnames(X)
    X <- na.omit(X_trim)
  }
  if (m == 2) {
    hpts <- chull(X)
    hpts <- c(hpts, hpts[1])
    hpts_vs <- as.matrix(X[hpts, ])
    m <- Polygon(hpts_vs)
    ps <- Polygons(list(m), 1)
    sps <- SpatialPolygons(list(ps))
    sp_grid_pts <- spsample(sps, n = num_grid_pts, type = grid_type)
    grid_pts <- coordinates(sp_grid_pts)
    colnames(grid_pts) <- colnames(X)
  }
  else {
    hpts <- geometry::convhulln(X)
    hpts_ind <- unique(c(hpts))
    hpts_vs <- X[hpts_ind, ]
    grid_pts <- NULL
  }
  return(list(hpts_vs = hpts_vs, grid_pts = grid_pts, X = X))
}


cat("checkpoint15...\n")
X_check=function (X) 
{
  X <- as.matrix(X)
  m <- ncol(X)
  if (!any(apply(X, 2, is.numeric))) 
    stop("X must be numeric", call. = FALSE)
  if (m < 2) 
    stop("Exposure is not multivariate", call. = FALSE)
  return(list(X = X, m = m))
}


####################################################################
#function to make cluster hulls filled with points that signal a correspondence of cells, with points jittered to minimise overplotting effects
####################################################################


clusterhullpoint=function(so, umap.df= NULL, query.df, clusvar, clus, linecolor="white" ,fillvar=NULL, colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_", concavity=2, downsample=1, rescale.cells=NULL){
  require(sf)

  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf %>% mutate(UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
    
  }else{
   umap.df= umap.df  %>% mutate(UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
  }
  
  
  destinyclusters=query.df %>% group_by( !!sym(mapfun(clusvar))) %>% summarise(counts=n()) %>% as.data.frame %>% arrange(-counts)
  
  ##collect all the destinies of each cell in each  together
  
  if(is.null(rescale.cells)){
  des=destinyclusters %>% mutate(temptot= sum(counts), frac=counts/temptot, newcounts= ceiling(frac*temptot) ) %>% col2names(., mapfun(clusvar))
  }else{
    des=destinyclusters %>% mutate(temptot= sum(counts), newtot=!!rescale.cells, frac=counts/temptot, newcounts= ceiling(frac*newtot) ) %>% col2names(., mapfun(clusvar))

  }
  
  
  fcat("printing normal des:")
  print(destinyclusters)
  if(downsample!=1){
   fcat("printing downsampled des:")
    des= des %>% mutate(counts=ceiling(counts*downsample))
    print(des)
  }
    getclusterfrac=function(x) des[x, "frac"]
  umap.df[["clusterfrac"]]=NA
  umap.df[["clusterfrac"]]=sapply(umap.df %>% pull(clusvar), getclusterfrac) %>% unname

  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% dplyr::select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ] #obsolete but keep in case we need to revert to convex hulls
    
    
    
     clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
     
     rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
     clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% mutate(null="one")
     
   
      dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
     ld= nrow(dists)
     lc= nrow(clhull)
     inds=lapply(rownames(clhull), function(r){ names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]}) %>% unlist
     clhull[[fillvar]]= cl[inds, ][[fillvar]]
     clhull[[mapfun(varr)]]= cl[inds, ][[mapfun(varr)]]
     clhull[[varr]]= cl[inds, ][[mapfun(varr)]]
       
    
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #at this point downsampling should already have been taken care of
    numcells=des %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% pull(newcounts)
    originalcells=des %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% pull(counts)
    if(length(numcells)==0){
      numcells=0
    }
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    subsample=(numcells*downsample) %>% round
    #fcat("numcells is" , numcells)
    if(!is.null(colorpoints) & !is.na(numcells) & numcells>2){
    #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
    
      polygon = sf::st_polygon(list(clhull[, c("UMAP_1","UMAP_2")] %>% data.matrix))
      
     # Sample 50 random points within the polygon
      hsamples = sf::st_sample(polygon, size=numcells)
      hsamples= sf::st_coordinates(hsamples) %>% givecolnames(., nms=c("UMAP_1","UMAP_2"))
                           
    
    fcat("samples computed...")
    fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
    #print(hsamples)
    ###making matrix of mock coordinates for each cell
    
    if(subsample==originalcells){
      rplce=F}else{rplce=T}
    pointsdf=query.df %>% mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% sample_n(subsample, replace=rplce) %>% mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2]) 

    pointsdf[[varr]] = pointsdf[[mapfun(varr)]]
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size),
                   geom_point_rast(data=pointsdf, inherit.aes=T, size=pointsize,aes(x=UMAP_1, y=UMAP_2, color=!!sym(colorpoints))))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size),
        geom_point_rast(data=pointsdf,inherit.aes=T,size=pointsize, aes(x=UMAP_1, y=UMAP_2, color=!!sym(colorpoints)))))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(append(pols, scale_color_manual(values=allcolors[[colorpoints]])))
      #return(pols)
    }else{
      return(catpols(cc+1, pols))}
    
    }else{###if there are no query cells in cluster
      
      
      if(is.null(pols)){
        fcat("creating poygon list...")
        pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
      }else{
        fcat("appending poygon to list...")
        pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
      }
      
      if(cc==length(clus)){
        #msg("finishing up. Attaching color scheme...")
        fcat("finishing up. Attaching color scheme...")
        return(append(pols, scale_color_manual(values=allcolors[[colorpoints]])))
      }else{
        return(catpols(cc+1, pols))}
      
      
      
      
      
      
      }
    
    
     
  }
  
  catpols(1)
}



####just make cluster hulls
clusterhull=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_"){
  require(sp)

  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf     
  }
  
 umap.df= umap.df %>% mutate(null=NA, UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))

  
    if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
    }
 
 if( is.null(allcolors[[fillvar]])){
  
   allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
 }

  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    clhull=cl[chull(cl), ]
    
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    #fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #numcells=0
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    
    
    #if(!is.null(colorpoints) && numcells>2){
    #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
    
    #fcat("samples computed...")
    #fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
    #print(hsamples)
    ###making matrix of mock coordinates for each cell
    #pointsdf=query.df %>% mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2])
    
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                    , scale_fill_manual(values=allcolors[[fillvar]]) 
                    ) %>% Reduce(append, .)
             )
      #return(pols)
    }else{
      return(catpols(cc+1, pols))}
  }
  catpols(1)
}

################################################################################
#uses concaveman for cluster hulls
################################################################################

clusterhull2=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_", concavity=2){
  require(sp)

  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf     
  }
  
 umap.df= umap.df %>% mutate(null=NA, UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))

  
    if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
    }
 
 if( is.null(allcolors[[fillvar]])){
  
   allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
 }

  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ]
    
     clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
     
     rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
     clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% mutate(null="one")
     # fillvar on the concave matrix gets assiged the color of the closest point
     
     

       
     if(fillvar != "null"){
      dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
     ld= nrow(dists)
     lc= nrow(clhull)
     clhull[[fillvar]]= cl[lapply(rownames(clhull), function(r){
       names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]  
       
       }) %>% unlist, ][[fillvar]]
       
       }

     
 
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    #fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #numcells=0
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    
    
    #if(!is.null(colorpoints) && numcells>2){
    #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
    
    #fcat("samples computed...")
    #fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
    #print(hsamples)
    ###making matrix of mock coordinates for each cell
    #pointsdf=query.df %>% mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2])
    
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                    , scale_fill_manual(values=allcolors[[fillvar]]) 
                    ) %>% Reduce(append, .)
             )
      #return(pols)
    }else{
      return(catpols(cc+1, pols))}
  }
  catpols(1)
}


################################################################################
# non recursive version of clusterhull2 to avoid weird R issues
################################################################################


clusterhull3=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_", concavity=2){
  require(sp)

  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf     
  }
  
 umap.df= umap.df %>% mutate(null=NA, UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))

  
    if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
    }
 
 if( is.null(allcolors[[fillvar]])){
  
   allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
 }

 pols=NULL
polygon=function (x, y = NULL, density = NULL, angle = 45, border = NULL, 
    col = NA, lty = par("lty"), ..., fillOddEven = FALSE) 
{
    ..debug.hatch <- FALSE
    xy <- xy.coords(x, y, setLab = FALSE)
    if (is.numeric(density) && all(is.na(density) | density < 
        0)) 
        density <- NULL
    if (!is.null(angle) && !is.null(density)) {
        polygon.onehatch <- function(x, y, x0, y0, xd, yd, ..debug.hatch = FALSE, 
            ...) {
            if (..debug.hatch) {
                points(x0, y0)
                arrows(x0, y0, x0 + xd, y0 + yd)
            }
            halfplane <- as.integer(xd * (y - y0) - yd * (x - 
                x0) <= 0)
            cross <- halfplane[-1L] - halfplane[-length(halfplane)]
            does.cross <- cross != 0
            if (!any(does.cross)) 
                return()
            x1 <- x[-length(x)][does.cross]
            y1 <- y[-length(y)][does.cross]
            x2 <- x[-1L][does.cross]
            y2 <- y[-1L][does.cross]
            t <- (((x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - 
                x1))/(xd * (y2 - y1) - yd * (x2 - x1)))
            o <- order(t)
            tsort <- t[o]
            crossings <- cumsum(cross[does.cross][o])
            if (fillOddEven) 
                crossings <- crossings%%2
            drawline <- crossings != 0
            lx <- x0 + xd * tsort
            ly <- y0 + yd * tsort
            lx1 <- lx[-length(lx)][drawline]
            ly1 <- ly[-length(ly)][drawline]
            lx2 <- lx[-1L][drawline]
            ly2 <- ly[-1L][drawline]
            segments(lx1, ly1, lx2, ly2, ...)
        }
        polygon.fullhatch <- function(x, y, density, angle, ..debug.hatch = FALSE, 
            ...) {
            x <- c(x, x[1L])
            y <- c(y, y[1L])
            angle <- angle%%180
            if (par("xlog") || par("ylog")) {
                warning("cannot hatch with logarithmic scale active")
                return()
            }
            usr <- par("usr")
            pin <- par("pin")
            upi <- c(usr[2L] - usr[1L], usr[4L] - usr[3L])/pin
            if (upi[1L] < 0) 
                angle <- 180 - angle
            if (upi[2L] < 0) 
                angle <- 180 - angle
            upi <- abs(upi)
            xd <- cos(angle/180 * pi) * upi[1L]
            yd <- sin(angle/180 * pi) * upi[2L]
            if (angle < 45 || angle > 135) {
                if (angle < 45) {
                  first.x <- max(x)
                  last.x <- min(x)
                }
                else {
                  first.x <- min(x)
                  last.x <- max(x)
                }
                y.shift <- upi[2L]/density/abs(cos(angle/180 * 
                  pi))
                x0 <- 0
                y0 <- floor((min(y) - first.x * yd/xd)/y.shift) * 
                  y.shift
                y.end <- max(y) - last.x * yd/xd
                while (y0 < y.end) {
                  polygon.onehatch(x, y, x0, y0, xd, yd, ..debug.hatch = ..debug.hatch, 
                    ...)
                  y0 <- y0 + y.shift
                }
            }
            else {
                if (angle < 90) {
                  first.y <- max(y)
                  last.y <- min(y)
                }
                else {
                  first.y <- min(y)
                  last.y <- max(y)
                }
                x.shift <- upi[1L]/density/abs(sin(angle/180 * 
                  pi))
                x0 <- floor((min(x) - first.y * xd/yd)/x.shift) * 
                  x.shift
                y0 <- 0
                x.end <- max(x) - last.y * xd/yd
                while (x0 < x.end) {
                  polygon.onehatch(x, y, x0, y0, xd, yd, ..debug.hatch = ..debug.hatch, 
                    ...)
                  x0 <- x0 + x.shift
                }
            }
        }
        if (missing(col) || is.null(col)) {
            col <- par("fg")
        }
        else if (any(is.na(col))) {
            col[is.na(col)] <- par("fg")
        }
        if (is.null(border)) 
            border <- col
        if (is.logical(border)) {
            if (!is.na(border) && border) 
                border <- col
            else border <- NA
        }
        start <- 1
        ends <- c(seq_along(xy$x)[is.na(xy$x) | is.na(xy$y)], 
            length(xy$x) + 1)
        num.polygons <- length(ends)
        col <- rep_len(col, num.polygons)
        if (length(border)) 
            border <- rep_len(border, num.polygons)
        if (length(lty)) 
            lty <- rep_len(lty, num.polygons)
        if (length(density)) 
            density <- rep_len(density, num.polygons)
        angle <- rep_len(angle, num.polygons)
        i <- 1L
        for (end in ends) {
            if (end > start) {
                if (is.null(density) || is.na(density[i]) || 
                  density[i] < 0) 
                  .External.graphics(C_polygon, xy$x[start:(end - 
                    1)], xy$y[start:(end - 1)], col[i], NA, lty[i], 
                    ...)
                else if (density[i] > 0) {
                  polygon.fullhatch(xy$x[start:(end - 1)], xy$y[start:(end - 
                    1)], col = col[i], lty = lty[i], density = density[i], 
                    angle = angle[i], ..debug.hatch = ..debug.hatch, 
                    ...)
                }
                i <- i + 1
            }
            start <- end + 1
        }
        .External.graphics(C_polygon, xy$x, xy$y, NA, border, 
            lty, ...)
    }
    else {
        if (is.logical(border)) {
            if (!is.na(border) && border) 
                border <- par("fg")
            else border <- NA
        }
        .External.graphics(C_polygon, xy$x, xy$y, col, border, 
            lty, ...)
    }
    invisible()
}
  for(cc in 1:length(clus)){
  
    cl=umap.df %>% mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% dplyr::select(UMAP_1, UMAP_2, !!sym(fillvar))

    
     clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
     
     rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
     clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% mutate(null="one")
     # fillvar on the concave matrix gets assiged the color of the closest point
     
     if(fillvar != "null"){
      dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
     ld= nrow(dists)
     lc= nrow(clhull)
     clhull[[fillvar]]= cl[lapply(rownames(clhull), function(r){
       names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]  
       
       }) %>% unlist, ][[fillvar]]
       
       }

    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                    , scale_fill_manual(values=allcolors[[fillvar]]) 
                    ) %>% Reduce(append, .)
             )
      #return(pols)
    }
 } 
    return(pols) 
}




randomcolors= function(n)
  
  lapply(1:n, function(x) paste0("#",sample(c(0:9,LETTERS[1:6]), 1),
                                sample(c(0:9,LETTERS[1:6]), 1),
                                sample(c(0:9,LETTERS[1:6]), 1),
                                sample(c(0:9,LETTERS[1:6]), 1),
                                sample(c(0:9,LETTERS[1:6]), 1),
                                sample(c(0:9,LETTERS[1:6]), 1))) %>% Reduce(c, .)




demapfun=function(x) strsplit(x, split="mapfun_")[[1]][2]
mapfun=function(x) x %>% addprefix(., "mapfun_")


#####################################

prefixids.smart= function(so){
  colnames(so@assays$RNA@counts)= paste_(Project(so),  colnames(so@assays$RNA@counts))
  so
}
prefixids= function(so, pref=""){
  if ( !("patientID" %in% (so %>% metadata %>% colnames)) ){
    if ( ("orig.ident" %in% (so %>% metadata %>% colnames)) ){
  RenameCells(so, new.names= so %>% metadata %>% mutate(newcellid=paste0(pref,orig.ident,"_", so %>% colnames)) %>% pull(newcellid) )
    }else{
      RenameCells(so, new.names= so %>% metadata %>% mutate(newcellid=paste0(pref, so %>% colnames)) %>% pull(newcellid) )
    }}else{
    RenameCells(so, new.names= so %>% metadata %>% mutate(newcellid=paste0(pref, patientID, "_", so %>% colnames)) %>% pull(newcellid) )
  }
  
}

##assuming there is a patientID
prefixids= function(so, pref=""){
      RenameCells(so, new.names= so %>% metadata %>% mutate(newcellid=paste0(pref, so %>% colnames)) %>% pull(newcellid) )
  
}


prefixids.matrix= function(mat, prefix){
  colnames(mat)= paste_(prefix,  colnames(mat))
  mat
}

prefixids.meta= function(mat, prefix){
  rownames(mat)= paste_(prefix,  rownames(mat))
  mat
}

alph=77
applyalpha=function(color) color %>% paste0(., alph) 


################################################################################
# SPECIFIC NCNB2 functions
################################################################################

# split the hto_demux name into different bits. 


extractname = Vectorize(function(x) {
  #cat(x, "\n")
  x=gsub("-", "_", x)
  if (!grepl("_", x)) {
    return(c(NA, NA, NA, NA) %>% givename(., c("cell.line", "condition", "stage", "replicate")))
  } else{
    #HTO-H7_17q_day3_rep1_CMO
    
    #small replacement on the condition section of the string
    #cat("1")
     if (!grepl("^H7", x)) {
     x0=paste0("H7_", x)
       }else{x0=x}

    x1 = gsub("H7WT", "H7_WT", x0)
    x2 = gsub("17q1qMYCN_", "17qx_", x1) #transform all the potential17q1qs into 17qx
    x3 = gsub("17qMYCN_", "17q1qMYCNdox_", x2) #transform all the potential mycns into mycndox
    x4 = gsub("MYCNdox_", "MYCN_", x3) #transform all dox into mycn
    x5 = gsub("17qx", "17q1q", x4) #tranform all qx into 17q1q


    x6 = gsub("_CMO", "", x5)
    x7= gsub("rep","R", x6)
    x8= gsub("day", "D", x7)
  
    xf=x8

    #make a list with all the parts thereafter
    xs1 = strsplit(xf, split = "_")[[1]]
    #print(x2)
    names(xs1) = c("cell.line", "condition", "stage", "replicate")
    xs2 = as.list(xs1)
    xs2$condition = paste0("c", xs2$condition)
    
    xs2 %>% unlist
  }
  
}, USE.NAMES=F)

################################################################################
# screen print shortcuts
################################################################################


#concatenate a message for the screen, add screen format

fcat= function(...) cat(paste(...,sep=" "), "\n")
msg=fcat
showdataset=function(so)
fcat("Dataset", so$dsname %>% unique %>% paste(., collapse="-")) 


################################################################################
# dataset filtering on the go
################################################################################

filterds=function(so, dbcutoff=NULL, filter.empty.drops=NULL, demuxvar="sample_name", qc.pars=NULL){
  if(!is.null(qc.pars)){
   qcpars=qc.pars 
  }
  
  nm=Project(so)
  if(!(nm %in% names(qcpars))){
    st=paste(nm, "does not have an assigned QC filter")
    stop(st)
  }else{
    
  fcat("Updating filter metadata...")  
  so=AddMetaData(so, so %>% metadata %>% mutate(counts.cutoff=qcpars[[nm]]$counts, features.cutoff= qcpars[[nm]]$features, mito.cutoff=qcpars[[nm]]$mito, qcc.pass=((nFeature_RNA>=features.cutoff & nCount_RNA >= counts.cutoff & percent.mt<= mito.cutoff )))  %>% dplyr::select(counts.cutoff, features.cutoff, mito.cutoff, qcc.pass))
  }
  if(demuxvar %in% (so %>% metadata %>% colnames)){
    so=so[, so %>% metadata %>% filter(!(!!sym(demuxvar) %in% c("Negative", "Doublet", "Multiplet", "Unassigned", "Blank", NA)), nCount_RNA>=qcpars[[nm]]$counts, nFeature_RNA>=qcpars[[nm]]$features, percent.mt<=qcpars[[nm]]$mito ) %>% rownames] 

  }else{
  so=so[, so %>% metadata %>% filter(nCount_RNA>=qcpars[[nm]]$counts, nFeature_RNA>=qcpars[[nm]]$features, percent.mt<=qcpars[[nm]]$mito ) %>% rownames] 
  }  

  if( !is.null(dbcutoff) & ("scdblscore" %in% (so %>% metadata %>% colnames))){
  so=so[, so %>% metadata %>% filter(scdblscore<=dbcutoff) %>% rownames]
  so=AddMetaData(so, metadata= so %>% metadata %>% mutate(scDblFinder.cutoff=dbcutoff) %>% select(scDblFinder.cutoff))
    }
  if( !is.null(filter.empty.drops) & ("isemptydroplet" %in% (so %>% metadata %>% colnames))){
  so=so[, so %>% metadata %>% filter(isemptydroplet==FALSE) %>% rownames]
    }

 
so   
}


filteremptydrops=function(so, filter.empty.drops=NULL){
  nm=Project(so)

  if( !is.null(filter.empty.drops) & ("isemptydroplet" %in% (so %>% metadata %>% colnames))){
  so=so[, so %>% metadata %>% filter(isemptydroplet==FALSE) %>% rownames]
    }

 
so   
}



stagename= Vectorize(function(nm) stages[[nm]], USE.NAMES=F)


################################################################################
#reorder cells in clusters according to different methods (e.g., seriation or signature, i.e. how trong the markers are)
################################################################################

cellinfo.getcells=function(so, clusvar="seurat_clusters",ncells=100, clusters=NULL, fullrecreate=T, replicatevar=NULL, fullreload=F, ...){ #replicatevar=NULL,  clusters=NULL){
 
  Idents(so)=clusvar
if(is.null(clusters)){
clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
}else{
 clusids=clusters 
}
   
  fcat("Getting cells for each cluster...")
cellist=lapply(as.character(clusids), function(x) {
  out=WhichCells(so, idents=x, downsample=ncells)
  #fcat("total number of cells:", length(out))
  out
  }
  
  ) %>% givename(., as.character(clusids))

cellist
}


cellinfo.getmarkers=function(so, clusvar="seurat_clusters", selecttop=NULL){
  
if(is.null(so$RNA@misc$top_markers)){  
so=seuratmarkers.delegate(so, minfc=1, group_column=clusvar, replicate_column=replicatevar, fullrecreate=fullrecreate, fullreload=fullreload, ...)
}
  
}


  


seriatecells=function(so, clusvar="seurat_clusters", clusters=NULL, meth="signature", markers=NULL, ncells=100, ass="SCT", groupvarname="group1",lfcvarname="log_fc", extended.output=T, deduped=T, nmarkers=NULL){
  
Idents(so)=clusvar
if(is.null(clusters)){
clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
}else{
 clusids=clusters 
}

  
   markers= so@assays$RNA@misc$top_markers
   
   
   if(is.null(markers)){
    stop("There are no markers! please add markers under @assays$RNA@misc$top_markers\n" ) 
   }else{
    fcat("Splitting groups...")
    markerlist0=lapply(clusids, function(x) markers %>% filter(!!sym(groupvarname)==x)) %>% givename(., clusids)
    
        hasgenes= Vectorize(function(tabb) as.logical(nrow(tabb)>0), USE.NAMES=F)
        
        
    fcat("current length of markerlist0", length(markerlist0))
    which.have.genes=hasgenes(markerlist0)
    markerlist0=markerlist0[which.have.genes]
    fcat("current length of markerlist0 after filtering", length(markerlist0))
    fcat("these clusters have genes:", paste(markerlist0 %>% names, collapse=" ")) 
    
    fcat("current length of clusids", length(clusids))
    clusids=as.character(clusids)[which.have.genes]
    fcat("current length of clusids after filter", length(clusids))
    names(markerlist0)= as.character(clusids) 
    
    
    if(deduped==T){
    fullgns= lapply(markerlist0, function(x) x %>% pull(feature)) %>% Reduce(c, .)

###when duplicates, keep the gene as marker in the cluster where it shows highest log_fold change
fcat("Assigning genes in more than one cluster to the highest expressing cluster...")    
deduped=lapply(unique(fullgns), function(x){
  lapply(markerlist0, function(y) y %>% filter(feature==x)) %>% Reduce(condrbind, .) %>% arrange(-!!sym(lfcvarname)) %>% head(n=1)  %>% select(feature,!!sym(lfcvarname), !!sym(groupvarname)) 
  
}) %>% Reduce(rbind, .)

deduped2= deduped %>% mutate(ncluster=sapply(deduped[[groupvarname]], function(x) generalreplace(x, clusids,1:length(clusids)), USE.NAMES=F )) %>% arrange(ncluster)

#arrange genes from highest log fold change to lowest per cluster
markers.revised= deduped2 %>% group_by(ncluster) %>% arrange(-!!sym(lfcvarname), .by_group=T) %>% ungroup()


    }else{
      markers.revised=markers # forwarding the original table
      }
    
    
markergaps= markers.revised %>% pull(ncluster) %>% as.integer %>% diff %>% as.logical %>% which
markervector=markers.revised %>% pull(feature) 

if(is.null(nmarkers)){    
markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature))
}else{
  markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature) %>% head(nmarkers))

  }
    names(markerlist)=as.character(clusids)  
    markergaps=lapply(1:length(markerlist), function(x) rep(x, length(markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
    
   }
   
  
    
fcat("Getting cells for each cluster...")
cellist=lapply(as.character(clusids), function(x) {
  out=WhichCells(so, idents=x, downsample=ncells)
  #fcat("total number of cells:", length(out))
  out
  }
  
  ) %>% givename(., as.character(clusids))


ordcellist=tryCatch({

if (meth=="rankcounts"){
  celltotals=list()
  ordcellist=lapply(1:length(cellist), function(x) {
    
    newcells= metadata(so)[cellist[[x]], ] %>% select(!!sym(clusvar), nCount_RNA) %>% arrange(-nCount_RNA) %>% slice_head(n=ncells.show) %>% rownames
    celltotals[[x]]<<- rep(x, length(newcells))
    newcells
  }
  )
  

}else{
  if (meth=="signature"){
    fcat("Using signature method to rank cells...")
    celltotals=list()
    ordcellist=lapply(1:length(cellist), function(x) {
      
      print(markerlist[[x]])
      
      sumcells=apply(so[markerlist[[x]],cellist[[x]]]@assays[[ass]]@counts, 2, sum)
      
      
      #arranged from weakest to strongest  
      newcells= metadata(so)[cellist[[x]], ] %>% mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    ) %>% givename(., clusids)
    
  }else{
    if (meth=="seriate"){
      library(seriation)
      celltotals=list()
      ordcellist=lapply(1:length(cellist), function(x) {
        
        
        ss=seriation::get_order(seriate(as.matrix(so[markerlist[[x]],cellist[[x]]]@assays$RNA@data), margin=2, method="BEA_TSP"))
        celltotals[[x]]<<- rep(x, length(ss))
               cellist[[x]][ss]
      }
      )  %>% givename(., clusids)
      

    }else{
     
      if(meth=="seurat"){
          fcat("Using Seurat Module Score method to rank cells...")
    celltotals=list()
    
    newnames=names(markerlist) %>% make.names
    so <- AddModuleScore(so, markerlist, name=newnames )
    
    
    ordcellist=lapply(1:length(clusids), function(x) {
      
      print(markerlist[[x]])
      
      modulename= paste0(newnames[x], x)
      
      #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
      
      newcells=so@meta.data[cellist[[x]], ] %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
      
      #arranged from weakest to strongest  
      #newcells= metadata(so)[cellist[[x]], ] %>% mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    ) %>% givename(., clusids)
    
          
        }else{
          
          if(meth=="seurat2"){
          fcat("Using Seurat Module Score method to rank cells and retrieveng strongest", ncells, "cells in clusters...")
    celltotals=list()
    
    so <- AddModuleScore(so, markerlist, name=names(markerlist) )
    
    
    ordcellist=lapply(1:length(clusids), function(x) {
      
      print(markerlist[[x]])
      
      modulename= paste0(clusids[x], x)
      
      #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
      newcells=so@meta.data %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
      newcells=newcells[1:ncells]
      
      #arranged from weakest to strongest  
      #newcells= metadata(so)[cellist[[x]], ] %>% mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    )
    
      
          
        }
          
        }
      
      
    }}}
  ordcellist} , error=function(e){fcat("There was trouble with seriating the cells. returning unseriated cells");
    print(e); 
    cellist
    })
  
  celltotals=lapply(1:length(ordcellist), function(x) rep(x, length(ordcellist[[x]])))
  
  fcat("making cell vector...")
    cellvector=Reduce(c, ordcellist)
    fcat("getting the gap positions...")
    cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which  

cell.reference=lapply(1:length(ordcellist), function(nn) rep(names(ordcellist)[nn], length(ordcellist[[nn]])   )) %>% Reduce(c, .)  
marker.reference=lapply(1:length(markerlist), function(nn) rep(names(markerlist)[nn], length(markerlist[[nn]])   )) %>% Reduce(c, .)

cell.annotation= as.data.frame(cell.reference) %>% giverownames(cellvector) %>% givecolnames(., 1, clusvar)
marker.annotation= as.data.frame(marker.reference) %>% giverownames(markervector) %>% givecolnames(., 1, clusvar)


    

if(extended.output==T){
  list(cells=cellvector,
       cell.list=ordcellist,
       gaps.cells=cellgaps,
       gaps.markers=markergaps,
       markers=markervector,
       markerlist=markerlist,
       cell.annotation=cell.annotation,
       marker.annotation=marker.annotation,
       clustering.variable=clusvar,
       input.cluster.ids=clusids,
       method=meth,
       cell.metadata=clusvar,
       marker.metadata=clusvar
       )}else{
cellvector
  }
}



seriatecells2=function(so, clusvar="seurat_clusters", clusters=NULL, meth="signature", markers=NULL, ncells=100, ass="SCT", groupvarname="group",lfcvarname="logFC", extended.output=F, deduped=T, groups.match.markers=T){
  
Idents(so)=clusvar
if(is.null(clusters)){
clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
}else{
 clusids=clusters 
}
markers.are.list=F
markers.are.table=F
## check if markers have been provided externally
  if(!is.null(markers) & class(markers)=="list"){
    markers.are.list=T
    markers.are.table=F
    fcat("External marker list has been provided.")
  }else{
   markers.are.list=F
   if(!is.null(markers) & (class(markers) %in% c("data.table", "matrix","data.frame"))){
   fcat("External marker table has been provided.")
   markers.are.table=T
     }
  }
### CASE WHEN EXTERNAL LIST OF MARKERS IS PROVIDED
if(markers.are.list){    
    markerlist0=markers
    
    if(!is.null(names(markerlist0))){# if markerlist has names
    clusids=names(markerlist0)
    
    #all(clusters %in% (ref.cellinfo$markerlist %>% names))
        
    }else{ ### if markerlist has no names just make some
    fcat("external markers have no names. Automatically generating names...")
      clusids=(1:length(markerlist0)) %>% paste0("C", .)
      names(markerlist0)= clusids
      
    }}

# assemble the markers into a simple dataframe to perform duplication operations
if(deduped==T){
  fcat("deduplicating marker list...")
fullgns=lapply(1:length(markerlist0), function(x) markerlist0[[x]]) %>% Reduce(c, .)
fullgns.groups=lapply(1:length(markerlist0), function(x) rep(names(markerlist0)[x], length(markerlist0[[x]])) ) %>% Reduce(c, .)

# in the absence of GE information, remove duplicates in the order they appear first time marker appears is considered the 
#most relevant. 
dupd= duplicated(fullgns)
markers.df=as.data.frame(fullgns.groups[!dupd]) %>% giverownames(fullgns[!dupd]) %>% givecolnames(., 1, nms=clusvarname) 

markerlist= lapply(as.character(clusids), function(x)  markers.df %>% filter(!!sym(clusvarname)==!!x) %>% rownames)
names(markerlist)=as.character(clusids)
}

### NO EXTERNAL MARKERS PROVIDED: MARKER TABLE EXPECTED WITHIN OBJECT
if(!markers.are.list & !markers.are.table & is.null(markers)){
   #if external markers are not provided at all, start from the top markers table
   markers= so@assays$RNA@misc$top_markers
   if(is.null(markers)){
     markers= so@assays$RNA@misc$markers
     
     if(is.null(markers)){
     #if top_markers is not there, then throw an error
    stop("There are no markers! please add markers under @assays$RNA@misc$top_markers\n" ) 
     }else{
      markers.are.table=T 
     }
       
   }else{
     markers.are.table=T
   }
   
   if(markers.are.table){
   # if there is a table, collect the DE markers for each group. 
    fcat("Splitting groups and retrieving top markers from table...")
    markerlist0=lapply(clusids, function(x) markers %>% filter(!!sym(groupvarname)==x))
    
    hasgenes= Vectorize(function(tabb) as.logical(nrow(tabb)>0), USE.NAMES=F)
    fcat("current length of markerlist0", length(markerlist0))
    which.have.genes=hasgenes(markerlist0)
    markerlist0=markerlist0[which.have.genes]
    fcat("current length of markerlist0 after filtering", length(markerlist0))
     fcat("current length of clusids", length(clusids))
    clusids=as.character(clusids)[which.have.genes]
    fcat("current length of clusids after filter", length(clusids))
    names(markerlist0)= clusids 
   
    
    fcat("Markers found for the following groups:", paste(clusids, collapse=","))
   
   } 
   
   #removing duplicates when a marker DE table is provided (not a list)
    if(deduped==T & markers.are.table){
    fullgns= lapply(markerlist0, function(x) x %>% pull(feature)) %>% Reduce(c, .)

###when duplicates, keep the gene as marker in the cluster where it shows highest log_fold change
fcat("Assigning genes in more than one cluster to the highest expressing cluster...")    
deduped=lapply(unique(fullgns), function(x){
  lapply(markerlist0, function(y) y %>% filter(feature==x)) %>% Reduce(condrbind, .) %>% arrange(-!!sym(lfcvarname)) %>% head(n=1)  %>% select(feature,!!sym(lfcvarname), !!sym(groupvarname)) 
  
}) %>% Reduce(rbind, .)

deduped2= deduped %>% mutate(ncluster=sapply(deduped[[groupvarname]], function(x) generalreplace(x, clusids,1:length(clusids)), USE.NAMES=F )) %>% arrange(ncluster)

#arrange genes from highest log fold change to lowest per cluster
markers.revised= deduped2 %>% group_by(ncluster) %>% arrange(-!!sym(lfcvarname), .by_group=T) %>% ungroup()


    }
    
    
markergaps= markers.revised %>% pull(ncluster) %>% as.integer %>% diff %>% as.logical %>% which
#markervector=markers.revised %>% pull(feature) 
    
markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature))
    
    names(markerlist)=as.character(clusids)  
        
   }
  
### COMMON DOWNSTREAM PROCESSING AFTER DEFINING THE MARKER LIST 

markervector= markerlist %>% Reduce(c, .)
 markergaps=lapply(1:length(markerlist), function(x) rep(x, length(markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
 
    
fcat("Getting cells for each cluster...")
cellclasses=so %>% metadata %>% pull(clusvar) %>% unique
fcat("cellclasses are", cellclasses)
cellist=lapply(cellclasses, function(x) {
  
#  if(!(x %in% ( cellclasses ))  ){
  #  NULL}else{
  out=WhichCells(so, idents=x, downsample=ncells)
  #fcat("total number of cells:", length(out))
  out
#}
}
  
  ) %>% givename(., cellclasses)

#fcat("removing clusters for which no cells were found")
#nonnulls=lapply(cellist, function(x) is.null(x)) %>% unlist %>% unname
#cellist= cellist[!nonnulls] %>% givename(., as.character(clusids[!nonnulls]))


if (meth=="rankcounts"){
  celltotals=list()
  ordcellist=lapply(1:length(cellist), function(x) {
    
    newcells= metadata(so)[cellist[[x]], ] %>% select(!!sym(clusvar), nCount_RNA) %>% arrange(-nCount_RNA) %>% slice_head(n=ncells.show) %>% rownames
    celltotals[[x]]<<- rep(x, length(newcells))
    newcells
  }
  )
  

}else{
  if (meth=="signature"){
    fcat("Using signature method to rank cells...")
    celltotals=list()
    ordcellist=lapply(1:length(cellist), function(x) {
      
      print(markerlist[[x]])
      
      sumcells=apply(so[markerlist[[x]],cellist[[x]]]@assays[[ass]]@counts, 2, sum)
      
      
      #arranged from weakest to strongest  
      newcells= metadata(so)[cellist[[x]], ] %>% mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    )
    
    cellvector=Reduce(c, ordcellist)
    fcat("getting the gap positions...")
    cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
  }else{
    if (meth=="seriate"){
      library(seriation)
      celltotals=list()
      ordcellist=lapply(1:length(cellist), function(x) {
        
        
        ss=seriation::get_order(seriate(as.matrix(so[markerlist[[x]],cellist[[x]]]@assays$RNA@data), margin=2, method="BEA_TSP"))
        celltotals[[x]]<<- rep(x, length(ss))
               cellist[[x]][ss]
      }
      ) 
      

    }else{
     
      if(meth=="seurat"){
          fcat("Using Seurat Module Score method to rank cells...")
    celltotals=list()
    
    newnames=names(markerlist) %>% make.names
    #names(newnames)=names(markerlist)
    fcat("newnames is ", newnames)
    so=tryCatch({ AddModuleScore(so, markerlist, name=newnames)}, error=function(e){
      fcat("trying 12 bins because of numerical issues...")
      AddModuleScore(so, markerlist, nbin=12, name=newnames)})
    
    
    if(groups.match.markers){
    ordcellist=lapply(1:length(markerlist), function(x) {
      
      print(names(markerlist)[x])
      
      modulename= paste0(newnames[x], x)
      markername=names(markerlist)[x]
      #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
      newcells=so@meta.data[cellist[[markername]], ] %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
      if(!is.null(ncells) & ncells<= length(newcells)){
        #get the ncells strongest cells
      newcells=newcells[(length(newcells)-ncells):length(newcells)]
      }
      #arranged from weakest to strongest  
      #newcells= metadata(so)[cellist[[x]], ] %>% mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    )
    
    cellvector=Reduce(c, ordcellist)
    fcat("getting the gap positions...")
    cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
    }else{ # if cell groups do not match markers, we just retrieve the ncells with the strongest signature
    fcat("cells and markers do not have corresponding ids. Getting cells with strongest signature")
         ordcellist=lapply(1:length(clusids), function(x) { #we evaluate the strength of all signatures regardless of the cell label
      
      print(names(markerlist)[x])
      
      modulename= paste0(newnames[clusids[x]], x)
      
      #arrange from biggest to smallest
      newcells=so@meta.data %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
      
           if(!is.null(ncells) & ncells<= length(newcells)){
        #get the ncells strongest cells
      newcells=newcells[(length(newcells)-ncells):length(newcells)]
      }
      #arranged from weakest to strongest  
      #newcells= metadata(so)[cellist[[x]], ] %>% mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    )
    
    cellvector=Reduce(c, ordcellist)
    fcat("getting the gap positions...")
    cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
      
      } 
        }else{ #another method different from above
          
          if(meth=="seurat2"){
          fcat("Using Seurat Module Score method to rank cells and retrieveng strongest", ncells, "cells in clusters...")
    celltotals=list()
    
    so <- AddModuleScore(so, markerlist, name=names(markerlist) )
    
    
    ordcellist=lapply(1:length(clusids), function(x) {
      
      print(markerlist[[x]])
      
      modulename= paste0(clusids[x], x)
      
      #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
      newcells=so@meta.data %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
      newcells=newcells[1:ncells]
      
      #arranged from weakest to strongest  
      #newcells= metadata(so)[cellist[[x]], ] %>% mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    )
    
    cellvector=Reduce(c, ordcellist)
    fcat("getting the gap positions...")
    cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
      
          
        }
          
        }
      
      
    }
  
    }
}
if(extended.output==T){
  list(cells=cellvector,cell.list=ordcellist , gaps.cells=cellgaps, gaps.markers=markergaps, markers=markervector, markerlist=markerlist)}else{
cellvector
  }
}

################################################################################
# classify cells into singlet, doublet or negative based on the hto_demux
################################################################################

classifydemux= Vectorize(function(x){
  if(x!="Negative" & x!= "Doublet" & x!= "Unassigned"& x!= "Blank" & x!= "Multiplet"){
  return("Singlet")
    }else{
    return(x)
      }
}, USE.NAMES=F)


################################################################################
# retrieve a doublet scoring filter for a dataset
################################################################################

getdbcutoff= Vectorize(function(nm){
 dbcutoffs[[nm]]
   
}, USE.NAMES=F)


################################################################################
# paste two strings together with a comma in between (great with Reduce())
################################################################################
  

pastec=function(a,b) paste0(c(a,b), collapse=",")


################################################################################
# a more comprehensive condition taking into account the demux state, negative or doublet. 
################################################################################

get_condition_demux=Vectorize(function(condition, demux){
  
  if(is.na(condition)){
    
    if(!(is.na(demux))){
    return(demux)
    }else{
    return("Blank")
      
      }
      
  }else{
  return(condition)
    }
  
}, USE.NAMES=F)


################################################################################
# remove the last character from a string
################################################################################
chopstring=function(st) substr(st, 1, nchar(st)-1)


################################################################################
#fix bad names in 
################################################################################

fixnames= function(x){nms= x %>% colnames %>% make.names; colnames(x)=nms; fcat("new names are", nms); x}


################################################################################
#get the replicate
################################################################################
getrep= Vectorize(function(x) repmap[[x]], USE.NAMES=F)


################################################################################
# trim final one in the cell names
###############################################################################

trimone= Vectorize(function(x) strsplit(x, split="-")[[1]][1]) 


################################################################################
#
################################################################################

is.loaded= function(pkgname) !requireNamespace(pkgname, quietly = TRUE)

################################################################################
# make square legend
################################################################################

squarelegend = function(varr, size=1, shape=0, colorlist=allcolors){
  
  arglist=list()
  arglist[[varr]]= guide_legend(
      override.aes = list(
        shape = shape,  # 0: square, 1: circle
        size = size,   # Customize point sizes
        color = allcolors[[varr]],  # Customize point colors
        label = allcolors[[varr]] %>% names # Customize legend labels
      )
    )
    
 do.call(guides, arglist)

}

################################################################################
# compare category changes for cells across several category variables
################################################################################

transferplot= function(so=NULL, df=NULL,  labelvars, colorlist=allcolors, facetvar="all_cells", facet.ncol=NULL, facet.nrow=NULL, column.width=0.5, alpha=0.5, scale=5){

  
sca=scale  
freqvar=facetvar
library(ggalluvial)
library(ggrepel)

#prepare colors if unavailable
  

#lapply(labelvars, function(x) allcolors[[k]])    
if(!is.null(df)){
  met=df %>% mutate(all_cells="All cells")
}else{
if(!is.null(so)){  
met= so %>% metadata %>% mutate(all_cells="All cells")  
}else{
error("Please add either a Seurat object or a data frame as input!")
}}


fcat("preparing colors")
for(k in 1:length(labelvars)){
  
  if(is.null(colorlist[[labelvars[k]]])){
    fcat(labelvars[k])
    cats=met %>% pull(!!sym(labelvars[k])) %>% unique
    fcat("unique categories:", paste(cats, collapse=" "))
   colorlist[[labelvars[k]]]<- randomcolors( cats %>% length) %>% givename(., cats)
   #allcolors[[labelvars[k]]]<<- colorlist[[labelvars[k]]]
     
  }
}
fcat("Colors used for plotting:")
cat(code.for.list.of.vectors(colorlist[labelvars]))



# factors suspended for now
for(k in 1:length(labelvars)){
 met = met %>% mutate(!!labelvars[k]:=factor(!!sym(labelvars[k]), levels=colorlist[[labelvars[k]]] %>% names ) )
   
}


t1=met %>% dplyr::select(c(freqvar,labelvars)) %>% group_by(!!!syms(c(freqvar, labelvars))) %>% summarise( counts=n())  %>% group_by(!!sym(freqvar)) %>% mutate(Freq=counts/sum(counts)) %>% ungroup() 
t2= t1 %>% mutate(ind= 1:nrow(t1)) %>% pivot_longer(., cols=labelvars, names_to="originalcolumn", values_to="cell.label")

tpdf(path=params$plotpath, paste0("alluvialplot_labeltransfer", paste(labelvars, collapse="-")), width=pw*5/2*sca/10, height=pw*sca*3/14)
ga=ggplot(t2,
          aes(y = Freq, x=factor(originalcolumn, levels=labelvars), stratum=cell.label, alluvium=ind, fill=cell.label))+ 
  geom_flow() +
  geom_stratum(width = column.width, alpha=alpha) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)))+#, min.segment.length=0.1) +
  scale_x_discrete(expand= c(.1, .1), labels = labelvars)+
  facet_wrap(~factor(get(freqvar), levels=colorlist[[freqvar]] %>% names), ncol=facet.ncol, nrow=facet.nrow)+
  scale_fill_manual(values=colorlist[labelvars] %>% Reduce(c,.))+
  theme_classic()
print(ga+NoLegend())
dev.off()

ga+NoLegend()
}

################################################################################
# get highest loadings from pca
################################################################################

gethighloadings=function(so, comps=1:5, nfeatures=10, reduction="pca", reduction.key="PC_"){
  rk=reduction.key
  compnms= paste0(rk, comps)
  flist= lapply(compnms, function(x){
    
    list(positive=
  so@reductions[[reduction]]@feature.loadings %>% as.data.frame %>% select(!!sym(x)) %>% arrange(-!!sym(x)) %>% head(nfeatures) %>% rownames
  , negative=
    so@reductions[[reduction]]@feature.loadings %>% as.data.frame %>% select(!!sym(x)) %>% arrange(!!sym(x)) %>% head(nfeatures) %>% rownames
    ) 
    
  }) %>% givename(compnms)
    flist
}
  

################################################################################

################################################################################
code.for.list.of.vectors <- function(list_of_vectors) {
  code0 <- "list(\n"
  codelist=list()
  nms=names(list_of_vectors)
  # Loop through each vector in the list
  for (i in seq_along(list_of_vectors)) {
    
    if(is.null(nms)){
    vector_name <- paste0("vector", i)
    }else{
     vector_name <- nms[i] 
    }
    code <-  paste0("  ", vector_name, " = c(\n")
    
    # Loop through each element in the vector
    for (element_name in names(list_of_vectors[[i]])) {
      element_value <- toString(list_of_vectors[[i]][element_name])
      code <- paste0(code, paste0("   ", element_name, " = '", element_value, "',\n"))
    }
    
    code <- paste0(code, ")\n")
    codelist[[i]]=code
  }
  code2= codelist %>% Reduce(c, .) %>% paste(., collapse=",")
  
  code3 <- paste0(code0, code2, ")\n") 
  
  return(code3 %>% gsub(",\n)", ")", .))
}

# Example usage:
#vector1 <- c("name1" = "value1", "name2" = "value2")
#vector2 <- c("name3" = "value3", "name4" = "value4")
#list_of_vectors <- list(vector1, vector2)

#code_string <- code.for.list.of.vectors(allcolors[labelvars])
#code_string
#cat(code_string)

################################################################################
# remove a field from a list
################################################################################

dropfield= function(lisst, x){ lisst[[x]]=NULL; lisst} 


################################################################################
# create an arbitrary colormap
################################################################################

breakColors = function(breaks, colors, center=0, tol=0.001)
{
    ## In case of explicit color definitions
    nbreaks = length(breaks)
    nclass  = nbreaks - 1
    if (!is.function(colors)) {
        ncolors = length(colors)
        if (ncolors > nclass) {
            warning("more colors than classes: ignoring ", ncolors-nclass, " last colors")
            colors = colors[1:nclass]
        } else if (nclass > ncolors) {
            stop(nclass-ncolors, " more classes than colors defined")
        }
    } else {
        ## Are the classes symmetric and of same lengths?
        clens = diff(breaks)
        aclen = mean(clens)
        if (aclen==0) stop("Dude, your breaks are seriously fucked up!")
        relerr = max((clens-aclen)/aclen)
        if ( (center %in% breaks) & (relerr < tol) ) { ## yes, symmetric
            ndxcen = which(breaks==center)
            kneg = ndxcen -1 
            kpos = nbreaks - ndxcen
            kmax = max(kneg, kpos)
            colors = colors(2*kmax)
            if (kneg < kpos) {
                colors = colors[ (kpos-kneg+1) : (2*kmax) ]
            } else if (kneg > kpos) {
                colors = colors[ 1 : (2*kmax - (kneg-kpos)) ]
            }
        } else {                                      ## no, not symmetric
            colors = colors(nclass)
        }
    }
    colors
}

################################################################################

removegene=function(genes,x) genes[genes!=x]

################################################################################
# format recuciton key based on an arbitrary id.
################################################################################


prepare.rk=function(x)  paste0(x, "_") %>% gsub("\\.", "", .)


################################################################################
#
################################################################################


sling_cell_df <- function(sling_ctr){
  raw_pseudotime <- slingPseudotime(sling_ctr)
  order_df <- apply(raw_pseudotime, 2, function(x){
    x[!is.na(x)] <- rank(x[!is.na(x)])
    x
  }) %>%
    as.data.frame()%>% 
    rename_with(.,~paste0("order","_",.)) %>% 
    mutate(mean_order = apply(., 1, mean, na.rm = TRUE))
  
  
  pseudotime_df <- apply(raw_pseudotime, 2, function(x){
    x1 <- x[!is.na(x)]
    x[!is.na(x)] <- (x1-min(x1))/(max(x1)-min(x1))*100
    x
  })%>%
    as.data.frame()%>% 
    rename_with(.,~paste0("pseudotime","_",.))%>% 
    mutate(mean_pseudotime = apply(., 1, mean, na.rm = TRUE))
  
  df <- bind_cols(raw_pseudotime,pseudotime_df, order_df) %>% 
    as.data.frame()
  
  return(df)
}


################################################################################
# ENTROPY FUNCTION
################################################################################

#https://rpubs.com/philjet/shannonentropy
#entropy function
entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}



################################################################################
# findmin
################################################################################

findmin = function(arr){
 
m= min(abs(arr))
return(which(arr==m))

}
  

################################################################################
# get stage number from a "DX" format
################################################################################
  
get.stage.number= Vectorize(
  function(x) substr(x, 2, 10) %>% as.numeric
)

################################################################################
#
################################################################################

prepare.rk=function(x)  paste0(x, "_") %>% gsub("\\.", "", .)



################################################################################
# dplyr - ready get rows
################################################################################


getrows=function(x, rows) x %>% names2col(., "nm") %>% filter(nm %in% rows) %>% select(-nm)


################################################################################
# generalreplace, grep version
################################################################################

#value may be longer than elements in reff and so may contain reff
generalreplace.grepl=function(vec, refs, targets, return.sums=T){
 subss=list()
 finalvector=vec
  for(n in 1:length(refs)){
    reff=refs[n]
    which.match=grepl(reff, vec)
    subss[[n]]=rep(FALSE, length(vec))
    subss[[n]][which.match]=TRUE
    finalvector=gsub(reff, targets[n], finalvector )
    
  }
  allsums=subss %>% Reduce('+', .)  
 if(any(allsums)>1){
 warning("Same value substituted multiple times")
 }
if(return.sums==T){
   list(final.vector=finalvector, number.replacements=allsums)
}else{
 finalvector 
}
  
  }
################################################################################
#
################################################################################

      reorderdf= function(df, varr, catorder, return.list=T){
        
        l=lapply(catorder, function(x){
          df %>% filter(!!sym(varr)==x)
          })
        if(return.list){
        }else{
        l=l %>% Reduce(rbind,.)
        }
        l
      }
    
    



################################################################################
# wrapper to take 2 seurat objects with a mapfun variable and show the points of one in another's clusters (groups)
################################################################################

#project.to.glasswork=function(ref, query, )
################################################################################
# functions that emulate matlab zeros and ones
################################################################################
zeros <- function(rows, cols) {
  return(matrix(0, nrow = rows, ncol = cols))
}

ones <- function(rows, cols) {
  return(matrix(1, nrow = rows, ncol = cols))
}


remove.zerorows= function(matt) matt[apply(matt, 1, sum)!=0,]
get.variant.rows= function(matt){

  #non.zero.rows=apply(matt, 1, sum)!=0
  variant.rows=apply(matt, 1, sd)!=0
  #print(cbind(non.zero.rows, variant.rows))
   variant.rows

  }


get.nonzerorows= function(matt) (apply(matt, 1, sum)!=0) %>% unname  
  

fill.missing.markers.matrix=function(matt, markernames, fill.value=0){
missing.ones=NULL 
  
  newlines=lapply(1:length(markernames), function(k){ 

    nl=NULL
    if(!(markernames[k] %in% rownames(matt))){
    
      
      
      #missing.ones=c(missing.ones, markernames[k])
      
      nc=length(colnames(matt))
      nl=matrix(data=0, nrow=1, ncol= nc, dimnames=list(markernames[k], colnames(matt)))
    }
  nl
  }) %>% Reduce(rbind, .)
  #fcat("number of columns of matt", length(colnames(matt)), "number of columns of newlines", ncol(newlines))
  matt=rbind(matt, newlines)
#if(!is.null(missing.ones)){
 
  #fcat("missing genes filled:", paste(missing.ones, collapse=",")) 
#}
matt
}



################################################################################
#for plotting purposes, slightly jitter one position of a row
################################################################################


random.glitch.row=function(mat, rownumber=x, val=0.00001){

  s=sample(1: ncol(mat), 1)
  mat[  rownumber,  s]=mat[  rownumber,  s]+val
  
  mat
  
}

################################################################################
# randomly glitch all invariant rows for plotting purposes
################################################################################

glitch.invariants=function(matt, val=0.00001){
 
invariant.rows= !get.variant.rows(matt)   


for(k in 1:length(invariant.rows)){
set.seed(invariant.rows[k])
  matt=random.glitch.row(matt, rownumber=invariant.rows[k], val=val)
  
    
}
set.seed(42)
matt
}
  
################################################################################
#create a heatmap based directly from a cellinfo object
################################################################################
  

cellinfo.heatmap= function(so,
                           cellinfo,
                           assay="RNA",
                           name=NULL,
                           pheatmap.params=NULL, ## currently not used
                           genes.to.label=NULL, 
                           color.labels.by=NULL,
                           colorlist=NULL,
                           cmap=NULL,
                           colorbar.max=2,
                           colorbar.spacing=0.1,
                           show.legend=TRUE,
                           return.everything=F,
                           return.cellinfo=F){

#function refresh cell info to update ubject
  

  
# heatmap colors

if(is.null(cmap)){
  cmap=c("#0B9988", "#f7f7f7", "#F9AD03") 
}

colorpalette=colorRampPalette(cmap, space="rgb")

####Heatmap color schemes
#petrol gold

br=seq(-colorbar.max,colorbar.max,colorbar.spacing)  ## numeric breaks of the color bins
colls=colorpalette(length(br))   

#cellinfo= refresh.cellinfo(cellinfo)





fcat("getting data matrix...")
mat=so[cellinfo$markers, cellinfo$cells ]@assays[[assay]]@data %>% as.matrix
fcat("finding invariant and absent genes...")
absent.genes= setdiff(cellinfo$markers, rownames(mat))
variant.rows=get.variant.rows(mat)
invariant.genes=mat[!variant.rows, ] %>% rownames 
fcat("Warning: the following invariant genes have been removed:\n", invariant.genes %>% dput)
cellinfo2=removemarkers(cellinfo, markers=union(invariant.genes, absent.genes))

if(!is.null(cellinfo$cell.metadata)){
fcat("adding cell annotation from variables in metadata...")
cellinfo2=cellinfo2 %>% add.cell.annotation(so, ., cellinfo2$cell.metadata)
}



if(!is.null(genes.to.label)){
  fcat("finding genes to label...")
 # original markers
gns=cellinfo2$markers


geneguide=  1:length(gns)
names(geneguide)=gns
#mrkrs are genes of interest



if(!is.null(color.labels.by) & !is.null(colorlist)){

gene.order.info=geneguide[genes.to.label]

  
getcol=Vectorize(function(x) 
  if(x %in% (allcolors[[color.labels.by]] %>% names)){
  allcolors[[color.labels.by]][x]}else{
    NA
    }, USE.NAMES=F)

getposition=Vectorize(function(x) geneguide[x], USE.NAMES=F)

colref= (cellinfo2$marker.annotation %>% names2col(., "gene") %>% filter(gene %in% genes.to.label) %>% mutate(color=getcol(!!sym(color.labels.by)), position= getposition(gene) )  )
genecols=colref %>% arrange(position) %>% pull(color) %>% unname

  
linesgp= grid::gpar(col=genecols)
labelsgp= grid::gpar(col=genecols)
}else{
  linesgp=NULL
  labelsgp=NULL
  }
ha = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = geneguide[genes.to.label] %>% unname %>% as.numeric, 
                                                                   labels = genes.to.label,
                                                                   lines_gp= linesgp, 
                                                                   labels_gp= labelsgp))


#geneguide[genes.to.label] %>% as.character %>% dput
   
}else{ha=NULL}

fcat("getting data matrix...")

mat2=so[cellinfo2$markers, cellinfo2$cells ]@assays[[assay]]@data %>% as.matrix

cats=cellinfo2$cell.metadata

for(ct in cats){
 
if(!is.null(allcolors[[ct]])){  
 cellinfo2$cell.annotation[[ct]]= factor(cellinfo2$cell.annotation[[ct]], levels=allcolors[[ct]] %>% names) 
}
}


ph= ComplexHeatmap::pheatmap(mat2,
                             scale="row",
                             cluster_row=FALSE,
                             show_colnames=FALSE,
                             show_rownames=FALSE,
                             cluster_col=FALSE, 
                             breaks=br,
                             col=colls, 
                             annotation_col=cellinfo2$cell.annotation,
                             annotation_row=cellinfo2$marker.annotation,
                             gaps_col= cellinfo2$gaps.cells,
                             gaps_row= cellinfo2$gaps.markers,
                             annotation_colors=colorlist,
                             border_color=NA,
                             right_annotation=ha, 
                             fontsize=5,
                             main=name 
                             #show_heatmap_legend=show.legend
                             #column_gap = unit(.2, "mm")
)


if(return.everything){
  
list(heatmap=ph, cellinfo=cellinfo2, matrix=mat2)
}else{

if(return.cellinfo){
 return(cellinfo2) 
}else{
  
  return(ph)
  }
  
}

}



################################################################################
#refresh cellinfo object
################################################################################


refresh.cellinfo=function(cellinfo, cell.group.label="cell.group", marker.group.label="marker.group"){
  fcat("readjusting")
  
 cellinfo$markers= cellinfo$markerlist %>% Reduce(c, .)
 cellinfo$cells= cellinfo$cell.list %>% Reduce(c, .)
 fcat("making references")
 cell.reference=lapply(1:length(cellinfo$cell.list), function(nn) rep(names(cellinfo$cell.list)[nn], length(cellinfo$cell.list[[nn]])   )) %>% Reduce(c, .)  
marker.reference=lapply(1:length(cellinfo$markerlist), function(nn) rep(names(cellinfo$markerlist)[nn], length(cellinfo$markerlist[[nn]])   )) %>% Reduce(c, .)

##making sure there arent deduped genes

tff=!duplicated(cellinfo$markers)
cellinfo$markers=cellinfo$markers[tff]
marker.reference=marker.reference[tff]


tfc=!duplicated(cellinfo$cells)
cellinfo$cells=cellinfo$cells[tfc]
cell.reference=cell.reference[tfc]
 
## reassemble deduplicated lists
markernames=cellinfo$markerlist %>% names
cellnames=cellinfo$cell.list %>% names

cellinfo$markerlist=lapply(names(cellinfo$markerlist), function(x) cellinfo$markers[marker.reference==x]) %>% givename(., markernames)
cellinfo$cell.list=lapply(names(cellinfo$cell.list), function(x) cellinfo$cells[cell.reference==x]) %>% givename(., cellnames)


 cellinfo$gaps.markers=lapply(1:length(cellinfo$markerlist), function(x) rep(x, length(cellinfo$markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
 cellinfo$gaps.cells= lapply(1:length(cellinfo$cell.list), function(x) rep(x, length(cellinfo$cell.list[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which

fcat("making annotations")
if(!is.null(cellinfo$clustering.variable)){
  cell.group.label=cellinfo$clustering.variable
  
  if(zsoverlap(names(cellinfo$markerlist), names(cellinfo$cell.list))==1 ){
  fcat("Labels of cells and markers are shared.\n Assuming", cellinfo$clustering.variable,"as common variable for markers and cells")
    marker.group.label=cellinfo$clustering.variable
    }
}else{
  fcat("No clustering variable. Checking metadata labels")
if(!is.null(cellinfo$cell.metadata)){
  fcat("Cell metadata found")
  cell.group.label=cellinfo$cell.metadata[1]
}

if(!is.null(cellinfo$marker.metadata)){
  fcat("Marker metadata found")
  marker.group.label=cellinfo$marker.metadata[1]
}
  }

#if(is.null(cellinfo$cell.annotation)){}
previous.cell.annotation=cellinfo$cell.annotation
fresh.cell.annotation=as.data.frame(cell.reference) %>% giverownames(., cellinfo$cells) %>% givecolnames(., 1, cell.group.label)



if(is.null(previous.cell.annotation) || length(intersect(rownames(previous.cell.annotation), cellinfo$cells))==0){
cellinfo$cell.annotation= fresh.cell.annotation
}else{
#bring in the previous cell information dat for the new cells.
othercols=setdiff(colnames(previous.cell.annotation), cell.group.label) %>% unique

added.cell.annotation=previous.cell.annotation %>% getrows(., rownames(fresh.cell.annotation)) %>% select(all_of(othercols))

cellinfo$cell.annotation=cbind(fresh.cell.annotation, added.cell.annotation)
   
}

## annotation of markers
previous.marker.annotation=cellinfo$marker.annotation
fresh.marker.annotation= as.data.frame(marker.reference) %>% giverownames(., cellinfo$markers) %>% givecolnames(., 1, marker.group.label)

if(is.null(previous.marker.annotation)|| length(intersect(rownames(previous.marker.annotation), cellinfo$markers))==0){
cellinfo$marker.annotation= fresh.marker.annotation
}else{
#bring in the previous marker information dat for the new cells, except of the first column
othercols.markers=setdiff(colnames(previous.marker.annotation), marker.group.label) %>% unique

added.marker.annotation=previous.marker.annotation %>% getrows(., rownames(fresh.marker.annotation)) %>% select(all_of(othercols.markers))
cellinfo$marker.annotation=cbind(fresh.marker.annotation, added.marker.annotation)
   
}


cellinfo$cell.metadata= colnames(cellinfo$cell.annotation)
cellinfo$marker.metadata= colnames(cellinfo$marker.annotation)

 cellinfo
 }
  
################################################################################
#update clustering variable
################################################################################
update.group.vars=function(cellinfo, global=NULL, cells.group=NULL, markers.group=NULL){
  if(!is.null(global)){
  cellinfo$clustering.variable=global
  cellinfo$cell.metadata[1]=global
  cellinfo$marker.metadata[1]=global
  subst=colnames(cellinfo$cell.annotation)
   subst[1]=global
   colnames(cellinfo$cell.annotation)=subst
   substm=colnames(cellinfo$marker.annotation)
   substm[1]=global
   colnames(cellinfo$marker.annotation)=substm
  
  }else{

    if(!is.null(cells.group)){
  cellinfo$clustering.variable=cells.group
  cellinfo$cell.metadata[1]=cells.group
  
  subst=colnames(cellinfo$cell.annotation)
   subst[1]=cells.group
   colnames(cellinfo$cell.annotation)=subst
  
    }
    
     if(!is.null(markers.group)){
  
 
  cellinfo$marker.metadata[1]=markers.group

   substm=colnames(cellinfo$marker.annotation)
   substm[1]=markers.group
   colnames(cellinfo$marker.annotation)=substm
  
  }
      
      }
  cellinfo

}

################################################################################
# seriate genes given a cellinfo object
################################################################################

seriategenes=function(so, cellinfo, assay="SCT", seriate.groups=F, link.markers.cells=T){
  
  #mat=tryCatch({so[cellinfo$markers, ][[assay]]@data %>% as.matrix}, function(e) {
   # warning("problem extracting assay data matrix. extracting counts instead")
    #so[st, dscells][["RNA"]]@counts %>% as.matrix
    #})
  
  mat=so[cellinfo$markers, ][[assay]]@data %>% as.matrix
  absent.genes= setdiff(cellinfo$markers, rownames(mat))
  vr= get.variant.rows(mat) 
 
  cellinfo= removemarkers(cellinfo, union(mat[!vr,] %>% rownames, absent.genes))
   mat=mat[vr, ] 
clusmeans=c()
  ct=1
  fcat("Seriating genes")
ordered.gene.list=lapply(1: length(cellinfo$markerlist), function(nn){
  
  st= cellinfo$markerlist[[nn]]
  if(!link.markers.cells){
    dscells=colnames(mat)
  }else{
  dscells=cellinfo$cell.list[[nn]]
  }
  if(length(st)>1){
  ss=seriation::get_order(
    #assumes an ordered correspondence between cell markers and cells in cellinfo
    seriation::seriate(scale((mat[st, dscells]) %>% t) %>%t)
                       , method="PCA")
  clusmeans[ct]<<- mean(colMeans(scale((mat[st, dscells]) %>% t) %>%t)*(sapply(1:length(dscells), function(x) x**2, USE.NAMES=F)))
  ct<<-ct+1
  st[ss %>% unname]}else{
    
    clusmeans[ct]<<- mean(scale((mat[st, dscells]) %>% as.numeric)*(sapply(1:length(dscells), function(x) x**2, USE.NAMES=F)))

    st
    }
  }) %>% givename(., cellinfo$markerlist %>% names)


if(seriate.groups){
names(clusmeans)= 1:nk
neworder=clusmeans %>% sort %>% names
newnames=names(cellinfo$markerlist)[neworder]
cellinfo$markerlist=lapply(neworder, function(tt) ordered.gene.list[[tt]]) %>% givename(., newnames)
}else{

  cellinfo$markerlist=ordered.gene.list 
}

cellinfo.final=refresh.cellinfo(cellinfo)

cellinfo.final
  }



################################################################################
# cellinfo related functions
################################################################################
removenas=function(x) x[!is.na(x)]
adjust.markers= function(cellinfo, n){
  cellinfo$markerlist=lapply(cellinfo$markerlist, function(x) removenas(x[1:n])) %>% givename(., cellinfo$markerlist %>% names) 
cellinfo=refresh.cellinfo(cellinfo)
  cellinfo
}


################################################################################
# adjust cells (get top)
################################################################################

adjust.cells= function(cellinfo, n){
  cellinfo$cell.list=lapply(cellinfo$cell.list, function(x){
    
    if(n<=length(x)){
    removenas(x[(length(x)-n):length(x)])
    }else{
      x
    }
    
    
    
    }) %>% givename(., cellinfo$cell.list %>% names) 
cellinfo=refresh.cellinfo(cellinfo)
  cellinfo
}


################################################################################
# cellinfo collect values for a variable inside the cellinfo
################################################################################

cellinfo.get.values= function(cellinfo, variable, so){
  
metacols= so %>% metadata  %>% colnames
genenames= rownames(so)

if(variable %in% metacols){
  
  for(xx in names(cellinfo$cell.list)){
   
    cellinfo$values[[variable]][[xx]]= (so %>% metadata)[cellinfo$cell.list[[xx]], variable] 
  }
  
}


if(variable %in% genenames){
 
  
  metamat=join_meta_exp(so, genes=variable, assay=DefaultAssay(so))
  for(xx in names(cellinfo$cell.list)){ 
  cellinfo$values[[variable]][[xx]]=  metamat[cellinfo$cell.list[[xx]],variable ]
     
     
  }
  
}
cellinfo
}


################################################################################
# cellinfo single dataset workflow
################################################################################


cellinfo.workflow.singleds=function(so, grouping.variable="seurat_clusters", assay="RNA", markers=NULL, recalculate.markers=F, marker.params=NULL, colorlist=NULL, ncells=300, return.seurat=F){

if(assay=="RNA"){
 fcat("Warning: assay is set as RNA. consider using a SCT assay for better visualisation") 
}
# housekeeping.
  #check that there is an overlap between categories of clusvar and categories of group1 in markers table  
###############################################################################
#Step A. locate or recalculate markers. 
###############################################################################

if(is.null(so$RNA@misc$top_markers) & is.null(markers)){
 fcat("no markers found. Please calculate markers and pass them to argunment markers or to so$RNA@misc$top_markers") 

  }
##############
# Step b.  generate cellinfo with marker signatures. 
#############
if(is.null(colorlist)){
clusterss=NULL  
}else{
 clusterss=names(colorlist[[grouping.variable]])
   
}
  
cellinfo= seriatecells(so, clusvar=grouping.variable, clusters= clusterss, meth="seurat", extended.output=T, deduped=T, ncells=ncells)

outs= cellinfo.heatmap(so, cellinfo,return.everything=T, genes.to.label=adjust.markers(cellinfo, 2)$markers, colorlist=colorlist)
if(return.seurat){
 outs[["seuratobject"]]=so 
}

  outs
  
  }  

################################################################################
# cellinfo- create ranked cell plot... single variable
################################################################################

cellinfo.ordered.cell.plot= function(cellinfo, variable, scale.x=T, clusters=names(cellinfo$cell.list), colorlist=NULL, groupvar=cellinfo$clustering.variable){

  if(is.null(colorlist)){
    
   colorlist=list()
   
   catnames=names(cellinfo$cell.list)
   colorlist[[groupvar]]= randomcolors(length(catnames)) %>% givename(., catnames)
  }

  # structure data in df for ggplot
  

plotdf=lapply(clusters, function(gr){  
numcells= length(cellinfo$cell.list[[gr]])

if(scale.x){
 tot=numcells 
}else{
 tot=1 
}


cellranks=(1:numcells)/numcells
vals=cellinfo$values[[variable]][[gr]]

mt1=cellranks %>% as.data.frame %>% givecolnames(., nms="relative.rank.in.group")
mt1[, variable]=vals
mt2=giverownames(mt1, cellinfo$cell.list[[gr]]) %>% mutate(group=!!gr)

mt2    
}) %>% Reduce(rbind, .)

ggplot(plotdf)+geom_path(aes(x=relative.rank.in.group, y=!!sym(variable), color=group))+theme_classic()+xlab("Cell's rank in group")+NoLegend()+scale_color_manual(values=colorlist[[groupvar]])#+ylab(!!sym(variable))

}




################################################################################
#get markers
################################################################################

get.markers=function(so){
 
  so$RNA@misc$markers 
}

get.top.markers=function(so){
 
  so$RNA@misc$top_markers 
}

get.significant.markers=function(so){
 
  so$RNA@misc$significant_markers
}


################################################################################
# format gene vector for enrichr
################################################################################

format.enrichr=function(x) x %>% paste(., collapse="\n") %>% cat


################################################################################
# function to subset a cellinfo to a few subgroups.
################################################################################
subset.cellinfo=function(cellinfo, labels){
  fcat("restricting to labels only present in list")
  tf=lapply(labels, function(x) x %in% (cellinfo$markerlist %>% names)) %>% Reduce(c, .)
  labels2=labels[tf]
  
  cellinfo$markerlist=cellinfo$markerlist[labels2]
  cellinfo$cell.list=cellinfo$cell.list[labels2]
  cellinfo$input.cluster.ids=labels2
refresh.cellinfo(cellinfo)
  }
  
fillmat=function(genemat, xx){ 
      if(is.null(genemat[[xx]]))
      {genemat[[xx]]=0}
      
      return(genemat)}

add.cell.annotation= function(so, cellinfo, vars, assay="SCT", overwrite=T){
  
  gene.vars=vars[vars %in% rownames(so)]
  non.gene.vars=vars[!(vars %in% rownames(so))]
  meta.vars= non.gene.vars[non.gene.vars %in% (metadata(so) %>% colnames)]
  non.meta.vars= non.gene.vars[!(non.gene.vars %in% (metadata(so) %>% colnames))]
  
  absent.vars=intersect(non.gene.vars, non.meta.vars)
  if(length(absent.vars)>0){
  fcat("Warning: the following variables are absent from the dataset: ", paste(absent.vars, collapse=" "))
    }
  
  clls=cellinfo$cell.annotation %>% rownames
  
  
 new.annotation=metadata(so)[clls,] %>% select(all_of(meta.vars))
  #%>% givecolnames(., nms=meta.vars)
  
 #if overwrite is false, we only incorporate new columns
 if(overwrite==F){
   addition.vars=setdiff(meta.vars, cellinfo$metadata)
 }else{
   #if overwrite is true then we allow replacement of old variables
  addition.vars=meta.vars 
 }
 
 for(mv in addition.vars){
 cellinfo$cell.annotation[[mv]]= new.annotation[[mv]]
 }
 
 
     if(length(gene.vars)>0){
    genemat=join_meta_exp(so, genes=gene.vars, cells=clls)
    for(gn in gene.vars){
    genemat=fillmat(genemat, gn)
    }
    #%>%  givecolnames(., nms=gene.vars)
     cellinfo$cell.annotation=cellinfo$cell.annotation %>% cbind(., genemat[clls,] %>% select(all_of(gene.vars)) )
    }
  
  cellinfo$cell.metadata=colnames(cellinfo$cell.annotation)
  cellinfo
  
  
}

remove.cell.annotation=function(cellinfo, variables){
  
cellinfo$cell.annotation= cellinfo$cell.annotation %>% select(-all_of(variables)) 

cellinfo$cell.metadata=  setdiff(colnames(cellinfo$cell.annotation), variables)
cellinfo 
}

remove.marker.annotation=function(cellinfo, variables){
  
cellinfo$marker.annotation= cellinfo$marker.annotation %>% select(-all_of(variables)) 

cellinfo$marker.metadata=  setdiff(colnames(cellinfo$marker.annotation), variables)
cellinfo 
}

################################################################################
# remove genes from a cell info object
################################################################################
removemarkers=function(cellinfo, markers){
 nms=names(cellinfo$markerlist)
  cellinfo$markerlist=lapply(1:length(cellinfo$markerlist), function(x){
   gns=cellinfo$markerlist[[x]]
   
   gns[!(gns %in% markers)]
   
 }) %>% givename(., nms)
 
  refresh.cellinfo(cellinfo)
}


################################################################################

#szymkiewicz_simpson_coefficient to calculate size of overlaps between vectors
zsoverlap <- function(vector1, vector2) {

  # Calculate the size of the intersection
  intersection_size <- length(intersect(vector1 , vector2))
#fcat(intersection_size)
  # Calculate the minimum size of the two sets
  if(length(vector1)==0 || length(vector1)==0){
    fcat("Warning: some of the sets are zero length")
  }
  min_size <- min(length(vector1), length(vector2))
#fcat(min_size)
  # Calculate the Szymkiewicz–Simpson coefficient
  szymkiewicz_simpson <- intersection_size / min_size

  return(szymkiewicz_simpson)
}


################################################################################
# filter cells based on a prediciton score for a mapping.
################################################################################

  filterps= function(so, varbl="seurat_clusters", th) so[, so %>% metadata %>% filter(!!sym(paste0(mapfun(varbl), "_predicted.id.score"))>=!!th ) %>% rownames]



################################################################################
# selectively remove heatmap legends from a ComlpexHeatmap object
################################################################################


removelegends= function(hm, legends, where="top"){

  annoslot=paste0(where, "_annotation")
  #for(ann in c("top_annotation", "right_annotation", "left_annotation", "bottom_annotation")){
  for(leg in legends){
  slot(hm, annoslot)@anno_list[[leg]]@show_legend=F
  }
   
  hm 
  #}
}


################################################################################
#   classify cells into tf, ligands and receptors
################################################################################


#simpleCache(paste_("genes", "db", "human_tfs"), {
 # m <- fread("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")
  #> m[,.N,by=`Is TF?`]
  # Is TF?    N
  # 1:    Yes 1639
  # 2:     No 1126
  #m <- m[`Is TF?`=="Yes", `HGNC symbol`] #paste0(`Ensembl ID`,"-",`HGNC symbol`)]
  #m
#}, assignToVar="tfList",  reload=T)
#istf= function(x) x %in% tfList 
#lrtable=readRDS(file.path(params$resourcepath, "cell_cell_interactions/CellTalkDB/human_lr_pair.rds"))

#isligand=function(x) x %in%  (lrtable %>% pull(ligand_gene_symbol))
#isreceptor=function(x) x %in%  (lrtable %>% pull(receptor_gene_symbol))
#islrpair=function(x) x %in%  (lrtable %>% pull(lr_pair))


colorgene=function(x) {
  
  if(istf(x)){return("black")}else{
    if(isligand(x)){return("blue")}else{
      if(isreceptor(x) || (x %in% c("CRABP1", "RARRES1", "RARRES2", "CRABP2"))){return("red")}else{return("grey")}
    }
  }}

colorgene2=function(x) {
  
  if(istf(x)){return("black")}else{return("grey")}}
    
    



#alosteric protein database
#https://mdl.shsmu.edu.cn/ASD2023Common/static_file/archive_2023/ASD_Release_202309_AS.tar.gz

################################################################################
#
################################################################################

divide_legend_columns <- function(plot, num_columns) {
  # Get the current legend
  current_legend <- plot + theme(legend.position = "bottom")  # Adjust the position if needed
  
  # Modify the legend to have x columns
  modified_legend <- current_legend +
    guides(color = guide_legend(nrow = 1, ncol = num_columns))

  return(modified_legend)
}



################################################################################
# No Legend title
################################################################################
NoLegendTitle <- function(plot) {
  # Modify the legend to remove the title
   guides(guide_legend(title = NULL))
  
}
################################################################################
cat("checkpointend-2...\n")
remove_x_axis <- function() {
  # Modify the theme to remove x-axis
   theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                axis.ticks.x = element_blank())
}

################################################################################
# scale text size
################################################################################

scale.text=function(scaling_factor){
    theme(
      text = element_text(size =  scaling_factor),
      axis.title = element_text(size =  scaling_factor),
      axis.text = element_text(size = scaling_factor),
      axis.title.x = element_text(size = scaling_factor),
      axis.text.x = element_text(size =  scaling_factor),
      axis.title.y = element_text(size =  scaling_factor),
      axis.text.y = element_text(size =  scaling_factor),
      legend.text = element_text(size =  scaling_factor),
      legend.title = element_text(size =  scaling_factor),
      plot.title = element_text(size =  scaling_factor),
      plot.subtitle = element_text(size =  scaling_factor),
      plot.caption = element_text(size =  scaling_factor)
    )
}



scale.text.general=function(plt, sca){
 plt$theme$text$size= plt$theme$text$size*sca
  plt
}


################################################################################
# masking types in kameneva et al cell types
################################################################################
masktype=Vectorize(function(x) if(x =="SCP" || x=="sympathoblasts" || x=="mesenchyme"){return(x)}else{return("other")}, USE.NAMES=F)
 
masktype3=Vectorize(function(x) if(x =="SCP" || x=="sympathoblasts" || x=="mesenchyme" || x=="chromaffin"){return(x)}else{return("other")}, USE.NAMES=F)
 
################################################################################
#replace markers within a seurat object
################################################################################

replacemarkers=function(so, newmarkers){
 if(!is.null(newmarkers)){
  so$RNA@misc$top_markers=newmarkers 
  return(so)
  }else{
 return(so)
  }}
cat("checkpointend-1...\n")
################################################################################
# whether a gene is in chr 17 or chr1
################################################################################
  
  if(!exists("gtftable")){ 
gtfpath=params$gtfpath

gtftable=rtracklayer::import(gtfpath)
}
genepositions=gtftable %>% as.data.frame %>% select(seqnames, start, gene_name, type)%>% filter(type=="gene")


get.seqnames=Vectorize(function(x) (genepositions %>% filter(gene_name==x))$seqnames %>% as.character, USE.NAMES=F)

is.in.17=Vectorize(function(x){ if(!is.na(x) && grepl("^chr17$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.1=Vectorize(function(x){ if(!is.na(x) && grepl("^chr1$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.2=Vectorize(function(x){ if(!is.na(x) && grepl("^chr2$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.3=Vectorize(function(x){ if(!is.na(x) && grepl("^chr3$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.4=Vectorize(function(x){ if(!is.na(x) && grepl("^chr4$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.5=Vectorize(function(x){ if(!is.na(x) && grepl("^chr5$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.6=Vectorize(function(x){ if(!is.na(x) && grepl("^chr6$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.7=Vectorize(function(x){ if(!is.na(x) && grepl("^chr7$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.8=Vectorize(function(x){ if(!is.na(x) && grepl("^chr8$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.9=Vectorize(function(x){ if(!is.na(x) && grepl("^chr9$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.10=Vectorize(function(x){ if(!is.na(x) && grepl("^chr10$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.11=Vectorize(function(x){ if(!is.na(x) && grepl("^chr11$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.12=Vectorize(function(x){ if(!is.na(x) && grepl("^chr12$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.13=Vectorize(function(x){ if(!is.na(x) && grepl("^chr13$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.14=Vectorize(function(x){ if(!is.na(x) && grepl("^chr14$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.15=Vectorize(function(x){ if(!is.na(x) && grepl("^chr15$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.16=Vectorize(function(x){ if(!is.na(x) && grepl("^chr16$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.17=Vectorize(function(x){ if(!is.na(x) && grepl("^chr17$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.18=Vectorize(function(x){ if(!is.na(x) && grepl("^chr18$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.19=Vectorize(function(x){ if(!is.na(x) && grepl("^chr19$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.20=Vectorize(function(x){ if(!is.na(x) && grepl("^chr20$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.21=Vectorize(function(x){ if(!is.na(x) && grepl("^chr21$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)
is.in.22=Vectorize(function(x){ if(!is.na(x) && grepl("^chr22$", get.seqnames(x))){return(TRUE)}else{return(FALSE)}}, USE.NAMES=F)



################################################################################
#
################################################################################

rmquery=Vectorize(function(x){if(grepl("query",x)){
strsplit(x, split="query_")[[1]][2]
  }else{x} }, USE.NAMES=F)



rmquery.matrix= function(x){
  x %>% names2col(., "cellid") %>% mutate(cellid=rmquery(cellid)) %>% col2names(., "cellid")
}


cat("checkpointend...\n")

################################################################################
#  remove duplicates from a vector
################################################################################

deduplicate= function(x) x[!duplicated(x)]


################################################################################
# quick print a matrix's contents
################################################################################

quickprint=function(mat, wide=F, limit=5) if(wide){ mat[1:limit, ]}else{mat[1:limit, 1:limit]}


################################################################################
# fit hill function 
################################################################################

hillfunction <- function(x, k, n, beta) beta*((x**n)/((k**n)+(x**n)))

#model1 <- nls(fluorI ~ eDecay(t,myA,myT), data=ExpData, start=list(myA=10,myT=5))



################################################################################
# gglabel (auxiliary side bar plot to make a legend
################################################################################

manual.legend=function(colorlist, pname, size=7, trailing.spaces=15, bar.width=.6){


  p=colorlist[[pname]] %>% rev
  
  pl=ggplot(1:length(p) %>% as.data.frame %>% mutate(h=.1, colors=p, name=names(p)), aes(x=factor(.), y=factor(h)))+
    geom_col(aes(fill=factor(.)))+
    scale_fill_manual(values=namenums(p))+
    coord_cartesian( ylim=c(0,10))+
    theme_classic()+
    #geom_text(y=0.5, color="black", angle=90)+
    NoAxes()+
    NoLegend()+
    ggtitle(paste0(paste(rep(" ", trailing.spaces), collapse=""), pname))
  
  return(pl+geom_text(inherit.aes=T,aes(label=name), y=1.05, size=size, color="black", hjust=0)+coord_flip(ylim=c(0,1/bar.width)))
}






