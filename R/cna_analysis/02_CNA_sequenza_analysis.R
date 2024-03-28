reinstall = F

if(reinstall) {

  if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  # Bioconductor version 3.16 (BiocManager 1.30.17), R 4.2.0 (2022-04-22)

  BiocManager::install("copynumber", update = FALSE) # 1.38.0
  BiocManager::install("ComplexHeatmap", update = FALSE) # 2.14.0
    
  library(devtools)

  install_version("circlize", version = "0.4.15", repos = "http://cran.us.r-project.org", upgrade = FALSE)
  
  install_version("sequenza", version = "3.0.0", repos = "http://cran.us.r-project.org", upgrade = FALSE)

}

################################################################################

options(scipen = 12)

library("sequenza")
library("ComplexHeatmap")
library("circlize")
library("GenomicRanges")

################################################################################

# Paths
input_path  = "" # Where did you output the sequenza results?
output_path = "" # Where do you want the results from this analysis?
sample_ids  = c("") # Which samples do you want to analyse?

################################################################################

# Perform sequenza workflow for each sample
for(sample in sample_ids) {

  cna_data <- sequenza.extract(paste0(input_path,"/NCNB_",sample,"_parental_normal_binned.seqz.gz"), 
                               chromosome.list = paste0("chr",c(1:22,"X")))
  
  CP.example <- sequenza.fit(cna_data, female = TRUE, XY = c(X = "chrX", Y = "chrY"),
                             cellularity = seq(1,1.0002,0.0001))
  
  sequenza.results(sequenza.extract = cna_data, cp.table = CP.example,
                   sample.id = paste0("NB_",sample), out.dir=paste0(output_path,"/NB_",sample),
                   female = TRUE, XY = c(X = "chrX", Y = "chrY"))
  
}

################################################################################
################ Code for the copy number matrix heatmap #######################
################################################################################

# Load up an extract object
load(paste0(output_path,"/NB_",sample_ids[1],"/NB_",sample_ids[1],"_sequenza_extract.RData"))

# Get around the varying object name
NB_sequenza_extract = get(paste0(sample_ids[1],"_sequenza_extract"))

# Get the bins
bins = lapply(1:length(NB_sequenza_extract$BAF), function(i) {

  res = NB_sequenza_extract$BAF[[i]]
  res$chr = names(NB_sequenza_extract$BAF)[i]

  # Bins are overlapping by half, we'll remove this
  len = nrow(res)
  res = res[-len,]
  res = res[seq(1,len-1,by = 2),]

  # Make it so end isn't same as next start
  res$end = res$end - 1
  
  return(res)
})

# Format
bins = do.call(rbind, bins)

# Important columns and in order
bins = bins[,c("chr", "start", "end")]

# Bins reference
bins_gr = GRanges(seqnames = bins$chr, ranges = IRanges(start = bins$start, end = bins$end))

# Add it to a data.frame
cn_data = bins

# Run through each sample
for(sample in sample_ids) {

  # Read it in!
  segs = read.table(paste0(output_path,"/NB_",sample,"/NB_",sample,"_segments.txt"), 
                    header = T, sep = "\t", stringsAsFactors = F)
  
  # Segments as a GRanges
  segs_gr = GRanges(seqnames = segs$chromosome, ranges = IRanges(start = segs$start.pos, segs$end.pos))
  
  # Get overlaps
  overlaps = findOverlaps(query = segs_gr, subject = bins_gr)
  queryhit = queryHits(overlaps)
  subjhit  = subjectHits(overlaps)
  
  # Entries
  cn = rep(NA, times = length(bins_gr))
  
  # for loop to add in CN states
  for(i in 1:length(overlaps)) {
    
    same_subjhit = which(subjhit[i] == subjhit)

    if(length(same_subjhit) == 1) {

      # if there is only one overlap the replacement is simple
      cn[subjhit[i]] = segs$CNt[queryhit[i]]

    }

    if(length(same_subjhit) > 1) {

      # if there are multiple overlaps find the widths of the intersects 
      bin_ovlps = data.frame(width = width(intersect(bins_gr[subjhit[i]], segs_gr[queryhit[same_subjhit]])), 
                             cn = segs$CNt[queryhit[same_subjhit]])

      # Order by copy number to help calculate the median
      bin_ovlps = bin_ovlps[order(bin_ovlps$cn),]

      # Calculate the cumsum in CN order
      bin_ovlps$cumsum = cumsum(bin_ovlps$width)

      # The median is the first CN that contains the middle value
      bin_ovlps_med_cn = bin_ovlps$cn[which(bin_ovlps$cumsum > (sum(bin_ovlps$width) / 2))[1]]

      # This value now represents the bin
      cn[subjhit[i]] = bin_ovlps_med_cn

    }
    
  }
  
  # Add the column
  cn_data = cbind(cn_data, cn)
  
  # Add the column name
  colnames(cn_data)[ncol(cn_data)] = sample

}

# Remove NAs
cn_data = na.omit(cn_data)
rownames(cn_data) = NULL

# Add levels
cn_data$chr = factor(cn_data$chr, levels = paste0("chr",c(1:22,"X")))

# Chromosome splits
chr.ends = cumsum(table(cn_data$chr))

# Label positions
chr.labs = round(c(0, chr.ends[-length(chr.ends)]) + 
                   rle(as.character(cn_data$chr))$lengths / 2)

# Chromosome selection
chr.sel = c(1:17,seq(18,22,by=2),23)

# Chromosome annotation
bot.column.anno = columnAnnotation(link = anno_mark(at = chr.labs[chr.sel], 
                                                    labels = c(1:22,"X")[chr.sel], 
                                                    side = "bottom",
                                                    link_width = unit(0.25, "cm")))

# Input data
plot_mat = t(cn_data[,4:ncol(cn_data)])

# Change row names
rownames(plot_mat) = gsub("NB_", "", rownames(plot_mat))

# Replace
rownames(plot_mat) = gsub("_", " ", rownames(plot_mat))

# This creates Supplementary Figure 9F
pdf(paste0(output_path,"/NCNB_copy_number_heatmap_sequenza.pdf"), width = 7, height = 7/3)

# Make a heatmap
Heatmap(plot_mat, cluster_rows = T, cluster_columns = F, name = "CN",
        show_row_names = T,
        bottom_annotation = bot.column.anno,
        col = colorRamp2(c(0, 2, 4, 8), c("blue", "white", "red", "darkred")))

#Add lines
for(boundary in c(0,1)) {
  
  #Add the lines
  decorate_heatmap_body("CN", {
    grid.lines(c(0, 1), c(1 - boundary, 1 - boundary), gp = gpar(lwd = 0.5))
  })
  
}

#Add lines
for(boundary in c(0,1)) {
  
  #Add the lines
  decorate_heatmap_body("CN", {
    grid.lines(c(boundary, boundary), c(0, 1), gp = gpar(lwd = 1))
  })
  
}

#Add lines
for(boundary in chr.ends / nrow(cn_data)) {
  
  #Add the lines
  decorate_heatmap_body("CN", {
    grid.lines(c(boundary, boundary), c(0, 1), gp = gpar(lty = "dotted", lwd = 1))
  })
  
}

dev.off()

################################################################################
######## This code generates the depth ratio data for profile plotting #########
################################################################################

depth_ratios = lapply(paste0("NB_",sample_ids), function(sample) {
  
  # Load sample results
  load(paste0(output_path,"/",sample,"/",sample,"_sequenza_extract.RData"))
  
  # Get around the varying object name
  obj = get(paste0(sample,"_sequenza_extract"))
  
  # Get the bins
  drs = lapply(1:length(obj$ratio), function(i) {
    res = obj$ratio[[i]]
    res$chr = names(obj$ratio)[i]
    return(res)
  })
  
  # Collect the data per chromosome
  drs = do.call(rbind, drs)
  
  # Simplify the data
  drs = drs[,c("chr", "start", "end", "mean")]
  
  # Rename the column as the sample
  colnames(drs)[4] = sample
  
  # Match the heatmap binning style
  drs$end = drs$end - 1
  
  # Take only those that overlap with the heatmap
  drs = drs[paste0(drs$chr,"_",drs$start,"_",drs$end) %in% paste0(cn_data$chr,"_",cn_data$start,"_",cn_data$end),]
  
  return(drs)
})

# For neatness
names(depth_ratios) = paste0("NB_",sample_ids)

# Get backbone of bins
dr_bb = depth_ratios[[1]][,1:3]

# Combine it with the data from each sample
depth_ratios_mat = cbind(dr_bb, do.call(cbind, lapply(depth_ratios, function(i) i[,4])))

# Carry over the sample names
colnames(depth_ratios_mat)[4:ncol(depth_ratios_mat)] = sample_ids

################################################################################

# This output is used for Supplementary Figure 4b
write.csv(depth_ratios_mat, file = paste0(output_path,"/NCNB_depth_ratio_matrix.csv"), 
          quote = F, row.names = F)

write.csv(cn_data, file = paste0(output_path,"/NCNB_copy_number_matrix.csv"),
          quote = F, row.names = F)

################################################################################
