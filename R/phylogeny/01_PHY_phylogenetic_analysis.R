################################################################################

reinstall = F

if(reinstall) {

  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  # Bioconductor version 3.16 (BiocManager 1.30.17), R 4.2.0 (2022-04-22)

  library(devtools)

  install_version("reshape2", version = "1.4.4", repos = "http://cran.us.r-project.org", upgrade = FALSE)
  install_version("phangorn", version = "2.11.1", repos = "http://cran.us.r-project.org", upgrade = FALSE)
  install_version("ggplot2", version = "3.6.2", repos = "http://cran.us.r-project.org", upgrade = FALSE)

  BiocManager::install("ComplexHeatmap", update = FALSE) # 2.14.0
  BiocManager::install("ggtree", update = FALSE) # 3.6.2

  install_version("circlize", version = "0.4.15", repos = "http://cran.us.r-project.org", upgrade = FALSE)
  install_version("adephylo", version = "1.1-16", repos = "http://cran.us.r-project.org", upgrade = FALSE)
  install_version("cowplot", version = "1.1.1", repos = "http://cran.us.r-project.org", upgrade = FALSE)
  
}

################################################################################

library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(cowplot)
library(adephylo)
library(GenomicRanges)

################################################################################

# Select parameters
# Minimum variant allele frequency used for analysis
vaf = 0.2

# Output paths
output_path = "" # Where would you like the analysis to go?
bed_file    = "metadata/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed" # Path to targeted bed file
input_tsv   = "" # ncnb_smlVars_all_tn_pass_tiers.tsv, tsv of called variants

################################################################################

sample_names = c("") # Which samples do we want to run this on?

################################################################################

# Create output dir
out_dir = paste0("minvaf_",vaf)
dir.create(paste0(output_path,"/",out_dir))

# Read in the bed file of the exome panel
twist_bed = read.table(bed_file, 
                       stringsAsFactors = F)

# Read in mutation output
muts = read.table(input_tsv,
                  sep = "\t", header = T)

# Make GRanges objects
twist_gr = GRanges(seqnames = twist_bed$V1, ranges = IRanges(start = twist_bed$V2, 
                                                             end = twist_bed$V3))
muts_gr = GRanges(seqnames = muts$seqnames, IRanges(start = muts$start, end = muts$end))

# Find overlaps and subset table
muts = muts[queryHits(findOverlaps(muts_gr, twist_gr)),]

# Use Mutect2 calls
muts = muts[muts$tool == "Mutect2",]

# Get just unique stuff
muts = muts[!duplicated(paste0(muts$varkey,"_",muts$sample)),]

# Simple filtering
muts = muts[muts$FILTER=="PASS" & muts$VARIANT_CLASS=="SNV",]

# Mutation matrix
mut_mat = dcast(muts, varkey ~ sample, value.var = "T_VAF")

# Make mut names the rows
rownames(mut_mat) = mut_mat$varkey
mut_mat$varkey = NULL

# Replace NAs with zeros
mut_mat[is.na(mut_mat)] = 0

# Remove filter based on sample names
mut_mat = mut_mat[,colnames(mut_mat) %in% sample_names]

# Make binary
mut_mat_bin = as.matrix(mut_mat >= vaf) + 0

# Add false root
mut_mat_bin = cbind(mut_mat_bin, 0)
colnames(mut_mat_bin)[ncol(mut_mat_bin)] = "Parental"

# Add another false root
mut_mat_bin = cbind(mut_mat_bin, 0)
colnames(mut_mat_bin)[ncol(mut_mat_bin)] = "For_rooting"

# Make phyData
phyHs = as.phyDat(t(mut_mat_bin), type = "USER", levels = c("0", "1"))

# Make the tree we'll go with
terry = pratchet(phyHs, maxit = 1000)

# Root it
terry = root(terry, outgroup = "Parental", resolve.root = TRUE)

# Get branch lengths
terry = acctran(terry, phyHs)

# Now remove the tip we added for the purpose of rooting
terry = drop.tip(terry, "For_rooting")

# Import gene
gene = expression(paste(italic("BCOR"), " L1673F")) # Requires manual input

# Replace underscores
terry$tip.label = gsub("_", " ", terry$tip.label)

# Define colours
fills = c("1q" = "grey", "wtMYCN" = "#000000", "wtMYCN D19" = "#000000", 
          "wtMYCN DOX D19" = "#000000", "17q1qMYCN" = "orange", 
          "17q1qMYCN DOX D19" = "magenta",
          "17q1qMYCN D19" = "orange", "17qMYCN" = "#069af3", 
          "17qMYCN D19" = "#069af3", "17q" = "#00AD67", 
          "17q1qMYCN DOX" = "magenta", "17qMYCN DOX" = "#00ffff", 
          "17qMYCN DOX D19" = "#00ffff", "Parental" = "#777777")

# Premake tip preferences
fills    = fills[terry$tip.label]
outlines = fills
outlines[c(grep("D19", terry$tip.label))] = "black"
shapes   = rep(21, times = length(fills))
shapes[grep("DOX", terry$tip.label)] = 22

# This creates Figure 5E
p = ggtree(terry, size = 1) + 
  geom_tiplab(hjust = -c(1 / nchar(terry$tip.label)), size = 6) +
  geom_tippoint(fill = fills, 
                col = outlines,
                shape = shapes,
                size = 6, stroke = 1.5) + 
  ggplot2::xlim(0, max(ape::node.depth.edgelength(terry)) * 1.5) + 
  annotate("text", x = 0.2 * (max(ape::node.depth.edgelength(terry)) * 1.5), 
           y = ceiling((8 / 11) * (terry$Nnode + 1)) - 0.25, label = gene) + # this is manually added
  geom_treescale(x = 0.1 * (max(ape::node.depth.edgelength(terry)) * 1.5), 
                 y = terry$Nnode + 1.5, linesize = 1.25, offset = 0.1, width = 10) +
  theme(title=element_text(size=20, colour = "gray20"), 
        plot.title = element_text(hjust = 0.5)) 

ggsave(p, filename = paste0(output_path,"/",out_dir,"/NCNB_ggtree.pdf"), 
       width = 8, height = 5)

# Get the tip distances
tip_dists = distTips(terry)

# Make it a matrix
tip_dists = as.matrix(tip_dists)

# Make it a barplot
bar_df = melt(sort(tip_dists["Parental",]))

# Plot data wrangling
bar_df$Sample = rownames(bar_df)
bar_df$Sample = factor(bar_df$Sample, levels = bar_df$Sample[order(bar_df$value)])
colnames(bar_df)[1] = "Distance"

# Do the stats test
wilcox.test(bar_df[c("wtMYCN", "17qMYCN", "17q1qMYCN"),], 
            bar_df[c("wtMYCN DOX D19", "17qMYCN DOX D19", "17q1qMYCN DOX D19"),], 
            paired = TRUE, alternative = "two.sided")

# This creates Supplementary Figure 9G
p = ggplot(bar_df, aes(x = Sample, y = Distance, fill = Sample)) +
  ylab("Mutation Distance form Parental") + 
  geom_bar(stat="identity", col = "black") + 
  scale_fill_manual(values = fills) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), legend.position = "none")

ggsave(p, filename = paste0(output_path,"/",out_dir,"/NCNB_Parental_dist.pdf"), 
       width = 5, height = 5)

################################################################################
