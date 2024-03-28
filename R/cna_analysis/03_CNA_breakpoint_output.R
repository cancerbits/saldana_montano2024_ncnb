################################################################################
######## To be run after sequenza_analysis.R is done to extract info ###########
################################################################################

# Paths
input_path = "" # Where did you perform the Sequenza analysis?
sample_ids = c("") # Which samples do you want to perform this on?

################################################################################
######## Subset just the CNAs affecting the key chromosomes in our study #######
################################################################################

# Chromosomes of interest
coordinates = do.call(rbind, lapply(sample_ids, function(sample) {
  
  # Read it in!
  segs = read.table(paste0(input_path,"/NB_",sample,"/NB_",sample,"_segments.txt"), 
                    header = T, sep = "\t", stringsAsFactors = F)
  
  segs$sample = sample
  
  res = segs[segs$chromosome %in% c("chr1", "chr17") & segs$CNt != 2,]
  
  return(res)
  
}))

write.table(coordinates, file = paste0(input_path,"/NCNB_1q_17q_aberrations_coordinates.csv"), 
            sep = ",", quote = F, row.names = F)

################################################################################
######## Code to determine if breakpoints lie in the telomeric regions #########
################################################################################

# Get sizes
sizes = read.table("https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes", 
                   stringsAsFactors = F)

# Get appropriate lengths
chr1_len  = sizes[sizes$V1=="chr1",2]
chr17_len = sizes[sizes$V1=="chr17",2]

# Add size for info
coordinates$size = coordinates$end.pos - coordinates$start.pos

# Make anything close to the telomeres (close = 1Mbp) equal to "tel"
coordinates$start_tel_p = coordinates$start.pos < 10^6
coordinates$end_tel_p   = coordinates$end.pos < 10^6
coordinates$start_tel_q = ifelse(coordinates$chromosome=="chr1", 
                                 coordinates$start.pos > (chr1_len - 10^6), 
                                 ifelse(coordinates$chromosome=="chr17", 
                                        coordinates$start.pos > (chr17_len - 10^6), NA))
coordinates$end_tel_q = ifelse(coordinates$chromosome=="chr1", 
                               coordinates$end.pos > (chr1_len - 10^6), 
                               ifelse(coordinates$chromosome=="chr17", 
                                      coordinates$end.pos > (chr17_len - 10^6), NA))

# Are both breakpoint near the telomere?
coordinates$total_tel_seg = (coordinates$start_tel_p & coordinates$end_tel_p) | (coordinates$start_tel_q & coordinates$end_tel_q)

# Coordinates for hg38 to hg19 conversion
coordinates$convert.start = ifelse(coordinates$start_tel_p, coordinates$end.pos, coordinates$start.pos)
coordinates$convert.end   = ifelse(coordinates$end_tel_q, coordinates$start.pos, coordinates$end.pos)

# What size is the chromosome in total?
coordinates$chr_total_size = ifelse(coordinates$chromosome=="chr1", chr1_len, ifelse(coordinates$chromosome=="chr17", chr17_len, NA))

# Simplify the output
coordinates = coordinates[,c("chromosome", "start.pos", "end.pos", "size", "CNt", 
                             "N.BAF", "sample", "start_tel_p", "end_tel_p", 
                             "start_tel_q", "end_tel_q", "total_tel_seg", 
                             "convert.start", "convert.end", "chr_total_size")]

# This output is used for Supplementary Figure 4C
write.table(coordinates, file = paste0(input_path,"/NCNB_1q_17q_aberrations_coordinates_simple.csv"), 
            sep = ",", quote = F, row.names = F)

################################################################################
################## Simply output all aberrant segments #########################
################################################################################

# All aberrant segments
coordinates = do.call(rbind, lapply(sample_ids, function(sample) {
  
  # Read it in!
  segs = read.table(paste0(input_path,"/NB_",sample,"/NB_",sample,"_segments.txt"), 
                    header = T, sep = "\t", stringsAsFactors = F)
  
  segs$sample = sample
  
  res = segs[segs$CNt != 2,]
  
  return(res)
  
}))

write.table(coordinates, file = input_path,"/NCNB_all_aberrations_coordinates.csv", 
            sep = ",", quote = F, row.names = F)
