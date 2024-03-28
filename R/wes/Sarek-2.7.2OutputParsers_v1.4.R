# to be run in ccribioinf/dockrstudio:4.0.3-v1
# dockerhub manifest digest: sha256:055a81938337a3fc2a15039476094058f5691e02b529578041af508352264323

unzip <- function(x) {
  if(endsWith(x, ".gz")){
    zz <- gzfile(x,'rt')
    zz
  }else{
    x
  }
}

vcf_col_parse <- readr::cols(
  # "1" annotation can break the automatic parsing of readr
  `#CHROM` = readr::col_character(), 
  # POS is always a double
  POS = readr::col_double(),
  # col_guess() (default) assumes the type by considering the first 1000 entries
  # however if it guesses col_numeric() then ","s get silently removed. Better go with col_character()
  .default = readr::col_character()
  )

# from https://github.com/biobenkj/StrelkaParser
row_tibble <- function(x, col_names) {
  tibble::as_tibble(rbind(setNames(x, col_names)))
}

parse_format <- function(format, sample) {
  purrr::map2(
    strsplit(sample, ":", fixed = TRUE),
    strsplit(format, ":", fixed = TRUE),
    row_tibble)
}

parse_info <- function(info) {
  strsplit(info, ";", fixed = TRUE) %>%
    purrr::map(~row_tibble(sub("^.*=(.*)", "\\1", .x), sub("^(.*)=.*", "\\1", .x)))
}

split_tiers <- function(x) {
  name <- names(x)
  x$`1` <- x[[1]] %>%
    strsplit(",", fixed = TRUE) %>%
    purrr::map(~as.integer(.x[1])) %>%
    purrr::flatten_int()
  x$`2` <- x[[1]] %>%
    strsplit(",", fixed = TRUE) %>%
    purrr::map(~as.integer(.x[2])) %>%
    purrr::flatten_int()
  names(x)[2:3] <- paste0(name, "_", names(x)[2:3])
  x[2:3]
}

parse_strelka_snv_con <- function(x) {
  print(paste0("Processing ", basename(x)))
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       info = parse_info(INFO),
                       normal = parse_format(FORMAT, NORMAL),
                       tumor = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, N = normal, T = tumor, .sep = "_")
  vcf <- purrr::lmap_at(vcf, dplyr::matches("[TN]_[ACTG]U", vars = names(vcf)), split_tiers)
  vcf <- readr::type_convert(vcf)
  vcf <- tidyr::gather(vcf, allele, count, dplyr::matches("[TN]_[ACTG]U_[12]"))
  vcf <- tidyr::separate(vcf, allele, c("sample", "base", "tier"), sep = "[_U]+", remove = FALSE)
  vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
  vcf <- dplyr::mutate(vcf, T_REF_COUNT = count[REF == base & sample == "T" & tier == TQSS_NT],
                       T_ALT_COUNT = count[ALT == base & sample == "T" & tier == TQSS_NT],
                       N_REF_COUNT = count[REF == base & sample == "N" & tier == TQSS_NT],
                       N_ALT_COUNT = count[ALT == base & sample == "N" & tier == TQSS_NT],
                       T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT),
                       N_VAF = N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT))
  vcf <- dplyr::select(vcf, -c(normal, tumor, sample, base, tier))
  vcf <- tidyr::spread(vcf, allele, count)  # Takes up most of the time
  vcf
}

# only select tier1 values, motivated by
# https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md
parse_strelka_snv_t1 <- function(x) {
  print(paste0("Processing ", basename(x)))
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       info = parse_info(INFO),
                       normal = parse_format(FORMAT, NORMAL),
                       tumor = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, N = normal, T = tumor, .sep = "_")
  vcf <- purrr::lmap_at(vcf, dplyr::matches("[TN]_[ACTG]U", vars = names(vcf)), split_tiers)
  vcf <- readr::type_convert(vcf)
  vcf <- tidyr::gather(vcf, allele, count, dplyr::matches("[TN]_[ACTG]U_[12]"))
  vcf <- tidyr::separate(vcf, allele, c("sample", "base", "tier"), sep = "[_U]+", remove = FALSE)
  vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
  vcf <- dplyr::mutate(vcf, 
                       T_REF_COUNT = count[REF == base & sample == "T" & tier == 1],
                       T_ALT_COUNT = count[ALT == base & sample == "T" & tier == 1],
                       N_REF_COUNT = count[REF == base & sample == "N" & tier == 1],
                       N_ALT_COUNT = count[ALT == base & sample == "N" & tier == 1],
                       T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT),
                       N_VAF = N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT))
  vcf <- dplyr::select(vcf, -c(normal, tumor, sample, base, tier))
  vcf <- tidyr::spread(vcf, allele, count)  # Takes up most of the time
  vcf
}
# there might be a nice function{x, arg} solution to set if VAF calculation is conditional or tier1.

parse_strelka_indel_con <- function(x) {
  print(paste0("Processing ", basename(x)))
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       info = parse_info(INFO),
                       normal = parse_format(FORMAT, NORMAL),
                       tumor = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, N = normal, T = tumor, .sep = "_")
  vcf <- purrr::lmap_at(vcf, dplyr::matches("[TN]_T[AI]R", vars = names(vcf)), split_tiers)
  vcf <- readr::type_convert(vcf)
  vcf <- tidyr::gather(vcf, allele, count, dplyr::matches("[TN]_T[AI]R_[12]"))
  vcf <- tidyr::separate(vcf, allele, c("sample", "indel", "tier"), sep = "[_]+", remove = FALSE)
  vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
  vcf <- dplyr::mutate(vcf, 
                       T_REF_COUNT = count[indel == "TAR" & sample == "T" & tier == TQSI_NT],
                       T_ALT_COUNT = count[indel == "TIR" & sample == "T" & tier == TQSI_NT],
                       N_REF_COUNT = count[indel == "TAR" & sample == "N" & tier == TQSI_NT],
                       N_ALT_COUNT = count[indel == "TIR" & sample == "N" & tier == TQSI_NT],
                       T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT),
                       N_VAF = N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT))
  vcf <- dplyr::select(vcf, -c(normal, tumor, sample, indel, tier))
  vcf <- tidyr::spread(vcf, allele, count)  # Takes up most of the time
  vcf
}

# only select tier1 values, motivated by
# https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md
parse_strelka_indel_t1 <- function(x) {
  print(paste0("Processing ", basename(x)))
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       info = parse_info(INFO),
                       normal = parse_format(FORMAT, NORMAL),
                       tumor = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, N = normal, T = tumor, .sep = "_")
  vcf <- purrr::lmap_at(vcf, dplyr::matches("[TN]_T[AI]R", vars = names(vcf)), split_tiers)
  vcf <- readr::type_convert(vcf)
  vcf <- tidyr::gather(vcf, allele, count, dplyr::matches("[TN]_T[AI]R_[12]"))
  vcf <- tidyr::separate(vcf, allele, c("sample", "indel", "tier"), sep = "[_]+", remove = FALSE)
  vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
  vcf <- dplyr::mutate(vcf, 
                       T_REF_COUNT = count[indel == "TAR" & sample == "T" & tier == 1],
                       T_ALT_COUNT = count[indel == "TIR" & sample == "T" & tier == 1],
                       N_REF_COUNT = count[indel == "TAR" & sample == "N" & tier == 1],
                       N_ALT_COUNT = count[indel == "TIR" & sample == "N" & tier == 1],
                       T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT),
                       N_VAF = N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT))
  vcf <- dplyr::select(vcf, -c(normal, tumor, sample, indel, tier))
  vcf <- tidyr::spread(vcf, allele, count)  # Takes up most of the time
  vcf
}
# there might be a nice function{x, arg} solution to set if VAF calculation is conditional or tier1.

parse_paired_mutect2 <- function(x) {
  print(paste0("Processing ", basename(x)))
  meta <- file(x, "r")
  
  normal_sample_id <- NULL
  tumor_sample_id <- NULL
  
  # Loop through the lines in the file
  while (is.null(normal_sample_id) || is.null(tumor_sample_id)) {
    line <- readLines(meta, n = 1)
    if (length(line) == 0) {
      print("Couldn't find normal_sample and tumor_sample lines in Mutect2 VCF meta")
      break
    }
    if (grepl("^##normal_sample", line) && is.null(normal_sample_id)) {
      normal_sample_id <- sub("##normal_sample=", "", line)
    } else if (grepl("^##tumor_sample", line) && is.null(tumor_sample_id)) {
      tumor_sample_id <- sub("##tumor_sample=", "", line)
    }
    if (!is.null(normal_sample_id) && !is.null(tumor_sample_id)) {
      print("Extracted normal_sample and tumor_sample IDs from Mutect2 VCF meta")
      break
    }
  }
  close(meta)
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  print(paste("set", normal_sample_id, "as NORMAL"))
  vcf <- dplyr::rename(vcf, NORMAL = normal_sample_id)
  print(paste("set", tumor_sample_id, "as TUMOR"))
  vcf <- dplyr::rename(vcf, TUMOR = tumor_sample_id)
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       info = parse_info(INFO),
                       normal = parse_format(FORMAT, NORMAL),
                       tumor = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, N = normal, T = tumor, .sep = "_")
  vcf <- dplyr::select(vcf, -normal, -tumor)
  vcf <- tidyr::separate(vcf, N_AD, c("N_REF_COUNT", "N_ALT_COUNT"), ",", convert = TRUE)
  vcf <- dplyr::mutate(vcf, N_VAF = N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT))
  vcf <- tidyr::separate(vcf, T_AD, c("T_REF_COUNT", "T_ALT_COUNT"), ",", convert = TRUE)
  vcf <- dplyr::mutate(vcf, T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT))
  vcf <- readr::type_convert(vcf)
  vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
  vcf
}

parse_nonpaired_mutect2 <- function(x) {
  print(paste0("Processing ", basename(x)))
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  vcf <- dplyr::rename(vcf, TUMOR = names(vcf)[10])
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       info = parse_info(INFO),
                       tumor = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, T = tumor, .sep = "_")
  vcf <- dplyr::select(vcf, -tumor)
  vcf <- tidyr::separate(vcf, MBQ, c("MBQ_REF", "MBQ_ALT"), ",", convert = TRUE)
  vcf <- tidyr::separate(vcf, MMQ, c("MMQ_REF", "MMQ_ALT"), ",", convert = TRUE)
  vcf <- tidyr::separate(vcf, T_AD, c("T_REF_COUNT", "T_ALT_COUNT"), ",", convert = TRUE)
  vcf <- dplyr::mutate(vcf, T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT))
  vcf <- readr::type_convert(vcf)
  vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
  vcf
}

parse_nonpaired_mutect2_or_strelka <- function(x) {
  print(paste0("Processing ", basename(x)))
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  vcf <- dplyr::rename(vcf, TUMOR = names(vcf)[10])
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       info = parse_info(INFO),
                       tumor = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, T = tumor, .sep = "_")
  vcf <- dplyr::select(vcf, -tumor)
  vcf <- tidyr::separate(vcf, T_AD, c("T_REF_COUNT", "T_ALT_COUNT"), ",", convert = TRUE)
  vcf <- dplyr::mutate(vcf, T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT))
  vcf <- readr::type_convert(vcf)
  vcf <- dplyr::group_by(vcf, CHROM, POS, REF, ALT)
  vcf
}

# function required by parse_manta
# splits "[TN]_[PS]R" into "[TN]_[PS]R_REF" and "[TN]_[PS]R_ALT"
# e.g."N_PR" into "N_PR_REF" and "N_PR_ALT"
split_readcounts <- function(x) {
  name <- names(x)
  x$REF <- x[[1]] %>%
    strsplit(",", fixed = TRUE) %>%
    purrr::map(~as.integer(.x[1])) %>%
    purrr::flatten_int()
  x$ALT <- x[[1]] %>%
    strsplit(",", fixed = TRUE) %>%
    purrr::map(~as.integer(.x[2])) %>%
    purrr::flatten_int()
  names(x)[2:3] <- paste0(name, "_", names(x)[2:3])
  x[2:3]
}

parse_paired_manta <- function(x) {
  print(paste0("Processing ", basename(x)))
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  names(vcf)[ncol(vcf)-1] <- "NORMAL"
  names(vcf)[ncol(vcf)] <- "TUMOR"
  # skip empty Manta
  if (nrow(vcf) == 0) {
    return(vcf)
  }
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       # diferrent SVs have different info fields, 
                       # but parse_info manages to populate only relevant fields 
                       # and assign NA to irrelevant ones.
                       info = parse_info(INFO),
                       normal = parse_format(FORMAT, NORMAL),
                       tumor = parse_format(FORMAT, TUMOR))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, N = normal, T = tumor, .sep = "_")
  vcf <- purrr::lmap_at(vcf, dplyr::matches("[TN]_[PS]R", vars = names(vcf)), split_readcounts)
  # make the parser robust against vcfs that do not have an SR value
  if (!("T_SR_REF" %in% names(vcf))) {
    for (col in c("T_SR_REF", "T_SR_ALT", "N_SR_REF", "N_SR_ALT")) {
      vcf[[col]] <- NA
    }
  }
  vcf <- dplyr::mutate(vcf,
                       # to get a rough estimate of the VAF of SVs I followed 
                       # the strategy described in this link:
                       # https://github.com/Illumina/manta/issues/119
                       T_REF_COUNT = T_PR_REF + ifelse(is.na(T_SR_REF), 0, T_SR_REF),
                       T_ALT_COUNT = T_PR_ALT + ifelse(is.na(T_SR_ALT), 0, T_SR_ALT),
                       N_REF_COUNT = N_PR_REF + ifelse(is.na(N_SR_REF), 0, N_SR_REF),
                       N_ALT_COUNT = N_PR_ALT + ifelse(is.na(N_SR_ALT), 0, N_SR_ALT),
                       T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT),
                       N_VAF = N_ALT_COUNT / (N_ALT_COUNT + N_REF_COUNT),
                       # T_DP and N_DP might only be rough estimates of the actual depths
                       T_DP = T_ALT_COUNT + T_REF_COUNT,
                       N_DP = N_ALT_COUNT + N_REF_COUNT)
  vcf <- dplyr::select(vcf, -c(normal, tumor))
  vcf
}

parse_nonpaired_manta <- function(x) {
  print(paste0("Processing ", basename(x)))
  vcf <- readr::read_tsv(x, col_types = vcf_col_parse, comment = "##")
  vcf <- dplyr::rename(vcf, CHROM = `#CHROM`)
  names(vcf)[10] <- "SAMPLE"
  vcf <- dplyr::mutate(vcf, ALT = strsplit(ALT, ",") %>% purrr::map_chr(1),
                       # diferrent SVs have different info fields, 
                       # but parse_info manages to populate only relevant fields 
                       # and assign NA to irrelevant ones.
                       info = parse_info(INFO),
                       format = parse_format(FORMAT, SAMPLE))
  vcf <- tidyr::unnest(vcf, info)
  vcf <- tidyr::unnest(vcf, format)
  vcf <- tidyr::separate(vcf, PR, c("PR_REF", "PR_ALT"), ",", convert = TRUE)
  vcf <- tidyr::separate(vcf, SR, c("SR_REF", "SR_ALT"), ",", convert = TRUE)
  vcf <- dplyr::mutate(vcf,
                       # to get a rough estimate of the VAF of SVs I followed 
                       # the strategy described in this link:
                       # https://github.com/Illumina/manta/issues/119
                       T_REF_COUNT = PR_REF + ifelse(is.na(SR_REF), 0, SR_REF),
                       T_ALT_COUNT = PR_ALT + ifelse(is.na(SR_ALT), 0, SR_ALT),
                       T_VAF = T_ALT_COUNT / (T_ALT_COUNT + T_REF_COUNT),
                       # T_DP might only be a rough estimates of the actual depth
                       T_DP = T_ALT_COUNT + T_REF_COUNT)
  vcf
}
