# to be run in ccribioinf/dockrstudio:4.0.3-v1
# dockerhub manifest digest: sha256:055a81938337a3fc2a15039476094058f5691e02b529578041af508352264323

addColsToEmptyGr <- function(gr, col_names) {
  if (length(gr) != 0) {
    stop("GRanges is not empty!")
  }
  n_old <- ncol(mcols((gr)))
  n_new <- length(col_names)
  new_names <- c(names(mcols((gr))), col_names)
  elm_meta <- data.frame(matrix(nrow=0,ncol=n_old+n_new))
  colnames(elm_meta) <- new_names
  gr_out <- gr
  gr_out@elementMetadata <- DataFrame(elm_meta) # DataFrame is an S4 invention
  return(gr_out)
}

addWlColumnToGr <- function(gr, wl_gr, col_name) {
  gr@elementMetadata[[col_name]] <- 0
  gr_overlaps <- subjectHits(findOverlaps(wl_gr, gr))
  if (length(gr_overlaps) != 0) {
    gr@elementMetadata[[col_name]][gr_overlaps] <- 1
  }
  return(gr)
}

safe_logical_comparison <- function(x, logic, y) {
  # Check if x is numeric or NA
  if (!is.numeric(x) && !is.na(x)) {
    stop("safe_logical_comparison(x,logic,y): x is not numeric")
  }
  
  # Check if logic is valid
  if (!logic %in% c("l", "g", "le", "ge")) {
    stop("safe_logical_comparison(x,logic,y): logic must be one of 'l', 'g', 'le', or 'ge'")
  }
  
  # Check if y is not missing and is numeric
  if (missing(y) || !is.numeric(y)) {
    stop("safe_logical_comparison(x,logic,y): y is not provided")
  }
  
  # Perform the comparison based on the logic
  if (is.na(x)) {
    return(FALSE)
  } else {
    if (logic == "l") {
      return(x < y)
    } else if (logic == "g") {
      return(x > y)
    } else if (logic == "le") {
      return(x <= y)
    } else if (logic == "ge") {
      return(x >= y)
    }
  }
}

MAPPYACTS_filter <- function(data) {
  library(dplyr)
  library(magrittr)
  if (is(data, "GRanges")) {
    data <- as.data.frame(data)
  }
  # population frequency > 0.01%
  pop_af <- 0.0001
  # a quality < 30
  qual <- 30
  data %>%
    rowwise() %>%
    mutate(
      # population frequency > 0.01%
      MPY_1kG_AF = ifelse("AF" %in% names(data), safe_logical_comparison(AF, "g", pop_af), FALSE),
      MPY_gnomAD_AF = ifelse("gnomAD_AF" %in% names(data), safe_logical_comparison(gnomAD_AF, "g", pop_af), FALSE),
      MPY_gnomADe_AF = ifelse("gnomADe_AF" %in% names(data), safe_logical_comparison(gnomADe_AF, "g", pop_af), FALSE),
      MPY_gnomADg_AF = ifelse("gnomADg_AF" %in% names(data), safe_logical_comparison(gnomADg_AF, "g", pop_af), FALSE),
      MPY_MAX_AF = ifelse("MAX_AF" %in% names(data), safe_logical_comparison(MAX_AF, "g", pop_af), FALSE),
      # a quality < 30
      MPY_QUAL = ifelse(QUAL == ".", FALSE, safe_logical_comparison(as.numeric(QUAL), "l", qual)),
      MPY_GERMQ = ifelse("GERMQ" %in% names(data), safe_logical_comparison(GERMQ, "l", qual), FALSE),
      MPY_ROQ = ifelse("ROQ" %in% names(data), safe_logical_comparison(ROQ, "l", qual), FALSE),
      MPY_SEQQ = ifelse("SEQQ" %in% names(data), safe_logical_comparison(SEQQ, "l", qual), FALSE),
      MPY_STRANDQ = ifelse("STRANDQ" %in% names(data), safe_logical_comparison(STRANDQ, "l", qual), FALSE),
      MPY_STRQ = ifelse("STRQ" %in% names(data), safe_logical_comparison(STRQ, "l", qual), FALSE),
      # intronic variants
      MPY_INT = ifelse("Consequence" %in% names(data), Consequence == "intron_variant", FALSE),
      # synonymous variants
      MPY_SYN = ifelse("Consequence" %in% names(data), Consequence == "synonymous_variant", FALSE),
      # except those with cosmic ID
      MPY_COSV = ifelse("Existing_variation" %in% names(data), ((MPY_INT | MPY_SYN) & grepl("COSV", Existing_variation) & !(MPY_1kG_AF | MPY_gnomAD_AF | MPY_MAX_AF | MPY_QUAL | MPY_GERMQ | MPY_ROQ | MPY_SEQQ | MPY_STRANDQ | MPY_STRQ)), FALSE),
      # SOM: a coverage >= 20x and supported by at least 2 reads in tumor or plasma and no reads in germline
      MPY_SOM = ifelse("N_DP" %in% names(data), safe_logical_comparison(N_DP, "ge", 20) && safe_logical_comparison(T_DP, "ge", 20) && safe_logical_comparison(T_ALT_COUNT, "ge", 2) && N_ALT_COUNT == 0, FALSE)
    ) %>%
    mutate(
      MAPPYACTS = paste(c("1kG_AF","gnomAD_AF","gnomADe_AF","gnomADg_AF","MAX_AF","QUAL","GERMQ","ROQ","SEQQ","STRANDQ","STRQ","INT","SYN","COSV","SOM")[c(MPY_1kG_AF,MPY_gnomAD_AF,MPY_gnomADe_AF,MPY_gnomADg_AF,MPY_MAX_AF,MPY_QUAL,MPY_GERMQ,MPY_ROQ,MPY_SEQQ,MPY_STRANDQ,MPY_STRQ,MPY_INT,MPY_SYN,MPY_COSV,MPY_SOM)], collapse=","),
      MAPPYACTS = case_when(
        MAPPYACTS == "" ~ "PASS",
        MAPPYACTS == "SOM" ~ "PASS,SOM",
        TRUE ~ MAPPYACTS
      )
    ) %>%
    dplyr::select(-starts_with("MPY_"))
}
