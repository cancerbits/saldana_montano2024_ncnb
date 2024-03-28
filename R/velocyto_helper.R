# use this script to generate the barcode files (after filtering and demultiplexing) for the individual libraries
# you need to run this before velocyto

# current working directory should be the one holding the project config file (main project directory)

library(dplyr)

config <- yaml::read_yaml(file = 'config_bob.yaml')

md_path <- file.path(config$out_root, 'input_data', 'full_dataset_metadata_with_UMAP_coordinates.rds')
md <- readRDS(md_path)

# write the barcodes to tsv file in the correct place
# tsv file should have just one column, a barcode example is AAACCCAAGAAACCAT-1
lib_md <- readr::read_csv(file = file.path(config$project_root, 'metadata', 'velocyto_libraries.csv'))
md <- left_join(md, lib_md, by = c('orig.ident' = 'orig_ident'))

# we write the files to the tmp directory on scratch
out_dir <- file.path(config$tmp_root, paste0(config$project_name, '_velocyto'))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (lib_name in unique(md$library_name)) {
  tsv_file <- file.path(out_dir, sprintf('%s_filtered_barcodes.tsv', lib_name))
  filter(md, library_name == lib_name) %>%
    pull(Barcode) %>%
    as.data.frame() %>%
    readr::write_tsv(file = tsv_file, col_names = FALSE)
}
