# the idea is to source this setup script at the beginning of 
# all rmarkdown documents 
# it will read the project config file and set default options

# load project-specific parameters
config <- yaml::read_yaml(file = 'config.yaml')

# load color specifications:
COLOR_PALETTES <- sapply(
  split(data.table::fread(file.path(config$project_root, "metadata", "colors.csv")),by="category"), 
  function(x) { 
    res <- x[,value]
    names(res) <- x[,term]
    return(res);  
  }, 
  simplify=F)


# set knitr options
knitr::opts_chunk$set(comment = NA, fig.width = 7, fig.height = 4, out.width = '70%',
                      warning = TRUE, error = FALSE, echo = TRUE, message = TRUE,
                      dpi = 100)

# set some other package-specific options
options(ggrepel.max.overlaps = Inf)

options(future.plan = 'sequential', 
        future.globals.maxSize = 8 * 1024 ^ 3)

ggplot2::theme_set(ggplot2::theme_bw(base_size = 12))

# set a random seed
set.seed(993751)

# common constants:
DAY_ORDER <- names(COLOR_PALETTES$day) 
CONDITION_ORDER <- names(COLOR_PALETTES$condition)
CHROM_ORDER <- paste0("chr", c(1:22, "X"))

# store the current CPU and real time
SETUP_TIME <- proc.time()

# source any scripts with commonly used functions
source('R/utilities.R')
