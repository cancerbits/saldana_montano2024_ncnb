# this top-level script may be used to run (knit) the individual
# Rmd files

# set up output path for the reports
config <- yaml::read_yaml("config.yaml")

# # example 1: knit one specific file
# rmarkdown::render(input = 'Rmd/example.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())
# 
# # example 2: knit one specific file and pass parameters
# rmarkdown::render(input = 'Rmd/example_with_parameters.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env(),
#                   params = list(mean = 100, sd = 100),
#                   output_file = 'example_with_parameters_100_100')


# report_type <- "atac"
# rmd_dir <- file.path("Rmd", report_type)
# report_dir <- file.path(config$out_root, "results", "revision", report_type, "reports")
# rmds <- list.files(rmd_dir, pattern=".Rmd")
# dir.create(report_dir, showWarnings=FALSE)
# for(rmd in rmds) {
#   rmarkdown::render(input = file.path(rmd_dir, rmd),
#                            output_dir = report_dir,
#                            knit_root_dir = config$project_root,
#                            envir = new.env())
# }

report_type <- "scrna"
rmd_dir <- file.path("Rmd", report_type)
report_dir <- file.path(config$out_root, "results", "revision", report_type, "reports")
rmds <- list.files(rmd_dir, pattern=".Rmd")
dir.create(report_dir, showWarnings=FALSE)
for(rmd in rmds[3]) {
  rmarkdown::render(input = file.path(rmd_dir, rmd),
                    output_dir = report_dir,
                    knit_root_dir = config$project_root,
                    envir = new.env())
}
