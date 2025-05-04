#!/usr/bin/env Rscript

### Load required libraries
library(rmarkdown)
library(glue)
library(doParallel)
source("PanClockHelperFunctions.R")

### Set options for running efficient batch job
args = commandArgs(trailingOnly=TRUE)
registerDoParallel(cores = strtoi(Sys.getenv("NSLOTS")) - 1)
Sys.setenv(R_MAX_STACK_SIZE = "50000000000")  # Increase to 50000MB
options(expressions = 500000)

### Define the Rmd file and the output directory
rmd_file <- "Many_models_one_cancer_test.Rmd"
if (length(args) > 1){ 
  test_name <- args[1]
  cancer_type <- args[2] } else { 
  test_name <- "ComboRidge005" 
  cancer_type <- "LGG"
}

cat("Running on cancer:", cancer_type, "\n")


base_dir <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/Eitan"
output_dir <- file.path(base_dir, test_name)
dir.create(output_dir, showWarnings = FALSE)

### Set parameters for models
model_type <- "ElasticNet"
significance_cutoff <- 0.05
getting_combinations <- FALSE
apply_bias_correction <- FALSE
iteration_count <- 1

### Run cancer
# Construct the output file path
output_file <- file.path(
  output_dir,
  paste0(tools::file_path_sans_ext(rmd_file), "_", cancer_type, "_", test_name, ".html")
)

# Render the R Markdown file with error handling
tryCatch({
  rmarkdown::render(
    input = rmd_file,
    params = list(
      cancer_type=cancer_type,
      model_type=model_type,
      significance_cutoff=significance_cutoff,
      getting_combinations=getting_combinations,
      output_dir=output_dir,
      test_name=test_name,
      iteration_count=iteration_count,
      apply_bias_correction=apply_bias_correction),
    output_file = output_file
  )
  message(paste("Successfully processed:", cancer_type))
}, error = function(e) {
  # Handle the error
  message(paste("Error occurred with cancer type:", cancer_type))
  message(paste("Error details:", e$message))
})

writeUltimateSummaryTable(outputDir=paste0(output_dir, "/"), suffix=paste0("_", model_type, "_omics_summary.csv"))
