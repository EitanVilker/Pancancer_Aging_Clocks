#!/usr/bin/env Rscript

# Load required library
library(rmarkdown)
library(glue)
library(doParallel)
source("PanClockHelperFunctions.R")

args = commandArgs(trailingOnly=TRUE)
registerDoParallel(cores = 4)
Sys.setenv(R_MAX_STACK_SIZE = "50000000000")  # Increase to 50000MB
options(expressions = 500000)

# Define the Rmd file and the output directory
rmd_file <- "Many_models_one_cancer_test.Rmd"
if (length(args) > 0){
  test_name <- args[1]
} else { test_name <- "Ridge005" }

base_dir <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/Eitan"
# output_dir <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/"
output_dir <- file.path(base_dir, test_name)
dir.create(output_dir)

# cancer_types <- getEligibleCancers("RNAseq", cutoff=0) # If you want all cancers
cancer_types <- getEligibleCancers("RNAseq", cutoff=250)
cancer_types_methylation <- getEligibleCancers("methylation", cutoff=250)
cancer_types <- cancer_types[cancer_types %in% cancer_types_methylation]

# cancer_types <- c("KIRC")

model_type <- "ridge"
for (cancer_type in cancer_types) {
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
        significance_cutoff=0.05,
        getting_combinations=FALSE,
        output_dir=output_dir,
        test_name=test_name),
      output_file = output_file
    )
    message(paste("Successfully processed:", cancer_type))
  }, error = function(e) {
    # Handle the error
    message(paste("Error occurred with cancer type:", cancer_type))
    message("Error details:", e$message)
  })
}
