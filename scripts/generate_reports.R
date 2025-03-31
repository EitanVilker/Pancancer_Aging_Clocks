# Load required library
library(rmarkdown)
library(glue)
library(doParallel)
source("PanClockHelperFunctions.R")

registerDoParallel(cores = 4)
Sys.setenv(R_MAX_STACK_SIZE = "50000000000")  # Increase to 50000MB
options(expressions = 500000)

# Define the Rmd file and the output directory
rmd_file <- "Many_models_one_cancer_test.Rmd"
# output_dir <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results"
output_dir <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/Eitan/ElasticNetCV005"

# cancer_types <- getEligibleCancers("RNAseq", cutoff=0) # If you want all cancers
cancer_types <- getEligibleCancers("RNAseq")

model_type <- "ElasticNet"
for (cancer_type in cancer_types) {
  # Construct the output file path
  output_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(rmd_file), "_", cancer_type, "_", model_type, ".html")
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
        test_name="RidgeCV005"),
      output_file = output_file
    )
    message(paste("Successfully processed:", cancer_type))
  }, error = function(e) {
    # Handle the error
    message(paste("Error occurred with cancer type:", cancer_type))
    message("Error details:", e$message)
  })
}
