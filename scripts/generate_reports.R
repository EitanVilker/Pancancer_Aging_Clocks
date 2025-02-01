# Load required library
library(rmarkdown)
library(glue)
setwd("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/scripts")

# Define the Rmd file and the output directory
rmd_file <- "Many_models_one_cancer_test.Rmd"
output_dir <- "restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/Multiomics_multiple_models_per_cancer"

# Fetch cancer types dynamically from file names
folder_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/TCGA-GDC/RNAseq/processed_data/filtered"
files <- list.files(folder_path)

# Extract the variable part (e.g., "HNSC") from the filenames
cancer_types <- unique(sub(".*TCGA-(.*)_RNAseq_filtered\\.rds", "\\1", files))

# Print cancer types to verify
print(cancer_types)

# Loop through each cancer type and render the report
#for (cancer_type in cancer_types) {
# Construct the output file path
#  output_file <- file.path(
#    output_dir,
#    paste0(tools::file_path_sans_ext(rmd_file), "_", cancer_type, ".html")
#  )

# Render the R Markdown file
#  rmarkdown::render(
#    input = rmd_file,
#    params = list(cancer_type = cancer_type),
#    output_file = output_file
#  )
#}
for (cancer_type in cancer_types) {
  # Construct the output file path
  output_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(rmd_file), "_", cancer_type, ".html")
  )
  
  # Render the R Markdown file with error handling
  tryCatch({
    rmarkdown::render(
      input = rmd_file,
      params = list(cancer_type = cancer_type),
      output_file = output_file
    )
    message(paste("Successfully processed:", cancer_type))
  }, error = function(e) {
    # Handle the error
    message(paste("Error occurred with cancer type:", cancer_type))
    message("Error details:", e$message)
  })
}
