# Load required library
library(rmarkdown)
library(glue)
#setwd("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/scripts")
setwd("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/apazhern/Pan_Cancer_Aging_Clocks/Pancancer_Aging_Clocks/scripts")
# Define the Rmd file and the output directory
rmd_file <- "Many_models_one_cancer_test.Rmd"
output_dir <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results"

# Fetch cancer types dynamically from file names
folder_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/TCGA-GDC/RNAseq/processed_data/filtered"
files <- list.files(folder_path)

#OGGGG Extract the variable part (e.g., "HNSC") from the filenames


cancer_types <- unique(sub(".*TCGA-(.*)_RNAseq_filtered\\.rds", "\\1", files))


#Temporary One by One:
cancer_types <- c("BRCA", "CHOL", "DLBC", "GBM", "LAML", "LGG", 
                  "LIHC", "OV", "PRAD", "SARC", "THCA", "UCEC")
#cancer_types1 <- c("LGG", "LIHC", "OV", "PRAD", "SARC", "THCA", "UCEC")

cancer_types <- c("GBM") #Methylation still fails to run, will wait till I can run with everything

cancer_list_filtered <- setdiff(cancer_types, cancer_types1)

cancer_types <- cancer_list_filtered

#Problems: LAML is the only cancer missing in RPPA, all cancers and omics exist
#GBM: Problems with miRNA (Report generated w/o miRNA)
#"LGG", "LIHC", "OV", "PRAD", "SARC", "THCA", "UCEC"(Reports generated w/o methylation) 
#"BRCA", "CHOL", "DLBC" (Reports generated w/o (methylation AND RNAseq alone))
#LAML (Problems with RPPA only, All other(RNAseq, miRNA, methylation) individually run fine alone, but don't work when run together
  #mi_RNA & RNAseq (runs fine)
  #mi_RNA & methylation (runs fine)
  #RNAseq & methylation (that's where it breaks)
  #"miRNA","RNAseq","methylation" (that's where it breaks too)
  #Ran with: layer_combinations <- list(c("miRNA"), c("RNAseq"), c("methylation"), c("miRNA","RNAseq"), c("miRNA","methylation"))




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
