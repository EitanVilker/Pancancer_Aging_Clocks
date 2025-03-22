#Put all the csv's together across cancers

###### Get the cancers present at every Omics Layer ############

folder_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/TCGA-GDC/RNAseq/processed_data/filtered"
files <- list.files(folder_path)
cancer_types_RNAseq <- unique(sub(".*TCGA-(.*)_RNAseq_filtered\\.rds", "\\1", files))

folder_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/TCGA-GDC/miR/processed_data/filtered"
files <- list.files(folder_path)
cancer_types_miRNA <- unique(sub(".*TCGA-(.*)_miRNA\\.rds", "\\1", files))

folder_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/TCGA-GDC/methylation/processed_data/filtered"
files <- list.files(folder_path)
cancer_types_methylation <- unique(sub(".*TCGA-(.*)_methylation\\.rds", "\\1", files))

folder_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/TCGA-GDC/RPPA/processed_data/filtered"
files <- list.files(folder_path)
cancer_types_RPPA <- unique(sub(".*TCGA-(.*)_RPPA\\.rds", "\\1", files))

############ Defining File Paths for csv files #######

file_paths_RNAseq <- paste0("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/", 
                           cancer_types_RNAseq, 
                    "_elastic_net_omics_combinations.csv")

# Read and combine all CSVs
RNAseq_tables_list <- lapply(file_paths_RNAseq, function(file) {
  if (file.exists(file)) {
    read.csv(file, row.names = NULL)  # Adjust row.names if needed
  } else {
    warning(paste("File not found:", file))
    return(NULL)
  }
})

file_paths_miRNA <- paste0("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/", 
                           cancer_types_miRNA, 
                           "_elastic_net_omics_combinations.csv")

# Read and combine all CSVs
miRNA_tables_list <- lapply(file_paths_miRNA, function(file) {
  if (file.exists(file)) {
    read.csv(file, row.names = NULL)  # Adjust row.names if needed
  } else {
    warning(paste("File not found:", file))
    return(NULL)
  }
})

file_paths_methylation <- paste0("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/", 
                                 cancer_types_methylation, 
                           "_elastic_net_omics_combinations.csv")

methylation_tables_list <- lapply(file_paths_methylation, function(file) {
  if (file.exists(file)) {
    read.csv(file, row.names = NULL)  # Adjust row.names if needed
  } else {
    warning(paste("File not found:", file))
    return(NULL)
  }
})

file_paths_RPPA <- paste0("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/", 
                           cancer_types_RPPA, 
                           "_elastic_net_omics_combinations.csv")

RPPA_tables_list <- lapply(file_paths_RPPA, function(file) {
  if (file.exists(file)) {
    read.csv(file, row.names = NULL)  # Adjust row.names if needed
  } else {
    warning(paste("File not found:", file))
    return(NULL)
  }
})
################### Combine all cancers in one omics ###########################
# Combine all dataframes into one
RNAseq_tables_combined <- do.call(rbind, RNAseq_tables_list)

# Combine all dataframes into one
miRNA_tables_combined <- do.call(rbind, miRNA_tables_list)

# Combine all dataframes into one
methylation_tables_combined <- do.call(rbind, methylation_tables_list)

# Combine all dataframes into one
RPPA_tables_combined <- do.call(rbind, RPPA_tables_list)

################### Combine all omics in one table (Ultimate Summary Table) ###########################
combined_tables <- rbind(RNAseq_tables_combined, miRNA_tables_combined, 
                         methylation_tables_combined, RPPA_tables_combined)

################### Adding the number of features used #######################################
cancer_types <- cancer_types_RNAseq

file_paths <- paste0("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/", 
                     cancer_types, 
                     "_elastic_net_omics_combinations_model_weights.rds")

# Read all .rds files into a list
RDS_of_all_cancers <- lapply(file_paths, function(file) {
  if (file.exists(file)) {
    readRDS(file)  # Read .rds file
  } else {
    warning(paste("File not found:", file))
    return(NULL)
  }
})

#BRCA missing RDS, might need to re-run

#Naming the list elements based on cancer types
names(RDS_of_all_cancers) <- cancer_types

# Flatten the RDS_of_all_cancers list
RDS_of_all_cancers_flat <- list()

# Loop over each cancer type in the original list
for (cancer_type in names(RDS_of_all_cancers)) {
  for (layer_combination in names(RDS_of_all_cancers[[cancer_type]])) {
    # Create a single key: "ACC_miRNA" instead of RDS_of_all_cancers$ACC$ACC_miRNA
    table_name <- layer_combination
    
    # Assign to the new flat list
    RDS_of_all_cancers_flat[[table_name]] <- RDS_of_all_cancers[[cancer_type]][[layer_combination]]
  }
}

RDS_of_all_cancers <- RDS_of_all_cancers_flat

# Ensure the new columns exist in combined_tables
combined_tables$features <- NA
combined_tables$non_zero_weight_features <- NA
combined_tables$zero_weight_features <- NA

# Loop through each row of combined_tables
for (i in seq_len(nrow(combined_tables))) {
  # Construct the correct table name from combined_tables
  table_name <- paste0(combined_tables$cancer_type[i], "_", combined_tables$layer_combination[i])
  
  # Check if the constructed table_name exists in RDS_of_all_cancers
  if (!is.null(RDS_of_all_cancers[[table_name]])) {
    # Extract the corresponding table
    table <- RDS_of_all_cancers[[table_name]]
    
    # Check if it's a dataframe
    if (is.data.frame(table)) {
      # Count total features
      total_features <- nrow(table)
      
      # Count non-zero and zero-weight features
      non_zero_features <- sum(table$Weight != 0, na.rm = TRUE)
      zero_features <- sum(table$Weight == 0, na.rm = TRUE)
      
      # Assign values to combined_tables
      combined_tables$features[i] <- total_features
      combined_tables$non_zero_weight_features[i] <- non_zero_features
      combined_tables$zero_weight_features[i] <- zero_features
    } else {
      warning(paste("Expected a dataframe at:", table_name))
    }
  } else {
    warning(paste("Data not found for:", table_name))
  }
}

#Moving the feature columns to the begining instead of the end
combined_tables <- combined_tables[, c(1:5, ncol(combined_tables), 6:(ncol(combined_tables) - 1))]
combined_tables <- combined_tables[, c(1:5, ncol(combined_tables), 6:(ncol(combined_tables) - 1))]
combined_tables <- combined_tables[, c(1:5, ncol(combined_tables), 6:(ncol(combined_tables) - 1))]

#Getting Rid of the duplicates and just keeping one copy of them
unique_summary_table <- combined_tables %>%
  distinct(test_name, cancer_type, layer_combination, .keep_all = TRUE)

# Standardize layer_combination
unique_summary_table$layer_combination <- sapply(strsplit(unique_summary_table$layer_combination, "_"), function(x) {
  paste(sort(x), collapse = "_")
})

################## Saving the combined table ################################################
file_path <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/Ultimate_summary_table.csv"
write.csv(unique_summary_table, file_path, row.names = FALSE)
