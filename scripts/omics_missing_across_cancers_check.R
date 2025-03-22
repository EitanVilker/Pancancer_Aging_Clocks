############
#Checking which Omics Layers are missing in each cancer
############

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

############
length(cancer_types_RPPA)
cancer_types_RPPA
cancer_types_methylation
missing_in_RPPA <- setdiff(cancer_types_methylation, cancer_types_RPPA)
missing_in_RPPA