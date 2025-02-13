# ---- PART 1: Process Annotations GTF File ----

# Load necessary libraries
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# Define path to the downloaded GTF file
gtf_file <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/scripts/dataset/gene_annotations/gencode.v47.chr_patch_hapl_scaff.annotation.gtf"

# Check if the file is compressed
if (grepl("\\.gz$", gtf_file)) {
  gtf_data <- fread(cmd = paste("zcat", gtf_file), sep = "\t", header = FALSE)
} else {
  # Remove comment lines manually before loading
  gtf_data <- fread(cmd = paste("grep -v '^#'", gtf_file), sep = "\t", header = FALSE)
}

# Extract relevant columns where feature type is "gene"
gtf_data <- gtf_data[V3 == "gene", .(V9)]

# Extract Ensembl Gene ID and HGNC Gene Name using regex
gtf_data[, Ensembl_ID := sub(".*gene_id \"([^\"]+)\".*", "\\1", V9)]
gtf_data[, Gene_Name := sub(".*gene_name \"([^\"]+)\".*", "\\1", V9)]

# Remove version numbers from Ensembl IDs (e.g., ENSG00000101040.19 → ENSG00000101040)
gtf_data[, Ensembl_ID := sub("\\..*", "", Ensembl_ID)]

# Keep only unique mappings
gene_map <- unique(gtf_data[, .(Ensembl_ID, Gene_Name)])

# Save processed annotation file
annotation_file <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/Features/gene_annotation.csv"
fwrite(gene_map, annotation_file, row.names = FALSE)

print(paste("✅ Annotation file saved:", annotation_file))


# ---- PART 2: Load and Filter Gene Weights ----

# Load necessary libraries
library(data.table)

# Load annotation file
gene_map <- fread(annotation_file)

# Load the filtered weights dataset
HNSC_elastic_net_omics_combinations_model_weights <- readRDS(
  "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/HNSC_elastic_net_omics_combinations_model_weights.rds"
)
HNSC_RNAseq <- HNSC_elastic_net_omics_combinations_model_weights$HNSC_RNAseq

# Ensure HNSC_RNAseq is a data.table
HNSC_RNAseq <- as.data.table(HNSC_RNAseq)

# Check column names before proceeding
print("🔍 Checking column names in HNSC_RNAseq:")
print(colnames(HNSC_RNAseq))

# Remove version numbers from Ensembl IDs in HNSC_RNAseq
HNSC_RNAseq[, Feature := sub("\\..*", "", Feature)]

# Remove near-zero weight features
threshold <- 1e-8
HNSC_RNAseq_filtered <- HNSC_RNAseq[abs(Weight) > threshold]

# Remove categorical features (Keep only Ensembl Gene IDs)
HNSC_RNAseq_filtered <- HNSC_RNAseq_filtered[grepl("^ENSG", Feature)]

# Merge with annotation file to convert Ensembl IDs to Gene Names
HNSC_RNAseq_named <- merge(HNSC_RNAseq_filtered, gene_map, by.x = "Feature", by.y = "Ensembl_ID", all.x = TRUE)

# Replace missing gene names with original Ensembl IDs
HNSC_RNAseq_named[, Gene_Name := ifelse(is.na(Gene_Name) | Gene_Name == "", Feature, Gene_Name)]

# Define output path for the final file
output_path <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/Features/HNSC_RNAseq_weights_named.csv"

# Save the final cleaned file
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(HNSC_RNAseq_named[, .(Gene_Name, Weight)], output_path, row.names = FALSE)

# Confirm the file is saved
print(paste("🎉 Final gene-mapped file saved at:", output_path))