---
title: "Gene Pathway Enrichment Analysis"
author: "Yousry"
date: "2025-02-26"
output: html_document
params:
  cancer_types: ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
  # Modify this list as needed
---



## Load Required Libraries


``` r
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "fgsea", "msigdbr"))
# 
# if (!requireNamespace("data.table", quietly = TRUE)) {
#   install.packages("data.table")
# }

# Load libraries
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(fgsea)
library(msigdbr)
knitr::knit("../Analysis.Rmd", output = "Analysis.R")
```

```
##   |                                                                                                                    |                                                                                                            |   0%  |                                                                                                                    |.....                                                                                                       |   5% [unnamed-chunk-9]   |                                                                                                                    |..........                                                                                                  |  10%                     |                                                                                                                    |...............                                                                                             |  14% [unnamed-chunk-10]  |                                                                                                                    |.....................                                                                                       |  19%                     |                                                                                                                    |..........................                                                                                  |  24% [unnamed-chunk-11]  |                                                                                                                    |...............................                                                                             |  29%                     |                                                                                                                    |....................................                                                                        |  33% [unnamed-chunk-12]  |                                                                                                                    |.........................................                                                                   |  38%                     |                                                                                                                    |..............................................                                                              |  43% [unnamed-chunk-13]  |                                                                                                                    |...................................................                                                         |  48%                     |                                                                                                                    |.........................................................                                                   |  52% [unnamed-chunk-14]  |                                                                                                                    |..............................................................                                              |  57%                     |                                                                                                                    |...................................................................                                         |  62% [unnamed-chunk-15]  |                                                                                                                    |........................................................................                                    |  67%                     |                                                                                                                    |.............................................................................                               |  71% [unnamed-chunk-16]  |                                                                                                                    |..................................................................................                          |  76%                     |                                                                                                                    |.......................................................................................                     |  81% [unnamed-chunk-17]  |                                                                                                                    |.............................................................................................               |  86%                     |                                                                                                                    |..................................................................................................          |  90% [unnamed-chunk-18]  |                                                                                                                    |.......................................................................................................     |  95%                     |                                                                                                                    |............................................................................................................| 100% [unnamed-chunk-19]
```

```
## [1] "Analysis.R"
```

## Define function to process GTF annotations file


``` r
process_gtf_annotations <- function(gtf_file, annotation_output) {
  print("🔍 Processing GTF annotation file...")

  # Check if the file is compressed and load accordingly
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
  dir.create(dirname(annotation_output), recursive = TRUE, showWarnings = FALSE)
  fwrite(gene_map, annotation_output, row.names = FALSE)

  print(paste("✅ Annotation file saved:", annotation_output))
}

# usage:
# gtf_file <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/scripts/dataset/gene_annotations/gencode.v47.chr_patch_hapl_scaff.annotation.gtf"
# annotation_output <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/features/gene_annotation.csv"
# process_gtf_annotations(gtf_file, annotation_output)
```

## Define function for cancer-type specific RNA-seq models Gene_name vs weights csv


``` r
process_gene_weights <- function(cancer_type, weights_file=NULL, modelOutput=NULL, annotation_file="/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/features/gene_annotation.csv", output_dir="/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/features") {
  
  print(paste0("🔍 Processing gene weights for ", cancer_type))

  # Define dynamic file paths based on cancer type
  if (is.null(weights_file)){
    weights_file <- paste0("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/", cancer_type, "_elastic_net_omics_combinations_model_weights.rds")
  }
  output_path <- file.path(output_dir, paste0(cancer_type, "_RNAseq_weights_named.csv"))

  # Load the gene annotation file
  if (is.null(modelOutput)){ gene_map <- fread(annotation_file)}
  else{ geneMap <- getGeneMap(modelOutput) }

  # Check if weights file exists
  if (!file.exists(weights_file)) {
    print(paste0("❌ Error: Weights file not found for ", cancer_type, " at ", weights_file))
    return(NULL)
  }
  
  # Load the filtered weights dataset
  model_weights <- readRDS(weights_file)
  
  if (!paste0(cancer_type, "_RNAseq") %in% names(model_weights)) {
    print(paste0("❌ Error: No weights found for ", cancer_type))
    return(NULL)
  }
  
  cancer_data <- model_weights[[paste0(cancer_type, "_RNAseq")]]
  
  # Ensure it is a data.table
  cancer_data <- as.data.table(cancer_data)

  # Check column names
  print("🔍 Checking column names in cancer data:")
  print(colnames(cancer_data))

  # Remove version numbers from Ensembl IDs
  cancer_data[, Feature := sub("\\..*", "", Feature)]

  # Remove near-zero weight features
  threshold <- 1e-8
  cancer_filtered <- cancer_data[abs(Weight) > threshold]

  # Remove categorical features (Keep only Ensembl Gene IDs)
  cancer_filtered <- cancer_filtered[grepl("^ENSG", Feature)]

  # Merge with annotation file to convert Ensembl IDs to Gene Names
  cancer_named <- merge(cancer_filtered, gene_map, by.x = "Feature", by.y = "Ensembl_ID", all.x = TRUE)

  # Replace missing gene names with original Ensembl IDs
  cancer_named[, Gene_Name := ifelse(is.na(Gene_Name) | Gene_Name == "", Feature, Gene_Name)]

  # Save the final cleaned file
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(cancer_named[, .(Gene_Name, Weight)], output_path, row.names = FALSE)

  print(paste("Final gene-mapped file saved at:", output_path))
}
```

## Define Pathway Enrichment Analysis Function


``` r
run_pathway_analysis <- function(cancer_type, input_file=NULL, output_dir=NULL) {
  print(paste0("🔍 Running pathway enrichment for ", cancer_type))
  if (is.null(input_file)){
    input_file <- paste0("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/features/", cancer_type, "_RNAseq_weights_named.csv")
  }
  
  if (!file.exists(input_file)) {
    print(paste0("❌ Error: File not found for ", cancer_type))
    return(NULL)
  }
  
  genes <- fread(input_file)

  # Remove Ensembl version numbers (if present)
  genes[, Gene_Name := sub("\\..*", "", Gene_Name)]

  # Convert Gene Symbols to Entrez IDs
  gene_symbols <- genes$Gene_Name
  entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Merge Entrez IDs with gene weights
  genes_mapped <- merge(genes, entrez_ids, by.x = "Gene_Name", by.y = "SYMBOL", all.x = TRUE)
  genes_mapped <- genes_mapped[!is.na(ENTREZID)]  # Remove unmapped genes

  cat("✅ Successfully mapped", nrow(genes_mapped), "out of", nrow(genes), "genes for", cancer_type, "\n")

  # Define Output Directory
  if (is.null(output_dir)){
    output_dir <- paste0("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/pathways/", cancer_type)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Function to Save ORA Results
  save_results <- function(results, name) {
    if (!is.null(results) && nrow(as.data.frame(results@result)) > 0) {
      results_dt <- as.data.table(results@result)
      fwrite(results_dt, file.path(output_dir, paste0(name, "_raw.csv")))
      
      print(paste0("✅ ", name, " analysis completed with ", nrow(results_dt), " pathways found."))
      
      raw_significant <- results_dt[pvalue < 0.05]
      if (nrow(raw_significant) > 0) {
        fwrite(raw_significant, file.path(output_dir, paste0(name, "_raw_significant.csv")))
        print(paste0("✅ ", nrow(raw_significant), " pathways with p-value < 0.05 found in ", name, "."))
      } else {
        print(paste0("⚠️ No pathways with p-value < 0.05 found in ", name, "."))
      }

      filtered_results <- results_dt[p.adjust < 0.05]
      if (nrow(filtered_results) > 0) {
        fwrite(filtered_results, file.path(output_dir, paste0(name, "_filtered.csv")))
        print(paste0("✅ ", nrow(filtered_results), " pathways with adjusted p-value < 0.05 found in ", name, "."))
      } else {
        print(paste0("⚠️ No significant adjusted pathways found in ", name, "."))
      }
    } else {
      print(paste0("❌ No pathways found in ", name, "."))
    }
  }

  # Function to Save GSEA Results
  save_results_gsea <- function(results, name) {
    if (!is.null(results) && nrow(results) > 0) {
      results_dt <- as.data.table(results)
      fwrite(results_dt, file.path(output_dir, paste0(name, "_raw.csv")))
      
      print(paste0("✅ ", name, " GSEA analysis completed with ", nrow(results_dt), " pathways found."))

      raw_significant <- results_dt[pval < 0.05]
      if (nrow(raw_significant) > 0) {
        fwrite(raw_significant, file.path(output_dir, paste0(name, "_raw_significant.csv")))
        print(paste0("✅ ", nrow(raw_significant), " pathways with p-value < 0.05 found in ", name, " GSEA."))
      } else {
        print(paste0("⚠️ No pathways with p-value < 0.05 found in ", name, " GSEA."))
      }

      filtered_results <- results_dt[padj < 0.05]
      if (nrow(filtered_results) > 0) {
        fwrite(filtered_results, file.path(output_dir, paste0(name, "_filtered.csv")))
        print(paste0("✅ ", nrow(filtered_results), " pathways with adjusted p-value < 0.05 found in ", name, " GSEA."))
      } else {
        print(paste0("⚠️ No significant adjusted pathways found in ", name, " GSEA."))
      }
    } else {
      print(paste0("❌ No pathways found in ", name, " GSEA."))
    }
  }

  # **KEGG Pathway Enrichment Analysis (ORA)**
  kegg_results <- enrichKEGG(
    gene = genes_mapped$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 1.0
  )
  save_results(kegg_results, "KEGG")

  # **Reactome Pathway Enrichment Analysis (ORA)**
  reactome_results <- enrichPathway(
    gene = genes_mapped$ENTREZID,
    organism = "human",
    pvalueCutoff = 1.0
  )
  save_results(reactome_results, "Reactome")

  # **Gene Ontology (GO) Enrichment Analysis (ORA)**
  go_results_bp <- enrichGO(
    gene = genes_mapped$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pvalueCutoff = 1.0,
    readable = TRUE
  )
  save_results(go_results_bp, "GO")

  
  # **Gene Set Enrichment Analysis (GSEA)**
  
  # Create named vector for GSEA ranking
  gene_ranks <- setNames(genes_mapped$Weight, genes_mapped$ENTREZID)
  
  # Remove missing values
  gene_ranks <- na.omit(gene_ranks)
  
  # Order gene ranks for GSEA
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  # Check if there are enough genes for GSEA
  if (length(gene_ranks) < 10) {
    print("⚠️ Too few genes for GSEA. Skipping GSEA for this cancer type.")
    return(NULL)
  }
  
  # Load MSigDB gene sets
  hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  c2_pathways <- msigdbr(species = "Homo sapiens", category = "C2")
  
  # Convert gene symbols to Entrez IDs in pathways
  convert_to_entrez <- function(msigdb_data) {
    mapped <- bitr(msigdb_data$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    merged <- merge(msigdb_data, mapped, by.x = "gene_symbol", by.y = "SYMBOL", all.x = TRUE)
  
    # Remove rows with missing ENTREZID values
    merged <- na.omit(merged)
  
    # Ensure ENTREZID is treated as a character (to avoid factor conversion issues)
    merged$ENTREZID <- as.character(merged$ENTREZID)
  
    # Split into pathway-specific lists
    split(merged$ENTREZID, merged$gs_name)
  }
  
  # Convert pathways
  pathways_h <- convert_to_entrez(hallmark)
  pathways_c2 <- convert_to_entrez(c2_pathways)
  
  # Ensure pathways contain genes in gene_ranks and have enough members
  pathways_h <- pathways_h[lengths(pathways_h) > 5 & lengths(pathways_h) <= 500]
  pathways_c2 <- pathways_c2[lengths(pathways_c2) > 5 & lengths(pathways_c2) <= 500]
  
  # Check if valid pathways exist
  if (length(pathways_h) == 0) {
    print("❌ No valid Hallmark pathways available for GSEA.")
  } 
  if (length(pathways_c2) == 0) {
    print("❌ No valid C2 pathways available for GSEA.")
  }
  if (length(pathways_h) == 0 && length(pathways_c2) == 0) {
    return(NULL)
  }
  
  # Run FGSEA Simple
  if (length(pathways_h) > 0) {
    fgsea_results_h <- fgseaSimple(
      pathways = pathways_h, 
      stats = gene_ranks, 
      minSize = 10, 
      maxSize = 500, 
      nperm = 1000  # Lower permutation count for faster runtime
    )
    save_results_gsea(fgsea_results_h, "GSEA_H")
  }
  
  if (length(pathways_c2) > 0) {
    fgsea_results_c2 <- fgseaSimple(
      pathways = pathways_c2, 
      stats = gene_ranks, 
      minSize = 10, 
      maxSize = 500, 
      nperm = 1000  # Lower permutation count for faster runtime
    )
    save_results_gsea(fgsea_results_c2, "GSEA_C2")
  }
  
  # ** Done **
  # Summary
  print(paste0("🎉 Pathway enrichment analysis completed for ", cancer_type, "!"))

}
```

## Run Pathways Enrichment Analysis Function


``` r
if (FALSE){
  cancer_types <- params$cancer_types  # List of cancer types
  
  annotation_file <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/features/gene_annotation.csv"
  
  # Process Gene Weights for Each Cancer Type
  for (cancer in cancer_types) {
    process_gene_weights(cancer, annotation_file)
  }
  
  # Run the Pathway Enrichment Analysis
  for (cancer in cancer_types) {
    run_pathway_analysis(cancer)
  }
}
```

## Summary


``` r
if (FALSE){
  # Initialize a list to store significant pathways for each cancer type
  significant_pathways_list <- list()
  
  # Define the base directory where results are stored
  base_path <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/pathways"
  
  # Loop through each cancer type
  for (cancer in params$cancer_types) {
    # Define paths to adjusted significant pathways
    files <- list(
      KEGG = file.path(base_path, cancer, "KEGG_filtered.csv"),
      Reactome = file.path(base_path, cancer, "Reactome_filtered.csv"),
      GO = file.path(base_path, cancer, "GO_filtered.csv"),
      GSEA_H = file.path(base_path, cancer, "GSEA_H_filtered.csv"),
      GSEA_C2 = file.path(base_path, cancer, "GSEA_C2_filtered.csv")
    )
  
    # Store significant pathways for this cancer type
    pathways <- list()
    
    for (type in names(files)) {
      file <- files[[type]]
      
      if (file.exists(file)) {
        data <- fread(file, select = "ID")
        if (nrow(data) > 0) {
          # Store pathways along with their type
          pathways <- append(pathways, lapply(data$ID, function(id) list(Pathway = id, Pathway_Type = type)))
        }
      } else {
        print(paste("⚠️ Missing:", file))
      }
    }
  
    # Store unique pathways for this cancer if any exist
    if (length(pathways) > 0) {
      significant_pathways_list[[cancer]] <- rbindlist(pathways, fill = TRUE)
      significant_pathways_list[[cancer]][, Cancer_Type := cancer]
    }
  }
  
  # Check if any pathways were found
  if (length(significant_pathways_list) == 0) {
    print("⚠️ No significant pathways found in any cancer type.")
  } else {
    # Convert list to long format table
    pathway_counts <- rbindlist(significant_pathways_list, fill = TRUE)
  
    # Count occurrences of each pathway across cancer types
    pathway_summary <- pathway_counts[, .(
      Pathway_Type = unique(Pathway_Type),
      Cancer_Types = paste(Cancer_Type, collapse = ", "),
      Count = .N
    ), by = Pathway]
  
    # Sort pathways by Count (greatest to smallest)
    pathway_summary <- pathway_summary[order(-Count)]
  
    # Save results
    output_path <- file.path(base_path, "common_pathways.csv")
    fwrite(pathway_summary, output_path)
  
    # Print summary
    print("✅ Pathways identified across cancer types (sorted by count):")
    print(pathway_summary)
  
    # If running interactively, display results in a table
    if (interactive()) {
      library(DT)
      datatable(pathway_summary)
    }
  }
}
```


``` r
# print("Pathway enrichment analysis completed for all specified cancer types!")
```
