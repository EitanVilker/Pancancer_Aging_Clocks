# --- Install Necessary Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor and necessary libraries
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "fgsea", "msigdbr"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

# --- Load Libraries ---
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(fgsea)
library(msigdbr)

# --- Load Gene List ---
genes <- fread("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/Features/HNSC_RNAseq_weights_named.csv")

# Remove Ensembl version numbers (if present)
genes[, Gene_Name := sub("\\..*", "", Gene_Name)]

# Convert Gene Symbols to Entrez IDs
gene_symbols <- genes$Gene_Name
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Merge Entrez IDs with gene weights
genes_mapped <- merge(genes, entrez_ids, by.x = "Gene_Name", by.y = "SYMBOL", all.x = TRUE)
genes_mapped <- genes_mapped[!is.na(ENTREZID)]  # Remove unmapped genes

# Check gene mapping success
cat("✅ Successfully mapped", nrow(genes_mapped), "out of", nrow(genes), "genes.\n")

# --- Define Output Directory ---
output_dir <- "/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/pathways"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Function to Save Results (Raw, p.value < 0.05, and p.adjust < 0.05) ---
save_results <- function(results, name) {
  if (!is.null(results) && nrow(as.data.frame(results@result)) > 0) {
    results_dt <- as.data.table(results@result)
    
    # Save raw results (all pathways)
    fwrite(results_dt, file.path(output_dir, paste0(name, "_raw.csv")))
    print(paste0("✅ Raw ", name, " results saved."))
    
    # Save pathways with raw p-value < 0.05
    raw_significant <- results_dt[pvalue < 0.05]
    if (nrow(raw_significant) > 0) {
      fwrite(raw_significant, file.path(output_dir, paste0(name, "_raw_significant.csv")))
      print(paste0("✅ Raw significant (p < 0.05) ", name, " results saved."))
    }
    
    # Save adjusted significant results (p.adjust < 0.05)
    filtered_results <- results_dt[p.adjust < 0.05]
    if (nrow(filtered_results) > 0) {
      fwrite(filtered_results, file.path(output_dir, paste0(name, "_filtered.csv")))
      print(paste0("✅ Filtered ", name, " results saved."))
    } else {
      print(paste0("⚠️ No significant adjusted pathways found in ", name, " (p.adjust < 0.05)."))
    }
  } else {
    print(paste0("❌ No pathways found in ", name, "."))
  }
}

# --- KEGG Pathway Enrichment Analysis (ORA) ---
kegg_results <- enrichKEGG(
  gene = genes_mapped$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 1.0  # Save all results
)
save_results(kegg_results, "KEGG")

# --- Reactome Pathway Enrichment (ORA) ---
reactome_results <- enrichPathway(
  gene = genes_mapped$ENTREZID,
  organism = "human",
  pvalueCutoff = 1.0  # Save all results
)
save_results(reactome_results, "Reactome")

# --- Gene Ontology (GO) Enrichment (ORA) ---
go_results_bp <- enrichGO(
  gene = genes_mapped$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # Biological Process
  pvalueCutoff = 1.0,  # Save all results
  readable = TRUE
)
save_results(go_results_bp, "GO")

# --- Gene Set Enrichment Analysis (GSEA) with MSigDB ---
# Prepare ranked list
gene_ranks <- setNames(genes_mapped$Weight, genes_mapped$ENTREZID)

# Debugging: Print top-ranked genes
print("🔍 Checking top-ranked genes for GSEA:")
print(head(gene_ranks, 10))

# Load multiple pathway categories
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
c2_pathways <- msigdbr(species = "Homo sapiens", category = "C2")

# Convert to list format
pathways_h <- split(hallmark$gene_symbol, hallmark$gs_name)
pathways_c2 <- split(c2_pathways$gene_symbol, c2_pathways$gs_name)

# Run FGSEA Multilevel for Hallmark and C2 pathways
fgsea_results_h <- fgseaMultilevel(pathways = pathways_h, stats = gene_ranks)
fgsea_results_c2 <- fgseaMultilevel(pathways = pathways_c2, stats = gene_ranks)

# --- Function to Save GSEA Results (Standardized with ORA) ---
save_results_gsea <- function(results, name) {
  if (!is.null(results) && nrow(results) > 0) {
    results_dt <- as.data.table(results)
    
    # Save raw results (all pathways)
    fwrite(results_dt, file.path(output_dir, paste0(name, "_raw.csv")))
    print(paste0("✅ Raw ", name, " results saved."))
    
    # Save pathways with raw p-value < 0.05
    raw_significant <- results_dt[pval < 0.05]
    if (nrow(raw_significant) > 0) {
      fwrite(raw_significant, file.path(output_dir, paste0(name, "_raw_significant.csv")))
      print(paste0("✅ Raw significant (p < 0.05) ", name, " results saved."))
    }
    
    # Save adjusted significant results (p.adjust < 0.05)
    filtered_results <- results_dt[padj < 0.05]
    if (nrow(filtered_results) > 0) {
      fwrite(filtered_results, file.path(output_dir, paste0(name, "_filtered.csv")))
      print(paste0("✅ Filtered ", name, " results saved."))
    } else {
      print(paste0("⚠️ No significant adjusted pathways found in ", name, " (p.adjust < 0.05)."))
    }
  } else {
    print(paste0("❌ No pathways found in ", name, "."))
  }
}

# Save GSEA results
save_results_gsea(fgsea_results_h, "GSEA_H")
save_results_gsea(fgsea_results_c2, "GSEA_C2")

# Debugging: Check if any GSEA pathways are found
print(paste0("🔍 Number of pathways found in raw GSEA results (HALLMARK): ", nrow(fgsea_results_h)))
print(paste0("🔍 Number of pathways found in raw GSEA results (C2): ", nrow(fgsea_results_c2)))

print("🎉 Comprehensive pathway enrichment analysis completed and saved!")