---
title: "R Notebook"
output: html_notebook
---













<<<<<<< HEAD


=======
``` r
setup1 <- function(layersVector, cancerName="HNSC"){
  paths <- getPathsByLayerAndCancer()
  layerPaths <- list()
  cancerPaths <- list()
  assayNames <- list()
  needsImputation <- list()
  cancerSpecificLayerPaths <- getLayerPathsForSpecificCancer(paths, cancerName)
  
  # Don't run model if missing a layer
  for (layer in layersVector){
    if (length(cancerSpecificLayerPaths[[layer]]) == 0){ return(NULL) }
  }
  
  if ("methylation" %in% layersVector) { 
    layerPaths$methylation <- cancerSpecificLayerPaths$methylation
    assayNames$methylation <- "M-values"
    needsImputation$methylation <- TRUE
  }
  if ("RNAseq" %in% layersVector) { 
    layerPaths$RNAseq <- cancerSpecificLayerPaths$RNAseq_filtered
    assayNames$RNAseq <- "DESeq2_log"
    needsImputation$RNAseq <- TRUE
  }
  if ("miRNA" %in% layersVector) {
    layerPaths$miRNA <- cancerSpecificLayerPaths$miRNA
    assayNames$miRNA <- "DESeq2_log"
    needsImputation$miRNA <- FALSE
  }
  if ("RPPA" %in% layersVector) { 
    layerPaths$RPPA <- cancerSpecificLayerPaths$RPPA
    assayNames$RPPA <- "counts"
    needsImputation$RPPA <- TRUE
  }
  if ("RNAseq_Normal" %in% layersVector) { 
    layerPaths$RNAseq_Normal <- cancerSpecificLayerPaths$RNAseq_NORMAL_filtered
    assayNames$RNAseq_Normal <- "DESeq2_log"
    needsImputation$RNAseq_Normal <- FALSE
  }
  
  return(list(layerPaths=layerPaths, assayNames=assayNames, needsImputation=needsImputation, cancerName=cancerName))
}
```

# Perform second part of setup, including imputing, getting significant features


``` r
setup2 <- function(setup1List, featuresToEnsure=c("Age"), existingExperimentsList=NULL){
  if (is.null(setup1List)) { return(NULL) }
  
  PATH <- file.path( Sys.getenv("PCANAGE","/restricted/projectnb/agedisease/projects/pancancer_aging_pbock")) 
  layerPaths <- setup1List$layerPaths
  assayNames <- setup1List$assayNames
  needsImputation <- setup1List$needsImputation
  cancerName <- setup1List$cancerName
  removeHPVPositive <- FALSE
  # Set up reduced feature set based on association with age
  simplifiedData <- list()
  if (cancerName == "HNSC"){
    featuresToEnsure <- c("Age", "HPV.status")
    removeHPVPositive <- TRUE
  }

  # Get significant features for RNAseq and methylation
  if (!is.null(layerPaths$methylation)) { 
    simplifiedData$methylation <- list(data=readRDS(paste(PATH,"/results/methylation/DMP_analysis/TCGA-",cancerName,"_CpG-Sites_hg19_limma.rds", sep="")), metric="adj.P.Val")
  }
  if (!is.null(layerPaths$RNAseq)) { 
    simplifiedData$RNAseq <- list(data=readRDS(paste(PATH,"/results/RNAseq/DGE_analysis/TCGA-",cancerName,"_filtered_DESeq2_results.rds", sep="")), metric="padj")
  }
  # Using same age associated features as base for normal data
  if (!is.null(layerPaths$RNAseq_Normal)) {
    simplifiedData$RNAseq_Normal <- list(data=readRDS(paste(PATH,"/results/RNAseq/DGE_analysis/TCGA-",cancerName,"_filtered_DESeq2_results.rds", sep="")), metric="padj")
  }
  # if (!is.null(layerPaths$RNAseq_Normal)) { simplifiedData$RNAseq_Normal <- NULL }
  if (!is.null(layerPaths$miRNA)) { simplifiedData$miRNA <- NULL }
  if (!is.null(layerPaths$RPPA)) { simplifiedData$RPPA <- NULL }
  
  # Get list of SummarizedExperiments for each omics layer of the cancer
  if (is.null(existingExperimentsList)){ 
    experimentsList <- getExperimentsList(layerPaths, featuresToEnsure=featuresToEnsure, removeHPVPositive=removeHPVPositive)
  }
  
  # Perform imputation (only if an imputed version is not found)
  if (is.null(existingExperimentsList)) { 
    experimentsList <- imputeMissingValues(experimentsList, assayNames, needsImputation, cancerName) 
  } else { 
    experimentsList <- existingExperimentsList 
  }

  return(list(experimentsList=experimentsList, assayNames=assayNames, needsImputation=needsImputation, simplifiedData=simplifiedData))
}
```

# Wrapper for the two setup functions
>>>>>>> 55f5cc066a29b6b3b73b83b46c1eec722fb59001



<<<<<<< HEAD
=======
# Wrapper for ModelBuilding with default params to make life easier. Just need to pass in vector of omics layer names and the cancer name abbreviation, but can adjust other parameters as desired like setting splitSize to 0.8 for a holdouts test


``` r
standardizedModelBuilding <- function(layersVector, cancerName="HNSC", splitSize=0.999, testOnCompleteData=FALSE, graphTitle="", methodNames=c("glmnet"), seed=43, columnsToKeep=c("Age", "gender", "race", "HPV.status", "submitter_id"), stratifying=FALSE, existingExperimentsList=NULL, iterationCount=NULL, applyBiasCorrection=FALSE, combiningNormal=FALSE, scaleByNormal=FALSE){
  if (splitSize > 0.98){ testOnCompleteData=TRUE}
  modelBuildingInput <- performSetup(layersVector, cancerName=cancerName, existingExperimentsList=existingExperimentsList)
  if (is.null(modelBuildingInput)) { return(NULL) }
  
  # Set up title
  if (nchar(graphTitle) == 0){
    string1 <- paste(layersVector, collapse=", ")
    string2 <- paste(methodNames, collapse=", ")
    if (testOnCompleteData){ string3 <- "Complete Set"}
    else{ string3 <- "Holdouts"}
    graphTitle <- paste(string1, "; ", string2, "; ", string3, sep="")
  }
  
  return(list(ModelBuilding=ModelBuilding(modelBuildingInput, cancerName=cancerName, splitSize=splitSize, testOnCompleteData=testOnCompleteData, graphTitle=graphTitle, methodNames=methodNames, seed=seed, columnsToKeep=columnsToKeep, stratifying=stratifying, iterationCount=iterationCount, applyBiasCorrection=applyBiasCorrection, combiningNormal=combiningNormal, scaleByNormal=scaleByNormal), experimentList=modelBuildingInput$experimentsList))
}
```
>>>>>>> 55f5cc066a29b6b3b73b83b46c1eec722fb59001
