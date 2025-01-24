# Function to get a list of RangedSummarizedExperiment and SummarizedExperiment objects based on provided paths
# Removes subjects with missing values for inputted list of features
# Get experiment out of list using syntax: experimentList[[index]]
getExperimentsList <- function(paths, featuresToEnsure=c("Age"), removeHPVPositive=FALSE) {
  
  experimentList <- list()
  # Add each experiment
  for (i in 1:length(paths)){
    name <- names(paths)[i]
    experiment <- readRDS(paths[[i]])
    # Remove subjects with certain missing metadata
    for (feature in featuresToEnsure){
      experiment <- experiment[, !is.na(experiment[[feature]])]
    }
    if (removeHPVPositive){ experiment <- experiment[, experiment$HPV.status=="Negative"] }
    experiment$Age <- as.numeric(experiment$Age)
    experimentList[[name]] <- experiment
  }
  return(experimentList)
}

# Function to get a list of vectors of file paths for each cancer and layer
getPathsByLayerAndCancer <- function() {
  
  layerNames <- c("miRNA", "RNAseq_filtered", "methylation", "RPPA", "RNAseq_NORMAL_filtered", "binary-mutation", "SCNA")
  layerPartialPaths <- c("miR", "RNAseq", "methylation", "RPPA", "normal/RNAseq", "mutation", "scna")
  # filteredOrAll <- c(1,1,1,1,1,2,2)
  
  commonPath <- "/restricted/projectnb/agedisease/CBMrepositoryData/TCGA-GDC/"
  filteredPath <- "/processed_data/filtered/TCGA-*_"
  allPath <- "/processed_data/all/TCGA-*"
  
  pathsByLayerAndCancer <- list()
  for (i in 1:length(layerNames)){
    
    name <- layerNames[i]
    
    # Concatenate main path and specific layer
    filePath <- paste(commonPath, layerPartialPaths[i], sep="") 
    
    # Add filtered data to the path if possible, all if not
    if (i < 6){ filePath <- paste(filePath, filteredPath, sep="") }
    else{ filePath <- paste(filePath, allPath, sep="") }
    
    # Get each cancer type for the layer and add to the list
    filePath <- paste(filePath, name, ".rds", sep="")
    pathsByLayerAndCancer[[name]] <- Sys.glob(filePath)
  }
  
  # Return list of vectors of path names, indexed by layer names
  return(pathsByLayerAndCancer)
}

getLayerPathsForSpecificCancer <- function(paths, cancerName){
  layerPaths <- list()
  for (layer in names(paths)){
    layerPaths[[layer]] <- grep(cancerName, paths[[layer]], value=TRUE)
  }
  return(layerPaths)
}

getExperimentByLayerAndCancer <- function(layerName, cancerName){
  revisedLayerName <- layerName
  if (layerName == "RNAseq"){ revisedLayerName <- "RNAseq_filtered" }
  paths <- getPathsByLayerAndCancer()
  cancerPaths <- getLayerPathsForSpecificCancer(paths, cancerName)
  experiment <- getExperimentsList(list(layer=cancerPaths[[revisedLayerName]]))
  return(experiment)
}

getAllCancers <- function(paths){
  PATH <- file.path( Sys.getenv("PCANAGE","/restricted/projectnb/agedisease/projects/pancancer_aging_pbock")) 
  simplifiedData <- list()
  cancerPaths <- list()
  
  for (layer in c("RNAseq_filtered")){
    stringBefore <- "RNAseq/processed_data/filtered/TCGA-"
    stringAfter <- "_RNAseq_filtered."
    
    for (i in 1:length(paths[[layer]])){
      currentPath <- paths[[layer]][i]
      cancerName <- stringr::str_split_i(currentPath, stringBefore, 2)
      cancerName <- stringr::str_split_i(cancerName, stringAfter, 1)
      cancerPaths[["RNAseq"]][cancerName] <- currentPath
      if (cancerName == "CHOL" || cancerName == "DLBC") { simplifiedData[[cancerName]] <- "" }
      else { simplifiedData[[cancerName]] <- list(data=readRDS(paste(PATH,"/results/RNAseq/DGE_analysis/TCGA-",cancerName, "_filtered_DESeq2_results.rds", sep="")), metric="padj") }
    }
  }
  return(list(cancerPaths=cancerPaths, simplifiedData=simplifiedData))
}

isRelativeVarianceHigh <- function(column, threshold=0.25){
  if (!is.numeric(column)){ return(column) }
  return(abs(var(column) / mean(column)) >= threshold)
}

checkIfNameMatch <- function(name, experiment=NULL){
  if (is.null(experiment)){ return(FALSE) } 
  return(any(experiment$submitter_id==name))
}
