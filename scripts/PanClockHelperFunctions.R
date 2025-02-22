# Function to get a list of RangedSummarizedExperiment and SummarizedExperiment objects based on provided paths
# Removes subjects with missing values for inputted list of features
# Get experiment out of list using syntax: experimentList[[index]]
getExperimentsList <- function(paths, featuresToEnsure=c("Age"), removeHPVPositive=FALSE, removeVitalStatusUnreported=TRUE) {
  
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
    if (removeVitalStatusUnreported){ experiment <- experiment[, experiment$vital_status != "Not Reported"]}
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
  else if (layerName == "RNAseq_Normal") (revisedLayerName <- "RNAseq_NORMAL_filtered")
  else { revisedLayerName = layerName }
  paths <- getPathsByLayerAndCancer()
  cancerPaths <- getLayerPathsForSpecificCancer(paths, cancerName)
  layerPathsList <- list()
  layerPathsList[[layerName]] <- cancerPaths[[revisedLayerName]]
  experiment <- getExperimentsList(layerPathsList)
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

graphAgeModel <- function(tests, title, caption, biasCorrection=FALSE) {
  if (biasCorrection) {
    predictions <- tests$predicted_corrected_age
  } else {
    predictions <- tests$predicted_age
  }
  
  DF <- data.frame(
    sampleID = rownames(tests),
    true_age = tests$Age,
    pred_age = predictions
  )
  
  g1 <- ggplot2::ggplot(DF, aes(x = true_age, y = pred_age)) +
    geom_point(color = "blue") +
    geom_smooth(method = "lm", se = FALSE, aes(color = "Regression Line")) +  # Adding linear model fit line with legend
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", aes(color = "Identity Line")) +  # Identity line with legend
    scale_color_manual(values = c("Regression Line" = "red")) +  # Define custom colors
    labs(
      title = paste("True Age vs Predicted Age (", title, ")", sep=""),
      x = "True Age",
      y = "Predicted Age",
      color = "Legend"
    ) +
    labs(caption = caption) +
    ylim(10, 100) +  # Fix y-axis range
    theme(
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0)
    )
  
  print(g1)
}

reduceExperimentBySubjects <- function(exp1, exp2){
  return(exp1[, exp1$submitter_id %in% exp2$submitter_id])
}

reduceExperimentByFeatures <- function(exp1, exp2, assayName1, assayName2=NULL){
  if (is.null(assayName2)) { assayName2 <- assayName1 }
  assay1 <- getTransposedFrame(SummarizedExperiment::assay(exp1, assayName1))
  assay2 <- getTransposedFrame(SummarizedExperiment::assay(exp2, assayName2))
  assay1Filtered <- assay1[, names(assay1) %in% names(assay2)]
  commonFeatures <- names(assay1)[names(assay1) %in% names(assay2)]
  return(exp1[commonFeatures, ])
}

getFold <- function(data, fold){
  meta_trn = data[-fold, ]
  meta_tst = data[fold, ]
  return(list(train = meta_trn, test = meta_tst))
}

getNormalSize <- function(cancerName){
  normExp <- getExperimentByLayerAndCancer("RNAseq_Normal", cancerName)
  tumorExp <- getExperimentByLayerAndCancer("RNAseq", cancerName)
  experimentsList <- list(RNAseq = reduceExperimentBySubjects(tumorExp[[1]], normExp[[1]]))
  needsImputation <- list(RNAseq = TRUE, RNAseq_Normal = FALSE)
  assayNames <- list(RNAseq = "DESeq2_log", RNAseq_Normal = "DESeq2_log")
  imputedExperiments <- imputeMissingValues(experimentsList, assayNames, needsImputation, cancerName)
  return(imputedExperiments)
}

getTransposedFrame <- function(matrix){
  return(as.data.frame(t(matrix)))
}

getCancerLayersPrefix <- function(cancerType, layers){
  return(paste0(cancerType, "_", paste(layers, collapse = "_")))
}