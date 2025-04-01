---
title: "Analysis"
output: html_notebook
---



``` r
library(survival)
library(survminer)
source("run_analysis_pipeline.R")
```

# Runs CoxPH given experiment and predictions and output object containing each CoxPH run

``` r
getSurvivalStats <- function(fullExperiment, predictionObject, applyBiasCorrection=FALSE, covariates_to_include=c("gender", "race")){
  if (applyBiasCorrection) { predictedAges <- predictionObject$predicted_corrected_age }
  else { predictedAges <- predictionObject$predicted_age }
  df <- data.frame(predicted_age=predictedAges, submitter_id=predictionObject$submitter_id)
  return(run_analysis_pipeline(fullExperiment, df, prediction_object=predictionObject, covariates_to_include=covariates_to_include))
}
```

# Gets likelihood ratio between two models

``` r
getTrueLikelihoodRatio <- function(baselineModel, modifiedModel){
  
  baselineLikelihood <- logLik(baselineModel)
  modifiedLikelihood <- logLik(modifiedModel)
  newLogRatio <- -2 * (baselineLikelihood - modifiedLikelihood)
  degreesOfFreedom <- length(coef(modifiedModel)) - length(coef(baselineModel))
  pValue <- 1 - pchisq(newLogRatio, degreesOfFreedom)
  return(list(logratio=newLogRatio, pval=pValue))
}
```

### Save age model features and summary stats as well as survival analysis. Important to set additionalDescriptor to unique name if you don't want to overwrite previous outputs

``` r
saveClockOutput <- function(cancerType, layers, modelOutput, dir, additionalDescriptor, crossValidating=FALSE, applyBiasCorrection=FALSE){
  
  # Modify output path
  cancerLayersPrefix <- getCancerLayersPrefix(cancerType, layers)
  primaryDir <- paste0(dir, cancerType)
  dir.create(primaryDir, recursive = TRUE, showWarnings = FALSE)
  filePath <- paste0(primaryDir, "/", cancerLayersPrefix)

  # Get coefficients
  coefficientFrame <- list()
  if (!crossValidating){
    coefficientFrame[[cancerLayersPrefix]] <- getCoefficientDF(modelOutput$ModelBuilding$model$finalModel)
  } else{
     coefficientFrame[[cancerLayersPrefix]] <- modelOutput$ModelBuilding$combinedWeights
  }
  
  # Get stats
  stats <- getSurvivalStats(modelOutput$experimentList[[1]], modelOutput$ModelBuilding$predicted, applyBiasCorrection=applyBiasCorrection, covariates_to_include=c("race"))
  statsTable <- getFormattedStatsTable(stats)
  summaryTable <- getSummaryTable(cancerType, layers, modelOutput, stats)
  print(statsTable)
  
  # Save outputs
  coefficientPath <- paste0(filePath, "_coef", additionalDescriptor, ".rds")
  saveRDS(coefficientFrame, file = coefficientPath)
  write.csv(summaryTable, paste0(filePath, "_summary", additionalDescriptor, ".csv"), row.names = FALSE)
  write.csv(statsTable, paste0(filePath, "_stats", additionalDescriptor, ".csv"), row.names = TRUE)
  return(list(primaryDir=primaryDir, coefficientPath=coefficientPath))
}
```

### Perform pathway analysis

``` r
performEnrichment <- function(cancerType, layers, additionalDescriptor, clockOutput){
  
  # Initialize files and folders
  cancerLayersPrefix <- getCancerLayersPrefix(cancerType, layers)
  outputDir=paste0(clockOutput$primaryDir,"/Enrichment/", additionalDescriptor)
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  formattedWeightFileName <- paste0(outputDir, "/", cancerLayersPrefix, "_weights_named.csv") 

  # Run pathway analysis
  process_gene_weights(cancerType, weights_file=clockOutput$coefficientPath, output_dir=outputDir)
  run_pathway_analysis(cancerType, input_file=formattedWeightFileName, output_dir=outputDir)
}
```

### Get coefficient df

``` r
getCoefficientDF <- function(finalModel){

  # Extract lambda values and coefficients
  lambda_values <- finalModel$lambda
  closest_lambda_index <- which.min(abs(lambda_values - finalModel$lambdaOpt))
  beta_matrix <- finalModel$beta
  coef_names <- rownames(beta_matrix)
  beta_column <- beta_matrix[, closest_lambda_index]

  coef_df <- data.frame(
    Feature = coef_names,
    Weight = as.numeric(beta_column)
  )
  rownames(coef_df) <- coef_df$Feature
  return(coef_df)
}
```

### Map RNAseq gene IDs to names

``` r
getGeneMap <- function(modelOutput){
  features <- modelOutput$ModelBuilding$model$coefnames
  geneIDs <- c()
  geneNames <- c()
  RNAseqExperiment <- modelOutput$experimentList$RNAseq
  for (name in features){
    if (substr(name, 1, 4) == "ENSG"){
      geneIDs <- c(geneIDs, name)
      geneIndex <- which(rowRanges(RNAseqExperiment)$gene_id %in% name)
      geneName <- rowRanges(RNAseqExperiment)$gene_name[geneIndex]
      geneNames <- c(geneNames, geneName)
    }
  }
  
  geneMap <- data.frame(Ensembl_ID=geneIDs, Gene_Name=geneNames)
  return(geneMap)
}
```

### Find coefficients of individual features used in final model

``` r
## For RNAseq data
# To find coefficients: coefficients <- coef(modelPredictionsPair$model$finalModel, modelPredictionsPair$model$bestTune$lambda)
# To find the index of a coefficient: index <- which(coefficients@x %in% COEFFICIENT_OF_CHOICE)
# To find the gene ID of an index: geneID <- coefficients@Dimnames[[1]][index]
# To find the index of the gene in the full experiment: fullIndex <- which(rowRanges(experimentsList[[1]])$gene_id %in% geneID) 
# To find the gene name: featureName <- rowRanges(experimentsList[[1]])$gene_name[fullIndex]
# RNASeqElNet <- read.csv("./ModelBuildingOutput/RNAseqElNet.csv")

organizeOutputFeatures <- function(modelOutput, RNAExperiment=NULL){

  coefficients <- coef(modelOutput$model$finalModel, modelOutput$model$bestTune$lambda)
  # coefficients2 <- coef(modelOutput$model$finalModel, unlist(modelOutput$model$bestTune))
  featureFrame <- data.frame(features=coefficients@Dimnames[[1]][coefficients@i + 1], coefficients=coefficients@x)
  featureNames <- c()
  origin <- c()
  for (name in featureFrame$features){
    featureIdx <- which(colnames(modelOutput$predicted) %in% name)
    if (!is.null(RNAExperiment) && substr(name, 1, 4) == "ENSG"){
      origin <- c(origin, "RNAseq")
      geneIndex <- which(rowRanges(experimentsList[[1]])$gene_id %in% name)
      geneName <- rowRanges(RNAExperiment)$gene_name[geneIndex]
      featureNames <- append(featureNames, rowRanges(modelOutput$experimentsList$RNAseq)$gene_name[geneIndex])
    }
    else{
      featureNames <- append(featureNames, name)
      if(substr(name, 2, 4) == "hsa"){origin <- c(origin, "miRNA")}
      else if(length(featureIdx) > 0 && featureIdx < 9){origin <- c(origin, "covariate")}
      else if(name == "(Intercept)"){origin <- c(origin, "Intercept")}
      else{origin <- c(origin, "RPPA")}
    }
  }
  
  featureFrame$names <- featureNames
  featureFrame$origin <- origin
  featureFrame$absCoefficients <- abs(featureFrame$coefficients)
  featureFrame <- featureFrame[order(featureFrame$absCoefficients, decreasing=TRUE),]
  return(featureFrame)
}

# featureFrame <- organizeOutputFeatures(modelPredictionsPair, RNAExperiment=experimentsList[["RNAseq"]])
# 
# outFile <- "./ModelBuildingOutput/ElasticNetRNAseqHoldoutsFeatures.csv"
# write.csv(featureFrame, outFile, row.names = FALSE)
```


``` r
intializeSummaryTable <- function(){
  cancer_summary_table <- data.frame(
    test_name = character(),
    cancer_type = character(),
    model = character(),
    age_association_adjusted_p_cutoff = numeric(),
    layer_combination = character(),
    
    Rsquared = numeric(),
    RMSE = numeric(),
    alpha = numeric(),
    lambda = numeric(),
    non_zero_weight_feature_count = numeric(),
    zero_weight_feature_count = numeric(),
    
    interaction_against_baseline_likelihood_ratio = numeric(),
    interaction_against_baseline_pval = numeric(),
    interaction_term_coeff_interaction = numeric(),
    interaction_term_exp_coeff_interaction = numeric(),
    interaction_term_pval_interaction = numeric(),
    delta_age_coeff_interaction	= numeric(),
    delta_age_exp_coeff_interaction = numeric(),
    delta_age_pval_interaction = numeric(),
    chronological_coef_interaction = numeric(),
    chronological_exp_coef_interaction = numeric(),
    chronological_pval_interaction = numeric(),
    Concordance_interaction = numeric(),
    se_concordance_interaction = numeric(),
    lik_ratio_test_interaction = numeric(),
    lik_ratio_test_df_interaction = numeric(),
    lik_ratio_test_pval_interaction = numeric(),
    Wald_test_interaction = numeric(),
    Wald_test_df_interaction = numeric(),
    Wald_test_pval_interaction = numeric(),
    logrank_score_test_interaction = numeric(),
    logrank_score_test_df_interaction = numeric(),
    logrank_score_test_pval_interaction = numeric(),
    
    non_interaction_against_baseline_likelihood_ratio = numeric(),
    non_interaction_against_baseline_pval = numeric(),
    delta_age_coeff_non_interaction	= numeric(),
    delta_age_exp_coeff_non_interaction = numeric(),
    delta_age_pval_non_interaction = numeric(),
    chronological_coef_non_interaction = numeric(),
    chronological_exp_coef_non_interaction = numeric(),
    chronological_pval_non_interaction = numeric(),
    Concordance_non_interaction = numeric(),
    se_concordance_non_interaction = numeric(),
    lik_ratio_test_non_interaction = numeric(),
    lik_ratio_test_df_non_interaction = numeric(),
    lik_ratio_test_pval_non_interaction = numeric(),
    Wald_test_non_interaction = numeric(),
    Wald_test_df_non_interaction = numeric(),
    Wald_test_pval_non_interaction = numeric(),
    logrank_score_test_non_interaction = numeric(),
    logrank_score_test_df_non_interaction = numeric(),
    logrank_score_test_pval_non_interaction = numeric(),
    
    chronological_coef_baseline = numeric(),
    chronological_exp_coef_baseline = numeric(),
    chronological_pval_baseline = numeric(),
    Concordance_baseline = numeric(),
    se_concordance_baseline = numeric(),
    lik_ratio_test_baseline = numeric(),
    lik_ratio_test_df_baseline = numeric(),
    lik_ratio_test_pval_baseline = numeric(),
    Wald_test_baseline = numeric(),
    Wald_test_df_baseline = numeric(),
    Wald_test_pval_baseline = numeric(),
    logrank_score_test_baseline = numeric(),
    logrank_score_test_df_baseline = numeric(),
    logrank_score_test_pval_baseline = numeric()
  )
  return(cancer_summary_table)
}
```


``` r
getSummaryTable <- function(comb, modelOutput, stats, params=NULL){
  interactionLikelihood <- getTrueLikelihoodRatio(stats$baseline_model, stats$interaction_model)
  nonInteractionLikelihood <- getTrueLikelihoodRatio(stats$baseline_model, stats$non_interaction_model)
  if ("Fold1" %in% names(modelOutput$ModelBuilding$model) || "Fold01" %in% names(modelOutput$ModelBuilding$model)){
    alpha <- modelOutput$ModelBuilding$model[[1]]$model$bestTune$alpha
    lambda <- modelOutput$ModelBuilding$model[[1]]$model$bestTune$lambda
  }
  else{
    alpha = modelOutput$ModelBuilding$model$bestTune$alpha
    lambda = modelOutput$ModelBuilding$model$bestTune$lambda
  }
  non_zero_weights <- sum(modelOutput$ModelBuilding$weights$Weight != 0, na.rm = TRUE)
  zero_weights <- sum(modelOutput$ModelBuilding$weights$Weight == 0, na.rm = TRUE)
  if (!is.null(params)){
    test_name <- params$test_name
    model <- params$model_type
    age_association_adjusted_p_cutoff <- params$significance_cutoff
    cancer_type <- params$cancer_type
  }
  else{
    test_name <- "None"
    model <- "None"
    age_association_adjusted_p_cutoff <- "None"
    cancer_type <- "None"
  }

  summaryTable <- data.frame(
                    test_name = test_name,
                    cancer_type = cancer_type,
                    model = model,
                    age_association_adjusted_p_cutoff = age_association_adjusted_p_cutoff,
                    layer_combination = paste(comb, collapse = "_"),
                    
                    Rsquared = modelOutput$ModelBuilding$model$R2,
                    RMSE = modelOutput$ModelBuilding$model$RMSE,
                    alpha = alpha,
                    lambda = lambda,
                    non_zero_weight_feature_count = non_zero_weights,
                    zero_weight_feature_count = zero_weights,
                    
                    interaction_against_baseline_likelihood_ratio = interactionLikelihood$logratio,
                    interaction_against_baseline_pval = interactionLikelihood$pval,
                    interaction_term_coeff_interaction = stats$interaction_summary$coefficients["chronological:delta_age", "coef"],
                    interaction_term_exp_coeff_interaction = stats$interaction_summary$coefficients["chronological:delta_age", "exp(coef)"],
                    interaction_term_pval_interaction = stats$interaction_summary$coefficients["chronological:delta_age", "Pr(>|z|)"],
                    delta_age_coeff_interaction = stats$interaction_summary$coefficients["delta_age", "coef"],
                    delta_age_exp_coeff_interaction = stats$interaction_summary$coefficients["delta_age", "exp(coef)"],
                    delta_age_pval_interaction = stats$interaction_summary$coefficients["delta_age", "Pr(>|z|)"],
                    chronological_coef_interaction = stats$interaction_summary$coefficients["chronological", "coef"],
                    chronological_exp_coef_interaction = stats$interaction_summary$coefficients["chronological", "exp(coef)"],
                    chronological_pval_interaction = stats$interaction_summary$coefficients["chronological", "Pr(>|z|)"],
                    Concordance_interaction = stats$interaction_summary$concordance[1],
                    se_concordance_interaction = stats$interaction_summary$concordance[2],
                    lik_ratio_test_interaction = stats$interaction_summary$logtest["test"],
                    lik_ratio_test_df_interaction = stats$interaction_summary$logtest["df"],
                    lik_ratio_test_pval_interaction = stats$interaction_summary$logtest["pvalue"],
                    Wald_test_interaction = stats$interaction_summary$waldtest["test"],
                    Wald_test_df_interaction = stats$interaction_summary$waldtest["df"],
                    Wald_test_pval_interaction = stats$interaction_summary$waldtest["pvalue"],
                    logrank_score_test_interaction = stats$interaction_summary$sctest["test"],
                    logrank_score_test_df_interaction = stats$interaction_summary$sctest["df"],
                    logrank_score_test_pval_interaction = stats$interaction_summary$sctest["pvalue"],
                    
                    non_interaction_against_baseline_likelihood_ratio = nonInteractionLikelihood$logratio,
                    non_interaction_against_baseline_pval = nonInteractionLikelihood$pval,
                    delta_age_coeff_non_interaction	= stats$non_interaction_summary$coefficients["delta_age", "coef"],
                    delta_age_exp_coeff_non_interaction = stats$non_interaction_summary$coefficients["delta_age", "exp(coef)"],
                    delta_age_pval_non_interaction = stats$non_interaction_summary$coefficients["delta_age", "Pr(>|z|)"],
                    chronological_coef_non_interaction = stats$non_interaction_summary$coefficients["chronological", "coef"],
                    chronological_exp_coef_non_interaction = stats$non_interaction_summary$coefficients["chronological", "exp(coef)"],
                    chronological_pval_non_interaction =  stats$non_interaction_summary$coefficients["chronological", "Pr(>|z|)"],
                    Concordance_non_interaction = stats$non_interaction_summary$concordance[1],
                    se_concordance_non_interaction = stats$non_interaction_summary$concordance[2],
                    lik_ratio_test_non_interaction = stats$non_interaction_summary$logtest["test"],
                    lik_ratio_test_df_non_interaction = stats$non_interaction_summary$logtest["df"],
                    lik_ratio_test_pval_non_interaction = stats$non_interaction_summary$logtest["pvalue"],
                    Wald_test_non_interaction = stats$non_interaction_summary$waldtest["test"],
                    Wald_test_df_non_interaction = stats$non_interaction_summary$waldtest["df"],
                    Wald_test_pval_non_interaction = stats$non_interaction_summary$waldtest["pvalue"],
                    logrank_score_test_non_interaction = stats$non_interaction_summary$sctest["test"],
                    logrank_score_test_df_non_interaction = stats$non_interaction_summary$sctest["df"],
                    logrank_score_test_pval_non_interaction = stats$non_interaction_summary$sctest["pvalue"],
                    
                    chronological_coef_baseline = stats$baseline_summary$coefficients["chronological", "coef"],
                    chronological_exp_coef_baseline = stats$baseline_summary$coefficients["chronological", "exp(coef)"],
                    chronological_pval_baseline = stats$baseline_summary$coefficients["chronological", "Pr(>|z|)"],
                    Concordance_baseline = stats$baseline_summary$concordance[1],
                    se_concordance_baseline = stats$baseline_summary$concordance[2],
                    lik_ratio_test_baseline = stats$baseline_summary$logtest["test"],
                    lik_ratio_test_df_baseline = stats$baseline_summary$logtest["df"],
                    lik_ratio_test_pval_baseline = stats$baseline_summary$logtest["pvalue"],
                    Wald_test_baseline = stats$baseline_summary$waldtest["test"],
                    Wald_test_df_baseline = stats$baseline_summary$waldtest["df"],
                    Wald_test_pval_baseline = stats$baseline_summary$waldtest["pvalue"],
                    logrank_score_test_baseline = stats$baseline_summary$sctest["test"],
                    logrank_score_test_df_baseline = stats$baseline_summary$sctest["df"],
                    logrank_score_test_pval_baseline = stats$baseline_summary$sctest["pvalue"]
  )

  return(summaryTable)
}
```


``` r
getFormattedStatsTable <- function(stats){
  
  cols <- c("interaction_term_pval", "delta_age_pval", "chronological_age_pval", "concordance", "logtest", "logtest pvalue", "waldtest", "waldtest pvalue", "sctest", "sctest pvalue")
  non_interaction <- c(
    "N/A",
    stats$non_interaction_summary$coefficients["delta_age", "Pr(>|z|)"],
    stats$non_interaction_summary$coefficients["chronological", "Pr(>|z|)"],
    stats$non_interaction_summary$concordance[1],
    stats$non_interaction_summary$logtest["test"],
    stats$non_interaction_summary$logtest["pvalue"],
    stats$non_interaction_summary$waldtest["test"],
    stats$non_interaction_summary$waldtest["pvalue"],
    stats$non_interaction_summary$sctest["test"],
    stats$non_interaction_summary$sctest["pvalue"]
  )
  
  interaction <- c(
    stats$interaction_summary$coefficients["chronological:delta_age", "Pr(>|z|)"],
    stats$interaction_summary$coefficients["delta_age", "Pr(>|z|)"],
    stats$interaction_summary$coefficients["chronological", "Pr(>|z|)"],
    stats$interaction_summary$concordance[1],
    stats$interaction_summary$logtest["test"],
    stats$interaction_summary$logtest["pvalue"],
    stats$interaction_summary$waldtest["test"],
    stats$interaction_summary$waldtest["pvalue"],
    stats$interaction_summary$sctest["test"],
    stats$interaction_summary$sctest["pvalue"]
  )
  
  baseline <- c(
    "N/A",
    "N/A",
    stats$baseline_summary$coefficients["chronological", "Pr(>|z|)"],
    stats$baseline_summary$concordance[1],
    stats$baseline_summary$logtest["test"],
    stats$baseline_summary$logtest["pvalue"],
    stats$baseline_summary$waldtest["test"],
    stats$baseline_summary$waldtest["pvalue"],
    stats$baseline_summary$sctest["test"],
    stats$baseline_summary$sctest["pvalue"]
  )
  statsFrame <- t(data.frame(interaction=interaction, non_interaction=non_interaction, baseline=baseline))
  colnames(statsFrame) <- cols
  return(statsFrame)
}
```

# Function to produce a summary table of layers and cancers for their ML model info and survival stats. Excludes cancers with insufficient sample sizes.

``` r
writeUltimateSummaryTable <- function(outputDir="/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/Eitan/RidgeCV/", suffix="_ridge_omics_summary.csv"){
  cancer_types <- getEligibleCancers("RNAseq")
  cancer_types <- cancer_types[!cancer_types == "OV"] # Too little in methylation
  summary_file_paths <- paste0(outputDir, 
                        cancer_types, 
                        suffix)
  

  # Read and combine all CSVs
  summaries_list <- lapply(summary_file_paths, function(file) {
    if (file.exists(file)) {
      read.csv(file, row.names = NULL)  # Adjust row.names if needed
    } else {
      warning(paste("File not found:", file))
      return(NULL)
    }
  })
  
  combined_tables <- do.call(rbind, summaries_list)
  if (!("cancer_type" %in% colnames(combined_tables))){
    rep_cancers <- rep(cancer_types, each=4)
    combined_tables <- cbind(cancer_type=rep_cancers, combined_tables)
  }
  write.csv(combined_tables, paste0(outputDir, "UltimateSummaryTable", suffix))
}
```


``` r
library(vip)
featureImportance <- function(result){
  # Calculate permutation feature importance
  vip_data <- vi(result$elnet, method = "permute", target = "Age", train = meta_age_associated, metric = "MAE",
                 pred_wrapper = function(object, newdata) predict(object, newdata))
  
  # Plot the variable importance
  vip(vip_data, geom = "point", num_features = 20, horiz = TRUE) +
    labs(title = "Permutation Feature Importance for Elastic Net Model (Full DataSet)")
}
```
