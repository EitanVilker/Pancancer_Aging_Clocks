# library(SummarizedExperiment)
# library(ggplot2)
# library(reactable)
# library(caret)
# library(glmnet)
# library(dplyr)
# library(impute)
# library(survival)
# library(survminer)
# library(vip)
# library(tidyverse)

# Step 1: Create Metadata Table
make_meta_df <- function(fullExperiment) {
  metadata_table <- data.frame(submitter_id = rownames(fullExperiment@colData))
  for (name in names(fullExperiment@colData@listData)) {
    metadata_table[[name]] <- fullExperiment@colData@listData[[name]]
  }
  metadata_table[is.na(metadata_table)] <- 0
  return(metadata_table)
}

# Step 2: Prepare Survival Data
create_DFsurv <- function(predicted_ages, metadata_table, delta_age_thresh = 0) {
  #if (!identical(rownames(predicted_ages), rownames(metadata_table))) {
  #  stop("Row names of predicted_ages and metadata_table must match.")
  #}
  metadata_table$predicted_age <- predicted_ages$predicted_age
  metadata_table$delta_age <- metadata_table$predicted_age - metadata_table$Age
  DFsurv <- metadata_table |>
    dplyr::select(submitter_id, Age, predicted_age, vital_status, delta_age, Survival_Time, gender, HPV.status) |>
    dplyr::rename(
      chronological = Age,
      predicted = predicted_age,
      vitals = vital_status,
      fu = Survival_Time
    ) |>
    dplyr::mutate(meta_age = factor(
      dplyr::case_when(
        delta_age > delta_age_thresh ~ "older",
        delta_age < delta_age_thresh ~ "younger",
        TRUE ~ "same"
      ), levels = c("younger", "older", "same")
    )) |>
    dplyr::mutate(age_strata = cut(
      chronological,
      breaks = c(0, 50, 70, 90, +Inf),
      labels = c("0-50", "50-70", "70-90", "90+")
    ))
  DFsurv$surv <- survival::Surv(time = DFsurv$fu / 365, event = as.numeric(factor(DFsurv$vitals)) - 1)
  DFsurv$chronological <- scale(DFsurv$chronological)
  DFsurv$chronological <- DFsurv$chronological - min(DFsurv$chronological) + 0.01
  DFsurv$chronological_duplicate <- DFsurv$chronological + rnorm(length(DFsurv$chronological), mean=0, sd=0.011)
  DFsurv$chronological_duplicate2 <- DFsurv$chronological + rnorm(length(DFsurv$chronological), mean=0, sd=0.012)
  DFsurv$delta_age <- scale(DFsurv$delta_age)
  DFsurv$delta_age <- DFsurv$delta_age - min(DFsurv$delta_age) + 0.01
  
  return(DFsurv)
}

# Step 3: Run the CoxPH analysis
run_analysis_pipeline <- function(fullExperiment, prediction_df, useCovariates=TRUE, useGender=TRUE, degreesFreedomNonLinear=8, useComplicatedNonlinear=TRUE) {
  
  # Prepare Metadata Table
  metadata_table <- make_meta_df(fullExperiment)
  rownames(metadata_table) <- metadata_table$submitter_id
  
  # Prepare Predicted Ages Data
  predicted_ages <- prediction_df %>%
    dplyr::select(submitter_id, predicted_age)
  rownames(predicted_ages) <- predicted_ages$submitter_id
  
  # Filter Metadata for Matching Submitter IDs
  filtered_metadata <- metadata_table %>%
    dplyr::semi_join(predicted_ages, by = "submitter_id")
  rownames(filtered_metadata) <- filtered_metadata$submitter_id
  
  # Create Survival Dataset
  experiment_prediction_object <- create_DFsurv(predicted_ages, filtered_metadata)

  experiment_prediction_object_meta <- experiment_prediction_object %>%
    inner_join(filtered_metadata, by = "submitter_id")

  # Set up formulas, accounting for covariates appropriate to the cancer
  formulaVectorList = list()
  nonlinear_age_term <- paste0("pspline(chronological, ", degreesFreedomNonLinear, ")")
  
  # Add core deta_age models
  formulaVectorList$non_interaction_linear <- c("chronological + delta_age")
  formulaVectorList$interaction_linear <- c("chronological * delta_age")
  formulaVectorList$non_interaction_nonlinear <- c(paste0(nonlinear_age_term, "+ delta_age"))
  formulaVectorList$interaction_additive_nonlinear <- c(paste0(nonlinear_age_term, "+ chronological * delta_age"))
  formulaVectorList$interaction_linear_duplicate <- c("chronological * chronological_duplicate + chronological * delta_age")
  
  # Add baseline models without delta_age terms
  formulaVectorList$baseline <- c("chronological")
  formulaVectorList$baseline_nonlinear <- c(nonlinear_age_term)
  formulaVectorList$baseline_interaction <- c("chronological * chronological_duplicate")
  formulaVectorList$baseline_interaction2 <- c("chronological * chronological_duplicate + chronological * chronological_duplicate2")
  formulaVectorList$baseline_interaction_nonlinear_additive <- c(paste0(nonlinear_age_term, "+ chronological * chronological_duplicate"))
  
  # Add models that can fail with high degrees of freedom
  if (useComplicatedNonlinear){
    formulaVectorList$interaction_multiplicative_nonlinear <- c(paste0(nonlinear_age_term, "* delta_age"))
    formulaVectorList$interaction_nonlinear_duplicate <- c(paste0(nonlinear_age_term, "* chronological_duplicate + ", nonlinear_age_term, "* delta_age"))
    formulaVectorList$baseline_interaction_nonlinear_multiplicative <- c(paste0(nonlinear_age_term, "* chronological_duplicate"))
    formulaVectorList$baseline_interaction_nonlinear_multiplicative2 <- c(paste0(nonlinear_age_term, "* chronological_duplicate + ", nonlinear_age_term, "* chronological_duplicate2"))
  }
  
  # Add appropriate covariates if using
  if (useCovariates){
    if (useGender){
      if (nlevels(as.factor(experiment_prediction_object_meta$gender.y)) > 1){
        for (formula in names(formulaVectorList)){
          formulaVectorList[[formula]] <- c(formulaVectorList[[formula]], "gender.y")
        }
      }
    }
    if (nlevels(as.factor(experiment_prediction_object_meta$subtype_selected)) > 1 && TRUE){
      for (formula in names(formulaVectorList)){
        formulaVectorList[[formula]] <- c(formulaVectorList[[formula]], "subtype_selected")
      }
    }
    if (nlevels(as.factor(experiment_prediction_object_meta$treatments_radiation_treatment_or_therapy)) > 1 && TRUE){
      for (formula in names(formulaVectorList)){
        formulaVectorList[[formula]] <- c(formulaVectorList[[formula]], "treatments_radiation_treatment_or_therapy")
      }
    }
  }
  
  # Run the CoxPH model for each test
  testList <- list()
  for (formula in names(formulaVectorList)){
    formulaObject <- as.formula(paste("surv ~ ", paste(formulaVectorList[[formula]], collapse= "+")))
    testList[[formula]] <- coxph(
      formulaObject,
      data = experiment_prediction_object_meta,
      control = coxph.control(iter.max = 50, outer.max = 20)
    )
  }
  
  # Output the prediction object and each model and summary
  outputList <- list(survDF=experiment_prediction_object_meta)
  for (test in names(testList)){
    model <- testList[[test]]
    outputList[[paste0(test, "_model")]] <- model
    outputList[[paste0(test, "_summary")]] <- summary(model)
  }
  return(outputList)
}