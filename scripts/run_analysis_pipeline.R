run_analysis_pipeline <- function(fullExperiment, prediction_df, prediction_object=NULL) {
  library(dplyr)
  library(survival)
  library(SummarizedExperiment)
  library(ggplot2)
  library(reactable)
  library(caret)
  library(glmnet)
  library(dplyr)
  library(impute)
  library(survival)
  library(survminer)
  library(vip)
  library(tidyverse)
  
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
    DFsurv$surv <- survival::Surv(time = DFsurv$fu, event = as.numeric(factor(DFsurv$vitals)) - 1)
    return(DFsurv)
  }
  
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
  formula_vector_non_interaction <- c("chronological + delta_age", "race")
  formula_vector_interaction <- c("chronological * delta_age", "race")
  if (is.null(prediction_object) || "gendermale" %in% colnames(prediction_object)){
    formula_vector_non_interaction <- c(formula_vector_non_interaction, "gender.y")
    formula_vector_interaction <- c(formula_vector_interaction, "gender.y")
  }
  covariates_formula_non_interaction <- as.formula(paste("surv ~ ", paste(formula_vector_non_interaction, collapse= "+")))
  covariates_formula_interaction <- as.formula(paste("surv ~ ", paste(formula_vector_interaction, collapse= "+")))

  # Run Non-Interaction CoxPH Model
  test1_non_interaction <- coxph(
    covariates_formula_non_interaction,
    data = experiment_prediction_object_meta
  )
  
  # Run Interaction CoxPH Model
  test1_interaction <- coxph(
    covariates_formula_interaction,
    data = experiment_prediction_object_meta
  )
  
  # Summarize Results
  non_interaction_summary <- summary(test1_non_interaction)
  interaction_summary <- summary(test1_interaction)
  
  # Return Results
  return(list(
    non_interaction_model = test1_non_interaction,
    non_interaction_summary = non_interaction_summary,
    interaction_model = test1_interaction,
    interaction_summary = interaction_summary
  ))
}