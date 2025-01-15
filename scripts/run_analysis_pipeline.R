run_analysis_pipeline <- function(HNSC, three_layers_elnet) {
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
  make_meta_df <- function(HNSC) {
    metadata_table <- data.frame(submitter_id = rownames(HNSC@colData))
    for (name in names(HNSC@colData@listData)) {
      metadata_table[[name]] <- HNSC@colData@listData[[name]]
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
  metadata_table <- make_meta_df(HNSC)
  rownames(metadata_table) <- metadata_table$submitter_id
  
  # Prepare Predicted Ages Data
  predicted_ages <- three_layers_elnet %>%
    dplyr::select(submitter_id, predicted_age)
  rownames(predicted_ages) <- predicted_ages$submitter_id
  
  # Filter Metadata for Matching Submitter IDs
  filtered_metadata <- metadata_table %>%
    dplyr::semi_join(predicted_ages, by = "submitter_id")
  rownames(filtered_metadata) <- filtered_metadata$submitter_id
  
  # Create Survival Dataset
  hnsc_three_layers_elnet <- create_DFsurv(predicted_ages, filtered_metadata)
  
  hnsc_three_layers_elnet_meta <- hnsc_three_layers_elnet %>%
    inner_join(filtered_metadata, by = "submitter_id")
  
  # Run Non-Interaction CoxPH Model
  test1_non_interaction <- coxph(
    surv ~ chronological + delta_age + gender.y + race,
    data = hnsc_three_layers_elnet_meta
  )
  
  # Run Interaction CoxPH Model
  test1_interaction <- coxph(
    surv ~ chronological * delta_age + gender.y + race,
    data = hnsc_three_layers_elnet_meta
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