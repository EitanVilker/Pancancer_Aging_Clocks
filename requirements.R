# requirements.R
packages <- c(
  "Biobase", "SummarizedExperiment", "ggplot2", "reactable", "caret",
  "glmnet", "dplyr", "impute", "survival", "survminer", "vip", "tidyverse",
  "randomForest", "e1071", "class", "caTools", "VGAM", "SEtools",
  "TimeSeriesExperiment"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, install_if_missing)
