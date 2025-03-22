library(readr)

# Correct file path including the filename
hnsc <- read_csv("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/HNSC_elastic_net_omics_combinations.csv")

setwd("/restricted/projectnb/agedisease/projects/pancancer_aging_clocks/scripts/GitCore/results/")

#Permutation Test Variable Importance - Vip (you have your elastic net, how much does your model decrease if it changes the weights) 
hnsc$test_name

hnsc$delta_age_pval_non_interaction
hnsc$chronological_pval_non_interaction
hnsc$delta_age_pval_interaction
hnsc$chronological_pval_interaction
hnsc$chronological_pval_baseline

####################################################SAMPLE WORKS BEST##########################################################
# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Select relevant columns
hnsc_filtered <- hnsc %>%
  select(layer_combination, delta_age_pval_non_interaction, chronological_pval_non_interaction,
         delta_age_pval_interaction, chronological_pval_interaction, 
         chronological_pval_baseline)

# Initialize empty lists to store data for plotting
p_values_list <- list()
labels_list <- list()
colors_list <- list()

# Define colors for categories
color_map <- c("Non-Interaction" = "blue", "Interaction" = "red", "Baseline" = "green")

# Loop through each layer_combination and collect p-values
for (i in 1:nrow(hnsc_filtered)) {
  # Extract p-values
  p_values <- as.numeric(hnsc_filtered[i, c("delta_age_pval_non_interaction", 
                                            "chronological_pval_non_interaction",
                                            "delta_age_pval_interaction", 
                                            "chronological_pval_interaction", 
                                            "chronological_pval_baseline")])
  
  # Create labels and colors
  pval_labels <- rep(hnsc_filtered$layer_combination[i], 5)  # Use layer_combination as y-axis label
  colors <- c(color_map["Non-Interaction"], color_map["Non-Interaction"], 
              color_map["Interaction"], color_map["Interaction"], 
              color_map["Baseline"])
  
  # Store in lists
  p_values_list[[i]] <- p_values
  labels_list[[i]] <- pval_labels
  colors_list[[i]] <- colors
}

# Combine all p-values, labels, and colors
all_p_values <- unlist(p_values_list)
all_labels <- unlist(labels_list)
all_colors <- unlist(colors_list)

par(mar = c(5, 12, 4, 2))  # Bottom, Left, Top, Right

# Plot horizontal bar plot for all test names
bar_positions <- barplot(all_p_values,
                         names.arg = all_labels,
                         horiz = TRUE,
                         col = all_colors,
                         main = "P-Values Across All Tests",
                         xlab = "P-Value",
                         las = 1,    # Rotate y-axis labels for readability
                         cex.names = 0.8)  # Reduce label size if needed

# Add a vertical threshold line at p-value = 0.05
abline(v = 0.05, lty = 2, col = "black", lwd = 2)  # Dotted line

# Add a legend
legend("topright", legend = c("Non-Interaction", "Interaction", "Baseline"),
       fill = c("blue", "red", "green"))

############################################################################################################################################################


#################################################### CHRONOLOGICAL ONLY ########################################################################################################
# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Select relevant columns
hnsc_filtered <- hnsc %>%
  select(layer_combination, chronological_pval_non_interaction, chronological_pval_interaction, chronological_pval_baseline)

# Apply element-wise inverse only to numeric columns
hnsc_filtered <- hnsc_filtered %>%
  mutate(across(where(is.numeric), ~ -log(.)))

# Initialize empty lists to store data for plotting
p_values_list <- list()
labels_list <- list()
colors_list <- list()

# Define colors for categories
color_map <- c("Non-Interaction" = "blue", "Interaction" = "red", "Baseline" = "green")

# Loop through each layer_combination and collect p-values
for (i in 1:nrow(hnsc_filtered)) {
  # Extract p-values
  p_values <- as.numeric(hnsc_filtered[i, c("chronological_pval_non_interaction",
                                            "chronological_pval_interaction", 
                                            "chronological_pval_baseline")])
  
  # Create labels and colors
  pval_labels <- rep(hnsc_filtered$layer_combination[i], 3)  # Use layer_combination as y-axis label
  colors <- c(color_map["Non-Interaction"], 
              color_map["Interaction"], 
              color_map["Baseline"])
  
  # Store in lists
  p_values_list[[i]] <- p_values
  labels_list[[i]] <- pval_labels
  colors_list[[i]] <- colors
}

# Combine all p-values, labels, and colors
all_p_values <- unlist(p_values_list)
all_labels <- unlist(labels_list)
all_colors <- unlist(colors_list)

# Create labels: One per layer_combination
unique_labels <- hnsc_filtered$layer_combination
label_positions <- seq(3, length(all_p_values), by = 3)  # Position labels in the middle of each group

par(mar = c(5, 12, 4, 2))  # Bottom, Left, Top, Right

bar_positions <- barplot(all_p_values,
                         horiz = TRUE,
                         col = all_colors,
                         main = "-log Chronological Significance Values",
                         xlab = "-log(P-Value)",
                         yaxt = "n",  # Suppress default y-axis labels
                         las = 1)

# Add custom y-axis labels (one per layer_combination)
axis(2, at = bar_positions[label_positions], labels = unique_labels, las = 1, cex.axis = 0.8)

# Add a vertical threshold line at p-value = 0.05
abline(v = 2.995732, lty = 2, col = "black", lwd = 2)  # Dotted line

# Add a legend
legend("topright", legend = c("Baseline","Interaction", "Non-Interaction"),
       fill = c("green", "red", "blue"))

#################################################### DELTA_AGE ONLY ########################################################################################################
# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Select relevant columns
hnsc_filtered <- hnsc %>%
  select(layer_combination, delta_age_pval_non_interaction, delta_age_pval_interaction)

# Apply element-wise inverse only to numeric columns
hnsc_filtered <- hnsc_filtered %>%
  mutate(across(where(is.numeric), ~ -log(.)))

# Initialize empty lists to store data for plotting
p_values_list <- list()
labels_list <- list()
colors_list <- list()

# Define colors for categories
color_map <- c("Non-Interaction" = "blue", "Interaction" = "red")

# Loop through each layer_combination and collect p-values
for (i in 1:nrow(hnsc_filtered)) {
  # Extract p-values
  p_values <- as.numeric(hnsc_filtered[i, c("delta_age_pval_non_interaction",
                                            "delta_age_pval_interaction"
                                            )])
  
  # Create labels and colors
  pval_labels <- rep(hnsc_filtered$layer_combination[i], 2)  # Use layer_combination as y-axis label
  colors <- c(color_map["Non-Interaction"], 
              color_map["Interaction"])
  
  # Store in lists
  p_values_list[[i]] <- p_values
  labels_list[[i]] <- pval_labels
  colors_list[[i]] <- colors
}

# Combine all p-values, labels, and colors
all_p_values <- unlist(p_values_list)
all_labels <- unlist(labels_list)
all_colors <- unlist(colors_list)

# Create labels: One per layer_combination
unique_labels <- hnsc_filtered$layer_combination
label_positions <- seq(0, length(all_p_values), by = 2)  # Position labels in the middle of each group

par(mar = c(5, 12, 4, 2))  # Bottom, Left, Top, Right

bar_positions <- barplot(all_p_values,
                         horiz = TRUE,
                         col = all_colors,
                         main = "-log delta age Significance Values",
                         xlab = "-log(P-Value)",
                         yaxt = "n",  # Suppress default y-axis labels
                         las = 1)

# Add custom y-axis labels (one per layer_combination)
axis(2, at = bar_positions[label_positions], labels = unique_labels, las = 1, cex.axis = 0.8)

# Add a vertical threshold line at p-value = 0.05
abline(v = 2.995732, lty = 2, col = "black", lwd = 2)  # Dotted line

# Add a legend
legend("topright", legend = c("Interaction","Non-Interaction"),
       fill = c("red","blue"))

#################################################### DELTA_AGE ONLY ########################################################################################################

#One Plot with the chronological only
#One Plot with the deltas only










test <- hnsc_filtered^(-1)

locator(1)

# Example: Custom coordinates (adjust as needed)
legend(x = 0.03060285, y = 48.68329, legend = c("Non-Interaction", "Interaction", "Baseline"),
       fill = c("blue", "red", "green"), bty = "n")  # bty = "n" removes box around legend








# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)

# Filter data for test = "RPPA"
hnsc_filtered <- hnsc %>%
  filter(layer_combination == "RPPA") %>%
  select(test_name, delta_age_pval_non_interaction, chronological_pval_non_interaction,
         delta_age_pval_interaction, chronological_pval_interaction, 
         lik_ratio_test_pval_baseline) %>%
  rename(chronological_pval_baseline = lik_ratio_test_pval_baseline)

# Create a vector of p-values
p_values <- as.numeric(hnsc_filtered[1, c("delta_age_pval_non_interaction", 
                                          "chronological_pval_non_interaction",
                                          "delta_age_pval_interaction", 
                                          "chronological_pval_interaction", 
                                          "chronological_pval_baseline")])

# Create labels for the bars
pval_labels <- c("Delta Age", 
                 "Chronological",
                 "Delta Age", 
                 "Chronological", 
                 "Chronological")

# Define colors
colors <- c("blue", "blue", "red", "red", "green")

# Plot horizontal bar plot
barplot(p_values,
        names.arg = pval_labels,
        horiz = TRUE,
        col = colors,
        main = "P-Values for RPPA",
        xlab = "P-Value",
        las = 1)  # Rotate y-axis labels for readability

# Add a legend
legend("topright", legend = c("Non-Interaction", "Interaction", "Baseline"),
       fill = c("blue", "red", "green"))















The following are the values I want to access, make a plot where in blue the non-interaction terms are shown, in red the interaction, and in green the baseline
Show these as 5 bars in the y axis, as horizontal bars, there should be 5 of these bars of the specified types per test name.

#Include bias corrected, and include at the end, a table with for each model, the coefficients and their significance

# Load required library
library(ggplot2)
library(tidyr)
library(dplyr)

# Select relevant columns for plotting
hnsc_plot <- hnsc %>%
  select(test_name, delta_age_pval_non_interaction, chronological_pval_non_interaction,
         delta_age_pval_interaction, chronological_pval_interaction, 
         lik_ratio_test_pval_baseline) %>%
  rename(chronological_pval_baseline = lik_ratio_test_pval_baseline)

# Reshape data to long format for ggplot
hnsc_long <- hnsc_plot %>%
  pivot_longer(
    cols = -test_name,
    names_to = "pval_type",
    values_to = "p_value"
  ) %>%
  mutate(
    category = case_when(
      grepl("non_interaction", pval_type) ~ "Non-Interaction",
      grepl("interaction", pval_type) ~ "Interaction",
      grepl("baseline", pval_type) ~ "Baseline"
    )
  )

# Define color mapping
color_map <- c("Non-Interaction" = "blue", "Interaction" = "red", "Baseline" = "green")

# Create the horizontal bar plot
ggplot(hnsc_long, aes(x = p_value, y = test_name, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = color_map) +
  labs(title = "P-Values by Test Name",
       x = "P-Value",
       y = "Test Name",
       fill = "Category") +
  theme_minimal() +
  coord_flip()






# Reshape data to long format
hnsc_long <- hnsc %>%
  pivot_longer(
    cols = -test_name,
    names_to = "pval_type",
    values_to = "p_value"
  )

# Define color mapping
hnsc_long <- hnsc_long %>%
  mutate(
    category = case_when(
      grepl("non_interaction", pval_type) ~ "Non-Interaction",
      grepl("interaction", pval_type) ~ "Interaction",
      grepl("baseline", pval_type) ~ "Baseline"
    )
  )

# Set colors
color_map <- c("Non-Interaction" = "blue", "Interaction" = "red", "Baseline" = "green")

# Plot
ggplot(hnsc_long, aes(x = p_value, y = test_name, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = color_map) +
  labs(title = "P-Values by Test Name",
       x = "P-Value",
       y = "Test Name",
       fill = "Category") +
  theme_minimal() +
  coord_flip()
