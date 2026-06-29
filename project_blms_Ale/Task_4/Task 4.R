# ------------------------------------------------------------
# Setup
# This block sets the seed, defines the selected model,
# creates the Task 4 output folders, and loads the needed data.
# ------------------------------------------------------------

# Clear the current R environment.
rm(list = ls())

# Set the seed to make the results reproducible.
set.seed(123)

# Selected mixture model.
# This must match the H used in Task 3.
H_selected <- 3

# Number of Fourier harmonics used for the hour effect.
K_fourier <- 2

# Check that the selected values are valid.
stopifnot(H_selected %in% 1:5)
stopifnot(H_selected >= 2)
stopifnot(K_fourier >= 1)

# Main project folder.
base_dir <- "/home/blackmamba/project_blms/project_blms_Ale"

# Input folders.
data_dir <- file.path(base_dir, "data")
task_3_dir <- file.path(base_dir, "Task_3")
task_3_tables_dir <- file.path(task_3_dir, "results", "tables")

# Task 4 output folders.
task_4_dir <- file.path(base_dir, "Task_4")
results_dir <- file.path(task_4_dir, "results")
fits_dir <- file.path(results_dir, "fits")
tables_dir <- file.path(results_dir, "tables")
plots_dir <- file.path(results_dir, "plots")

# Create the output folders if they do not already exist.
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Path to the input dataset.
data_path <- file.path(data_dir, "hourly_load_factor.csv")

# Path to the MAP component assignments from Task 3.
map_assignments_path <- file.path(
  task_3_tables_dir,
  paste0("map_component_assignments_H", H_selected, ".rds")
)

# Load the original dataset.
data <- read.csv(data_path)

# Load MAP component assignments from Task 3.
map_assignments <- readRDS(map_assignments_path)

# Check that the hour variable exists in the original dataset.
if (!("hour" %in% names(data))) {
  stop("The dataset must contain a variable named 'hour'.")
}

# Check that the MAP component variable exists.
if (!("map_component" %in% names(map_assignments))) {
  stop("The MAP assignment file must contain a variable named 'map_component'.")
}

# Check that the two datasets have the same number of observations.
stopifnot(nrow(data) == nrow(map_assignments))

# Add hour to the MAP assignment table.
map_assignments$hour <- as.numeric(data$hour)

# Check that hour is correctly defined.
if (any(is.na(map_assignments$hour))) {
  stop("The hour variable contains missing or non-numeric values.")
}

# Convert the response variable to a factor.
map_assignments$map_component <- factor(
  map_assignments$map_component,
  levels = seq_len(H_selected)
)

# Check that all components are represented in the data.
component_counts <- table(map_assignments$map_component)

if (any(component_counts == 0)) {
  stop("At least one MAP component has zero observations. Multinomial regression cannot estimate all classes.")
}

# ------------------------------------------------------------
# Classification dataset with Fourier predictors
# This block builds the dataset used for the multinomial
# logistic regression.
# ------------------------------------------------------------

# Start from the MAP assignments from Task 3.
classification_data <- data.frame(
  observation = map_assignments$observation,
  map_component = map_assignments$map_component,
  hour = map_assignments$hour
)

# Add Fourier-transformed hour predictors.
for (k in seq_len(K_fourier)) {
  
  classification_data[[paste0("sin_", k)]] <- sin(
    k * 2 * pi * classification_data$hour / 24
  )
  
  classification_data[[paste0("cos_", k)]] <- cos(
    k * 2 * pi * classification_data$hour / 24
  )
}

# Save the classification dataset as RDS.
saveRDS(
  classification_data,
  file = file.path(
    tables_dir,
    paste0("classification_data_H", H_selected, "_K", K_fourier, ".rds")
  )
)

# Save the classification dataset as CSV.
write.csv(
  classification_data,
  file = file.path(
    tables_dir,
    paste0("classification_data_H", H_selected, "_K", K_fourier, ".csv")
  ),
  row.names = FALSE
)

# Print a short summary.
cat("\nClassification dataset completed.\n")
cat("Selected model: H =", H_selected, "\n")
cat("Fourier harmonics: K =", K_fourier, "\n")
cat("Number of observations:", nrow(classification_data), "\n")
cat("Variables included:\n")
print(names(classification_data))

# ------------------------------------------------------------
# Fit multinomial logistic regression
# This block fits the classification model using the
# Fourier-transformed hour predictors.
# ------------------------------------------------------------

# Load the package needed for multinomial logistic regression.
library(nnet)

# Build the model formula automatically from the Fourier predictors.
fourier_terms <- c()

for (k in seq_len(K_fourier)) {
  fourier_terms <- c(
    fourier_terms,
    paste0("sin_", k),
    paste0("cos_", k)
  )
}

multinom_formula <- as.formula(
  paste(
    "map_component ~",
    paste(fourier_terms, collapse = " + ")
  )
)

# Fit the multinomial logistic regression model.
multinom_fit <- multinom(
  formula = multinom_formula,
  data = classification_data,
  trace = FALSE
)

# Save the fitted model.
saveRDS(
  multinom_fit,
  file = file.path(
    fits_dir,
    paste0("multinomial_logistic_fit_H", H_selected, "_K", K_fourier, ".rds")
  )
)

# Save the model summary.
multinom_summary <- summary(multinom_fit)

saveRDS(
  multinom_summary,
  file = file.path(
    fits_dir,
    paste0("multinomial_logistic_summary_H", H_selected, "_K", K_fourier, ".rds")
  )
)

# Predict fitted class probabilities for the observed data.
fitted_probabilities_raw <- predict(
  multinom_fit,
  type = "probs"
)

# Convert fitted probabilities to a data frame.
if (H_selected == 2 && is.null(dim(fitted_probabilities_raw))) {
  
  fitted_probabilities <- data.frame(
    prob_component_1 = 1 - fitted_probabilities_raw,
    prob_component_2 = fitted_probabilities_raw
  )
  
} else {
  
  fitted_probabilities <- as.data.frame(fitted_probabilities_raw)
  colnames(fitted_probabilities) <- paste0("prob_component_", seq_len(H_selected))
}

# Add observation, hour, and observed MAP component.
fitted_probabilities <- data.frame(
  observation = classification_data$observation,
  hour = classification_data$hour,
  map_component = classification_data$map_component,
  fitted_probabilities
)

# Save fitted probabilities.
saveRDS(
  fitted_probabilities,
  file = file.path(
    tables_dir,
    paste0("fitted_component_probabilities_H", H_selected, "_K", K_fourier, ".rds")
  )
)

write.csv(
  fitted_probabilities,
  file = file.path(
    tables_dir,
    paste0("fitted_component_probabilities_H", H_selected, "_K", K_fourier, ".csv")
  ),
  row.names = FALSE
)

# Print a short summary.
cat("\nMultinomial logistic regression completed.\n")
cat("Selected model: H =", H_selected, "\n")
cat("Fourier harmonics: K =", K_fourier, "\n")
cat("Model formula:\n")
print(multinom_formula)
cat("\nModel summary:\n")
print(multinom_summary)

# ------------------------------------------------------------
# Predicted component probabilities by hour
# This block predicts the probability of each component
# on an hourly grid from 0 to 23.
# ------------------------------------------------------------

# Create the hourly prediction grid.
hour_grid <- data.frame(
  hour = 0:23
)

# Add Fourier-transformed hour predictors to the prediction grid.
for (k in seq_len(K_fourier)) {
  
  hour_grid[[paste0("sin_", k)]] <- sin(
    k * 2 * pi * hour_grid$hour / 24
  )
  
  hour_grid[[paste0("cos_", k)]] <- cos(
    k * 2 * pi * hour_grid$hour / 24
  )
}

# Predict component probabilities for each hour.
predicted_probabilities_raw <- predict(
  multinom_fit,
  newdata = hour_grid,
  type = "probs"
)

# Convert predicted probabilities to a data frame.
if (H_selected == 2 && is.null(dim(predicted_probabilities_raw))) {
  
  predicted_probabilities <- data.frame(
    prob_component_1 = 1 - predicted_probabilities_raw,
    prob_component_2 = predicted_probabilities_raw
  )
  
} else {
  
  predicted_probabilities <- as.data.frame(predicted_probabilities_raw)
  colnames(predicted_probabilities) <- paste0("prob_component_", seq_len(H_selected))
}

# Add the hour variable.
predicted_probabilities <- data.frame(
  hour = hour_grid$hour,
  predicted_probabilities
)

# Save predicted probabilities as RDS.
saveRDS(
  predicted_probabilities,
  file = file.path(
    tables_dir,
    paste0("predicted_component_probabilities_by_hour_H", H_selected, "_K", K_fourier, ".rds")
  )
)

# Save predicted probabilities as CSV.
write.csv(
  predicted_probabilities,
  file = file.path(
    tables_dir,
    paste0("predicted_component_probabilities_by_hour_H", H_selected, "_K", K_fourier, ".csv")
  ),
  row.names = FALSE
)

# Print a short summary.
cat("\nPredicted component probabilities by hour completed.\n")
cat("Selected model: H =", H_selected, "\n")
cat("Fourier harmonics: K =", K_fourier, "\n")
cat("Prediction grid: hours 0 to 23\n")
cat("\nPredicted probabilities:\n")
print(round(predicted_probabilities, 4))


# ------------------------------------------------------------
# Plot predicted component probabilities by hour
# This block saves one plot for each component and one combined
# plot with all predicted component probabilities.
# ------------------------------------------------------------

# Fixed colors for all possible components.
component_colors_all <- c(
  "1" = "red",
  "2" = "orange",
  "3" = "forestgreen",
  "4" = "purple",
  "5" = "blue"
)

# Colors used for the selected model.
component_colors <- component_colors_all[as.character(seq_len(H_selected))]

# Common plot size.
plot_width <- 1800
plot_height <- 1100
plot_res <- 160

# Common plot limits.
hour_xlim <- range(predicted_probabilities$hour)
probability_ylim <- c(0, 1)

# Common graphical settings.
set_plot_style <- function() {
  par(
    mar = c(4.8, 5.0, 3.4, 1.2),
    mgp = c(2.8, 0.8, 0),
    cex.main = 1.25,
    cex.lab = 1.15,
    cex.axis = 1.05,
    cex = 1.05
  )
}

# ------------------------------------------------------------
# Individual predicted probability plots
# This block saves one plot for each component.
# ------------------------------------------------------------

for (h in seq_len(H_selected)) {
  
  prob_col <- paste0("prob_component_", h)
  color_h <- component_colors_all[as.character(h)]
  
  png(
    filename = file.path(
      plots_dir,
      paste0(
        "predicted_probability_by_hour_H",
        H_selected,
        "_K",
        K_fourier,
        "_component_",
        h,
        ".png"
      )
    ),
    width = plot_width,
    height = plot_height,
    res = plot_res
  )
  
  set_plot_style()
  
  plot(
    predicted_probabilities$hour,
    predicted_probabilities[[prob_col]],
    type = "n",
    xlim = hour_xlim,
    ylim = probability_ylim,
    xlab = "Hour of day",
    ylab = "Predicted probability",
    main = paste0("Predicted probability of component ", h)
  )
  
  grid()
  
  lines(
    predicted_probabilities$hour,
    predicted_probabilities[[prob_col]],
    type = "b",
    pch = 19,
    lwd = 2.6,
    col = color_h
  )
  
  legend(
    "topright",
    legend = paste0("Component ", h),
    col = color_h,
    lwd = 2.6,
    pch = 19,
    bty = "n"
  )
  
  dev.off()
}

# ------------------------------------------------------------
# Combined predicted probability plot
# This block saves one plot with all components together.
# ------------------------------------------------------------

png(
  filename = file.path(
    plots_dir,
    paste0(
      "predicted_probabilities_by_hour_H",
      H_selected,
      "_K",
      K_fourier,
      "_all_components.png"
    )
  ),
  width = plot_width,
  height = plot_height,
  res = plot_res
)

set_plot_style()

plot(
  NA,
  NA,
  xlim = hour_xlim,
  ylim = probability_ylim,
  xlab = "Hour of day",
  ylab = "Predicted probability",
  main = paste0(
    "Predicted MAP component probabilities, H = ",
    H_selected,
    ", K = ",
    K_fourier
  )
)

grid()

for (h in seq_len(H_selected)) {
  
  prob_col <- paste0("prob_component_", h)
  color_h <- component_colors_all[as.character(h)]
  
  lines(
    predicted_probabilities$hour,
    predicted_probabilities[[prob_col]],
    type = "b",
    pch = 19,
    lwd = 2.6,
    col = color_h
  )
}

legend(
  "topright",
  legend = paste0("Component ", seq_len(H_selected)),
  col = component_colors,
  lwd = 2.6,
  pch = 19,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# Comparison data from Task 3
# This block loads the empirical component fractions from Task 3
# so that the Task 3 and Task 4 plots can be compared.
# ------------------------------------------------------------

# Path to the empirical component fractions from Task 3.
empirical_fractions_path <- file.path(
  task_3_tables_dir,
  paste0("component_fractions_by_hour_H", H_selected, ".rds")
)

# Load empirical component fractions.
empirical_fractions <- readRDS(empirical_fractions_path)

# Save a copy inside Task 4 results for convenience.
saveRDS(
  empirical_fractions,
  file = file.path(
    tables_dir,
    paste0("task3_empirical_component_fractions_H", H_selected, ".rds")
  )
)

write.csv(
  empirical_fractions,
  file = file.path(
    tables_dir,
    paste0("task3_empirical_component_fractions_H", H_selected, ".csv")
  ),
  row.names = FALSE
)

# Print a short summary.
cat("\nPredicted probability plots completed.\n")
cat("Selected model: H =", H_selected, "\n")
cat("Fourier harmonics: K =", K_fourier, "\n")
cat("Saved individual and combined plots in Task_4/results/plots.\n")
cat("Loaded Task 3 empirical fractions for comparison.\n")


# ------------------------------------------------------------
# Task 3 vs Task 4 comparison plots
# This block compares empirical component fractions from Task 3
# with predicted component probabilities from Task 4.
# ------------------------------------------------------------

# Check that the required objects exist.
stopifnot(exists("empirical_fractions"))
stopifnot(exists("predicted_probabilities"))

# Rename Task 3 empirical columns for clarity.
empirical_comparison <- empirical_fractions[, c(
  "hour",
  "map_component",
  "fraction"
)]

names(empirical_comparison) <- c(
  "hour",
  "component",
  "empirical_fraction"
)

# Convert component column to numeric.
empirical_comparison$component <- as.numeric(
  as.character(empirical_comparison$component)
)

# Build Task 4 predicted comparison table in long format.
predicted_comparison <- data.frame()

for (h in seq_len(H_selected)) {
  
  prob_col <- paste0("prob_component_", h)
  
  predicted_h <- data.frame(
    hour = predicted_probabilities$hour,
    component = h,
    predicted_probability = predicted_probabilities[[prob_col]]
  )
  
  predicted_comparison <- rbind(
    predicted_comparison,
    predicted_h
  )
}

# Merge empirical fractions and predicted probabilities.
comparison_data <- merge(
  empirical_comparison,
  predicted_comparison,
  by = c("hour", "component")
)

# Order comparison data.
comparison_data <- comparison_data[
  order(comparison_data$hour, comparison_data$component),
]

# Save comparison data as RDS.
saveRDS(
  comparison_data,
  file = file.path(
    tables_dir,
    paste0("task3_vs_task4_comparison_H", H_selected, "_K", K_fourier, ".rds")
  )
)

# Save comparison data as CSV.
write.csv(
  comparison_data,
  file = file.path(
    tables_dir,
    paste0("task3_vs_task4_comparison_H", H_selected, "_K", K_fourier, ".csv")
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Individual comparison plots
# This block saves one Task 3 vs Task 4 plot for each component.
# ------------------------------------------------------------

for (h in seq_len(H_selected)) {
  
  component_h <- comparison_data[
    comparison_data$component == h,
  ]
  
  color_h <- component_colors_all[as.character(h)]
  
  png(
    filename = file.path(
      plots_dir,
      paste0(
        "task3_vs_task4_H",
        H_selected,
        "_K",
        K_fourier,
        "_component_",
        h,
        ".png"
      )
    ),
    width = plot_width,
    height = plot_height,
    res = plot_res
  )
  
  set_plot_style()
  
  plot(
    component_h$hour,
    component_h$empirical_fraction,
    type = "n",
    xlim = hour_xlim,
    ylim = probability_ylim,
    xlab = "Hour of day",
    ylab = "Probability / empirical fraction",
    main = paste0("Task 3 vs Task 4, component ", h)
  )
  
  grid()
  
  # Task 4 predicted probability as a line.
  lines(
    component_h$hour,
    component_h$predicted_probability,
    type = "l",
    lwd = 3.0,
    col = color_h
  )
  
  # Task 3 empirical fraction as points.
  points(
    component_h$hour,
    component_h$empirical_fraction,
    pch = 19,
    cex = 1.2,
    col = color_h
  )
  
  legend(
    "topright",
    legend = c(
      paste0("Task 4 predicted probability, component ", h),
      paste0("Task 3 empirical fraction, component ", h)
    ),
    col = c(color_h, color_h),
    lwd = c(3.0, NA),
    pch = c(NA, 19),
    bty = "n"
  )
  
  dev.off()
}

# ------------------------------------------------------------
# Combined comparison plot
# This block saves one Task 3 vs Task 4 plot with all components.
# ------------------------------------------------------------

png(
  filename = file.path(
    plots_dir,
    paste0(
      "task3_vs_task4_H",
      H_selected,
      "_K",
      K_fourier,
      "_all_components.png"
    )
  ),
  width = plot_width,
  height = plot_height,
  res = plot_res
)

set_plot_style()

plot(
  NA,
  NA,
  xlim = hour_xlim,
  ylim = probability_ylim,
  xlab = "Hour of day",
  ylab = "Probability / empirical fraction",
  main = paste0(
    "Task 3 empirical fractions vs Task 4 predicted probabilities, H = ",
    H_selected,
    ", K = ",
    K_fourier
  )
)

grid()

for (h in seq_len(H_selected)) {
  
  component_h <- comparison_data[
    comparison_data$component == h,
  ]
  
  color_h <- component_colors_all[as.character(h)]
  
  # Task 4 predicted probability as a line.
  lines(
    component_h$hour,
    component_h$predicted_probability,
    type = "l",
    lwd = 3.0,
    col = color_h
  )
  
  # Task 3 empirical fraction as points.
  points(
    component_h$hour,
    component_h$empirical_fraction,
    pch = 19,
    cex = 1.1,
    col = color_h
  )
}

legend(
  "topright",
  legend = c(
    paste0("Component ", seq_len(H_selected), " - Task 4 line"),
    paste0("Component ", seq_len(H_selected), " - Task 3 points")
  ),
  col = c(component_colors, component_colors),
  lwd = c(rep(3.0, H_selected), rep(NA, H_selected)),
  pch = c(rep(NA, H_selected), rep(19, H_selected)),
  bty = "n"
)

dev.off()

# Print a short summary.
cat("\nTask 3 vs Task 4 comparison completed.\n")
cat("Selected model: H =", H_selected, "\n")
cat("Fourier harmonics: K =", K_fourier, "\n")
cat("Saved comparison table and plots in Task_4/results.\n")