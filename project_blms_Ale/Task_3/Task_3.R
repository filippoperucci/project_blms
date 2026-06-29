# ------------------------------------------------------------
# Setup
# This block sets the seed, defines the selected model,
# creates the Task 3 output folders, and loads the data.
# ------------------------------------------------------------

# Clear the current R environment.
rm(list = ls())

# Set the seed to make the results reproducible.
set.seed(123)

# Selected mixture model.
# Choose one value manually between 1 and 5.
H_selected <- 3

# Check that the selected value is valid.
stopifnot(H_selected %in% 1:5)

# Main project folder.
base_dir <- "/home/blackmamba/project_blms/project_blms_Ale"

# Input folders.
data_dir <- file.path(base_dir, "data")
task_1_dir <- file.path(base_dir, "Task_1")
task_1_fits_dir <- file.path(task_1_dir, "results", "fits")

# Task 3 output folders.
task_3_dir <- file.path(base_dir, "Task_3")
results_dir <- file.path(task_3_dir, "results")
tables_dir <- file.path(results_dir, "tables")
plots_dir <- file.path(results_dir, "plots")

# Create the output folders if they do not already exist.
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Path to the input dataset.
data_path <- file.path(data_dir, "hourly_load_factor.csv")

# Load the dataset.
data <- read.csv(data_path)

# Extract the response variable.
x <- data$load_factor
n <- length(x)

# Check that all observations are inside the Beta support.
stopifnot(all(x > 0), all(x < 1))

# Check that the hour variable exists.
if (!("hour" %in% names(data))) {
  stop("The dataset must contain a variable named 'hour'.")
}

# Select the fitted model file from Task 1.
if (H_selected == 1) {
  fit_path <- file.path(task_1_fits_dir, "fit_beta_single_H1.rds")
} else {
  fit_path <- file.path(task_1_fits_dir, paste0("fit_beta_mixture_H", H_selected, ".rds"))
}

# Load the fitted model.
fit_selected <- readRDS(fit_path)

# Extract posterior samples.
samps_selected <- fit_selected$samples

# ------------------------------------------------------------
# Definition of MAP component assignment function
# This function computes posterior component probabilities
# and assigns each observation to its MAP component.
# Components are ordered by increasing component mean.
# ------------------------------------------------------------

compute_map_components <- function(x, samps, H) {
  
  # Convert MCMC samples to a matrix.
  mat <- as.matrix(samps)
  
  # Number of observations.
  n <- length(x)
  
  # Handle the single Beta model.
  if (H == 1) {
    
    posterior_component_prob <- matrix(1, nrow = n, ncol = 1)
    colnames(posterior_component_prob) <- "component_1"
    
    out <- data.frame(
      observation = seq_len(n),
      x = x,
      map_component = rep(1, n),
      posterior_component_prob
    )
    
    return(out)
  }
  
  # Number of posterior draws.
  S <- nrow(mat)
  
  # Create matrices for mixture parameters.
  alpha <- matrix(NA, nrow = S, ncol = H)
  beta  <- matrix(NA, nrow = S, ncol = H)
  p     <- matrix(NA, nrow = S, ncol = H)
  
  # Extract alpha, beta, and mixture weights for each component.
  for (h in 1:H) {
    alpha[, h] <- mat[, paste0("alpha[", h, "]")]
    beta[, h]  <- mat[, paste0("beta[", h, "]")]
    p[, h]     <- mat[, paste0("p[", h, "]")]
  }
  
  # Order components within each posterior draw by increasing component mean.
  alpha_ord <- alpha
  beta_ord  <- beta
  p_ord     <- p
  
  component_mean <- alpha / (alpha + beta)
  
  for (s in seq_len(S)) {
    ord <- order(component_mean[s, ])
    alpha_ord[s, ] <- alpha[s, ord]
    beta_ord[s, ]  <- beta[s, ord]
    p_ord[s, ]     <- p[s, ord]
  }
  
  # Store posterior probabilities P(z_i = h | x_i, data).
  posterior_component_prob <- matrix(
    NA,
    nrow = n,
    ncol = H
  )
  
  colnames(posterior_component_prob) <- paste0("component_", seq_len(H))
  
  # Compute posterior component probabilities for each observation.
  for (i in seq_len(n)) {
    
    # Store unnormalized probabilities for each posterior draw.
    prob_s <- matrix(NA, nrow = S, ncol = H)
    
    for (h in seq_len(H)) {
      prob_s[, h] <- p_ord[, h] * dbeta(
        x[i],
        shape1 = alpha_ord[, h],
        shape2 = beta_ord[, h]
      )
    }
    
    # Normalize probabilities across components for each posterior draw.
    prob_s <- prob_s / rowSums(prob_s)
    
    # Average probabilities over posterior draws.
    posterior_component_prob[i, ] <- colMeans(prob_s)
  }
  
  # Assign each observation to its MAP component.
  map_component <- max.col(
    posterior_component_prob,
    ties.method = "first"
  )
  
  # Store results.
  out <- data.frame(
    observation = seq_len(n),
    x = x,
    map_component = map_component,
    posterior_component_prob
  )
  
  out
}

# ------------------------------------------------------------
# MAP component assignment
# This block computes posterior component probabilities
# and assigns each observation to its MAP component.
# ------------------------------------------------------------

# Compute posterior component probabilities and MAP assignments.
map_assignments <- compute_map_components(
  x = x,
  samps = samps_selected,
  H = H_selected
)

# Save MAP assignments.
saveRDS(
  map_assignments,
  file = file.path(
    tables_dir,
    paste0("map_component_assignments_H", H_selected, ".rds")
  )
)

# Save MAP assignments also as CSV.
write.csv(
  map_assignments,
  file = file.path(
    tables_dir,
    paste0("map_component_assignments_H", H_selected, ".csv")
  ),
  row.names = FALSE
)

# Print a short summary.
cat("\nMAP component assignment completed.\n")
cat("Selected model: H =", H_selected, "\n")
cat("Number of observations:", nrow(map_assignments), "\n")

cat("\nNumber of observations assigned to each component:\n")
print(table(map_assignments$map_component))

cat("\nEmpirical fraction assigned to each component:\n")
print(round(prop.table(table(map_assignments$map_component)), 4))

# ------------------------------------------------------------
# Empirical component fractions by hour
# This block computes the fraction of observations assigned
# to each MAP component for every hour.
# ------------------------------------------------------------

# Add the hour variable to the MAP assignments.
map_assignments$hour <- data$hour

# Create all hour-component combinations.
hour_component_grid <- expand.grid(
  hour = sort(unique(data$hour)),
  map_component = seq_len(H_selected)
)

# Count observations by hour and MAP component.
component_counts <- as.data.frame(
  table(
    hour = map_assignments$hour,
    map_component = map_assignments$map_component
  )
)

# Convert columns to numeric.
component_counts$hour <- as.numeric(as.character(component_counts$hour))
component_counts$map_component <- as.numeric(as.character(component_counts$map_component))

# Merge with the full grid to include missing combinations.
component_counts <- merge(
  hour_component_grid,
  component_counts,
  by = c("hour", "map_component"),
  all.x = TRUE
)

# Replace missing counts with zero.
component_counts$Freq[is.na(component_counts$Freq)] <- 0

# Compute total observations for each hour.
hour_totals <- aggregate(
  Freq ~ hour,
  data = component_counts,
  FUN = sum
)

names(hour_totals)[2] <- "hour_total"

# Add hourly totals.
component_fractions <- merge(
  component_counts,
  hour_totals,
  by = "hour"
)

# Compute empirical fractions.
component_fractions$fraction <- component_fractions$Freq / component_fractions$hour_total

# Order the table.
component_fractions <- component_fractions[
  order(component_fractions$hour, component_fractions$map_component),
]

# Save the table as RDS.
saveRDS(
  component_fractions,
  file = file.path(
    tables_dir,
    paste0("component_fractions_by_hour_H", H_selected, ".rds")
  )
)

# Save the table as CSV.
write.csv(
  component_fractions,
  file = file.path(
    tables_dir,
    paste0("component_fractions_by_hour_H", H_selected, ".csv")
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot settings
# This block defines common plot settings for all component plots.
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
hour_xlim <- range(component_fractions$hour)
fraction_ylim <- c(0, 1)

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
# Individual component fraction plots
# This block saves one plot for each component.
# ------------------------------------------------------------

for (h in seq_len(H_selected)) {
  
  component_h <- component_fractions[
    component_fractions$map_component == h,
  ]
  
  color_h <- component_colors_all[as.character(h)]
  
  png(
    filename = file.path(
      plots_dir,
      paste0("component_fraction_by_hour_H", H_selected, "_component_", h, ".png")
    ),
    width = plot_width,
    height = plot_height,
    res = plot_res
  )
  
  set_plot_style()
  
  plot(
    component_h$hour,
    component_h$fraction,
    type = "n",
    xlim = hour_xlim,
    ylim = fraction_ylim,
    xlab = "Hour of day",
    ylab = "Fraction of observations",
    main = paste0("Hourly fraction assigned to component ", h)
  )
  
  grid()
  
  lines(
    component_h$hour,
    component_h$fraction,
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
# Combined component fraction plot
# This block saves one plot with all components together.
# ------------------------------------------------------------

png(
  filename = file.path(
    plots_dir,
    paste0("component_fractions_by_hour_H", H_selected, "_all_components.png")
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
  ylim = fraction_ylim,
  xlab = "Hour of day",
  ylab = "Fraction of observations",
  main = paste0("Hourly MAP component fractions, H = ", H_selected)
)

grid()

for (h in seq_len(H_selected)) {
  
  component_h <- component_fractions[
    component_fractions$map_component == h,
  ]
  
  color_h <- component_colors_all[as.character(h)]
  
  lines(
    component_h$hour,
    component_h$fraction,
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

# Print a short summary.
cat("\nComponent fractions by hour completed.\n")
cat("Selected model: H =", H_selected, "\n")
cat("Saved table and plots in Task_3/results.\n")