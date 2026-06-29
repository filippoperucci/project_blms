# ------------------------------------------------------------
# Setup
# This block sets the seed, defines the needed folders,
# creates the Task 2 output folders, and loads the data.
# ------------------------------------------------------------

# Clear the current R environment.
rm(list = ls())

# Set the seed to make the results reproducible.
set.seed(123)

# Models to compare.
# Choose values manually between 1 and 5.
# Example: c(1, 3, 5) compares only H = 1, H = 3, and H = 5.
H_to_compare <- c(1, 2, 3, 4, 5)

# Check that the selected values are valid.
stopifnot(all(H_to_compare %in% 1:5))

# Remove possible duplicated values.
H_to_compare <- unique(H_to_compare)

# Main project folder.
base_dir <- "/home/blackmamba/project_blms/project_blms_Ale"

# Input folders.
data_dir <- file.path(base_dir, "data")
task_1_dir <- file.path(base_dir, "Task_1")
task_1_fits_dir <- file.path(task_1_dir, "results", "fits")

# Task 2 output folders.
task_2_dir <- file.path(base_dir, "Task_2")
results_dir <- file.path(task_2_dir, "results")
fits_dir <- file.path(results_dir, "fits")
plots_dir <- file.path(results_dir, "plots")

# Create the output folders if they do not already exist.
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
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

# Common grid used for all predictive densities.
grid_x <- seq(0.001, 0.999, length.out = 300)

# Fixed colors for all possible models.
model_colors_all <- c(
  "1" = "red",
  "2" = "orange",
  "3" = "forestgreen",
  "4" = "purple",
  "5" = "blue"
)

# Colors used for the selected models.
model_colors <- model_colors_all[as.character(H_to_compare)]

# Common plot size.
plot_width <- 1200
plot_height <- 800
plot_res <- 120

# ------------------------------------------------------------
# Loading fitted models
# This block loads the posterior samples from Task 1.
# ------------------------------------------------------------

# Empty list to store the fitted models.
fitted_models <- list()

# Load each selected fitted model.
for (H in H_to_compare) {
  
  # Select the correct Task 1 fit file.
  if (H == 1) {
    fit_path <- file.path(task_1_fits_dir, "fit_beta_single_H1.rds")
  } else {
    fit_path <- file.path(task_1_fits_dir, paste0("fit_beta_mixture_H", H, ".rds"))
  }
  
  # Load and store the fitted model.
  fitted_models[[paste0("H=", H)]] <- readRDS(fit_path)
}

# ------------------------------------------------------------
# Posterior predictive density
# This section computes posterior predictive densities
# for the selected fitted models.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Definition of density function
# This function computes the posterior predictive density
# for a Beta model or a Beta mixture model with H components.
# ------------------------------------------------------------

posterior_predictive_density <- function(grid, samps, H) {
  
  # Convert MCMC samples to a matrix.
  mat <- as.matrix(samps)
  
  # Handle the single Beta model.
  if (H == 1) {
    
    predictive_density <- sapply(grid, function(g) {
      mean(dbeta(g, mat[, "alpha"], mat[, "beta"]))
    })
    
    return(predictive_density)
  }
  
  # Create empty matrices for mixture parameters.
  alpha <- matrix(NA, nrow = nrow(mat), ncol = H)
  beta  <- matrix(NA, nrow = nrow(mat), ncol = H)
  p     <- matrix(NA, nrow = nrow(mat), ncol = H)
  
  # Extract alpha, beta, and mixture weights for each component.
  for (h in 1:H) {
    alpha[, h] <- mat[, paste0("alpha[", h, "]")]
    beta[, h]  <- mat[, paste0("beta[", h, "]")]
    p[, h]     <- mat[, paste0("p[", h, "]")]
  }
  
  # Average the mixture density over posterior draws.
  predictive_density <- sapply(grid, function(g) {
    
    dens_s <- rep(0, nrow(mat))
    
    for (h in 1:H) {
      dens_s <- dens_s + p[, h] * dbeta(g, alpha[, h], beta[, h])
    }
    
    mean(dens_s)
  })
  
  predictive_density
}

# ------------------------------------------------------------
# Computation of predictive densities
# This block computes and saves the predictive density
# for each selected model.
# ------------------------------------------------------------

# Empty list to store all predictive densities.
predictive_densities <- list()

# Compute the predictive density for each selected model.
for (H in H_to_compare) {
  
  # Extract the fitted model.
  fit_H <- fitted_models[[paste0("H=", H)]]
  
  # Compute the posterior predictive density.
  density_H <- posterior_predictive_density(
    grid = grid_x,
    samps = fit_H$samples,
    H = H
  )
  
  # Store grid and density.
  predictive_densities[[paste0("H=", H)]] <- data.frame(
    H = H,
    grid_x = grid_x,
    predictive_density = density_H
  )
  
  # Save the predictive density for this model.
  saveRDS(
    predictive_densities[[paste0("H=", H)]],
    file = file.path(fits_dir, paste0("posterior_predictive_density_H", H, ".rds"))
  )
}

# Save all predictive densities together.
saveRDS(
  predictive_densities,
  file = file.path(fits_dir, "posterior_predictive_densities_selected_H.rds")
)

# ------------------------------------------------------------
# Posterior predictive data simulation
# This section simulates replicated datasets from the
# posterior predictive distribution of the selected models.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Definition of simulation function
# This function simulates one replicated dataset from a
# Beta model or a Beta mixture model with H components.
# ------------------------------------------------------------

simulate_posterior_predictive_data <- function(samps, H, n, output_path = NULL) {
  
  # Convert MCMC samples to a matrix.
  mat <- as.matrix(samps)
  
  # Randomly select one posterior draw.
  s <- sample(seq_len(nrow(mat)), size = 1)
  
  # Handle the single Beta model.
  if (H == 1) {
    
    x_rep <- rbeta(
      n = n,
      shape1 = mat[s, "alpha"],
      shape2 = mat[s, "beta"]
    )
    
  } else {
    
    # Extract parameters from the selected posterior draw.
    alpha_s <- numeric(H)
    beta_s  <- numeric(H)
    p_s     <- numeric(H)
    
    for (h in 1:H) {
      alpha_s[h] <- mat[s, paste0("alpha[", h, "]")]
      beta_s[h]  <- mat[s, paste0("beta[", h, "]")]
      p_s[h]     <- mat[s, paste0("p[", h, "]")]
    }
    
    # Simulate mixture components.
    z_rep <- sample(
      x = 1:H,
      size = n,
      replace = TRUE,
      prob = p_s
    )
    
    # Simulate replicated observations from the selected components.
    x_rep <- numeric(n)
    
    for (i in 1:n) {
      x_rep[i] <- rbeta(
        n = 1,
        shape1 = alpha_s[z_rep[i]],
        shape2 = beta_s[z_rep[i]]
      )
    }
  }
  
  # Store the replicated dataset.
  x_rep_df <- data.frame(
    H = H,
    observation = seq_len(n),
    x_rep = x_rep
  )
  
  # Save the replicated dataset if an output path is provided.
  if (!is.null(output_path)) {
    saveRDS(x_rep_df, file = output_path)
  }
  
  x_rep_df
}

# ------------------------------------------------------------
# Computation of replicated datasets
# This block simulates and saves one replicated dataset
# for each selected model.
# ------------------------------------------------------------

# Empty list to store replicated datasets.
replicated_datasets <- list()

# Simulate one replicated dataset for each selected model.
for (H in H_to_compare) {
  
  # Extract the fitted model.
  fit_H <- fitted_models[[paste0("H=", H)]]
  
  # Simulate and save one replicated dataset.
  replicated_datasets[[paste0("H=", H)]] <- simulate_posterior_predictive_data(
    samps = fit_H$samples,
    H = H,
    n = n,
    output_path = file.path(fits_dir, paste0("posterior_predictive_sample_H", H, ".rds"))
  )
}

# Save all replicated datasets together.
saveRDS(
  replicated_datasets,
  file = file.path(fits_dir, "posterior_predictive_samples_selected_H.rds")
)

# ------------------------------------------------------------
# Plot preparation
# This block defines empirical and replicated densities
# and sets common plot scales.
# ------------------------------------------------------------

# Estimate the empirical kernel density.
empirical_density <- density(x, from = 0, to = 1)

# Estimate kernel densities for replicated datasets.
replicated_densities <- lapply(replicated_datasets, function(d) {
  density(d$x_rep, from = 0, to = 1)
})

# Define one common y-axis limit for all plots.
common_max_y <- max(
  empirical_density$y,
  unlist(lapply(predictive_densities, function(d) d$predictive_density)),
  unlist(lapply(replicated_densities, function(d) d$y))
)

# Common plot limits.
common_xlim <- c(0.4, 1)
common_ylim <- c(0, 1.05 * common_max_y)

# Larger and cleaner plot size.
plot_width <- 1800
plot_height <- 1100
plot_res <- 160

# Common graphical settings.
set_plot_style <- function() {
  par(
    mar = c(4.5, 4.8, 3.2, 1.2),
    mgp = c(2.6, 0.8, 0),
    cex.main = 1.25,
    cex.lab = 1.15,
    cex.axis = 1.05,
    cex = 1.05
  )
}

# Common histogram style.
plot_empirical_histogram <- function(main_title) {
  hist(
    x,
    breaks = 35,
    freq = FALSE,
    col = "grey85",
    border = "white",
    xlim = common_xlim,
    ylim = common_ylim,
    xlab = "Hourly load factor",
    ylab = "Density",
    main = main_title
  )
  
  lines(
    empirical_density,
    col = "black",
    lwd = 2.5
  )
  
  rug(x, col = "grey30")
}

# ------------------------------------------------------------
# Empirical distribution plot
# This block saves the empirical histogram and density.
# ------------------------------------------------------------

png(
  filename = file.path(plots_dir, "empirical_distribution.png"),
  width = plot_width,
  height = plot_height,
  res = plot_res
)

set_plot_style()

plot_empirical_histogram(
  main_title = "Observed distribution of hourly load factor"
)

legend(
  "topleft",
  legend = "Empirical kernel density",
  col = "black",
  lwd = 2.5,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# Individual posterior predictive density plots
# This block saves one plot for each selected model.
# ------------------------------------------------------------

for (H in H_to_compare) {
  
  density_H <- predictive_densities[[paste0("H=", H)]]
  color_H <- model_colors_all[as.character(H)]
  
  png(
    filename = file.path(plots_dir, paste0("posterior_predictive_density_H", H, ".png")),
    width = plot_width,
    height = plot_height,
    res = plot_res
  )
  
  set_plot_style()
  
  plot_empirical_histogram(
    main_title = paste0("Predictive density check, H = ", H)
  )
  
  lines(
    density_H$grid_x,
    density_H$predictive_density,
    col = color_H,
    lwd = 2.8
  )
  
  legend(
    "topleft",
    legend = c(
      "Empirical kernel density",
      paste0("Posterior predictive density, H = ", H)
    ),
    col = c("black", color_H),
    lwd = c(2.5, 2.8),
    bty = "n"
  )
  
  dev.off()
}

# ------------------------------------------------------------
# Combined posterior predictive density plot
# This block saves the plot with all selected predictive densities.
# ------------------------------------------------------------

png(
  filename = file.path(plots_dir, "posterior_predictive_densities_all_H.png"),
  width = plot_width,
  height = plot_height,
  res = plot_res
)

set_plot_style()

plot_empirical_histogram(
  main_title = "Predictive density comparison across models"
)

for (H in H_to_compare) {
  
  density_H <- predictive_densities[[paste0("H=", H)]]
  color_H <- model_colors_all[as.character(H)]
  
  lines(
    density_H$grid_x,
    density_H$predictive_density,
    col = color_H,
    lwd = 2.8
  )
}

legend(
  "topleft",
  legend = c(
    "Empirical kernel density",
    paste0("Posterior predictive density, H = ", H_to_compare)
  ),
  col = c("black", model_colors),
  lwd = c(2.5, rep(2.8, length(H_to_compare))),
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# Individual replicated dataset plots
# This block saves one plot for each selected replicated dataset.
# ------------------------------------------------------------

for (H in H_to_compare) {
  
  density_rep_H <- replicated_densities[[paste0("H=", H)]]
  color_H <- model_colors_all[as.character(H)]
  
  png(
    filename = file.path(plots_dir, paste0("posterior_predictive_replicated_dataset_H", H, ".png")),
    width = plot_width,
    height = plot_height,
    res = plot_res
  )
  
  set_plot_style()
  
  plot_empirical_histogram(
    main_title = paste0("Replicated dataset check, H = ", H)
  )
  
  lines(
    density_rep_H,
    col = color_H,
    lwd = 2.8
  )
  
  legend(
    "topleft",
    legend = c(
      "Empirical kernel density",
      paste0("Replicated dataset density, H = ", H)
    ),
    col = c("black", color_H),
    lwd = c(2.5, 2.8),
    bty = "n"
  )
  
  dev.off()
}

# ------------------------------------------------------------
# Combined replicated dataset plot
# This block saves the plot with all selected replicated densities.
# ------------------------------------------------------------

png(
  filename = file.path(plots_dir, "posterior_predictive_replicated_datasets_all_H.png"),
  width = plot_width,
  height = plot_height,
  res = plot_res
)

set_plot_style()

plot_empirical_histogram(
  main_title = "Replicated dataset comparison across models"
)

for (H in H_to_compare) {
  
  density_rep_H <- replicated_densities[[paste0("H=", H)]]
  color_H <- model_colors_all[as.character(H)]
  
  lines(
    density_rep_H,
    col = color_H,
    lwd = 2.8
  )
}

legend(
  "topleft",
  legend = c(
    "Empirical kernel density",
    paste0("Replicated dataset density, H = ", H_to_compare)
  ),
  col = c("black", model_colors),
  lwd = c(2.5, rep(2.8, length(H_to_compare))),
  bty = "n"
)

dev.off()