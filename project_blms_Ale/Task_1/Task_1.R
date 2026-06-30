# ------------------------------------------------------------
# Package setup
# This block installs missing R packages and checks that JAGS
# is available on the system.
# ------------------------------------------------------------

# List of packages needed for the analysis.
required_packages <- c("rjags", "coda", "loo")

# Check which packages are not installed.
missing_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

# Install missing packages, if any.
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Stop the script if JAGS is not installed.
if (Sys.which("jags") == "") {
  stop(
    "JAGS is not installed on this system. ",
    "Please install it before running the script. ",
    "Ubuntu: sudo apt install jags. ",
    "Windows: install JAGS from https://sourceforge.net/projects/mcmc-jags/. ",
    "macOS: install JAGS from https://sourceforge.net/projects/mcmc-jags/ ",
    "or use Homebrew: brew install jags."
  )
}

# Load the required packages without startup messages.
suppressPackageStartupMessages(library(rjags))
suppressPackageStartupMessages(library(coda))
suppressPackageStartupMessages(library(loo))

# Clear the current R environment.
rm(list = ls())

# Set the seed to make the results reproducible.
set.seed(123)

# ------------------------------------------------------------
# Run mode
# If RERUN_JAGS = TRUE, the models are fitted again using JAGS.
# If RERUN_JAGS = FALSE, the previously saved fitted models are loaded.
# ------------------------------------------------------------

RERUN_JAGS <- FALSE

# ------------------------------------------------------------
# Project paths
# This block defines the folder structure used to load data
# and save the Task 1 outputs.
# ------------------------------------------------------------

# Main project folder.
base_dir <- "/home/blackmamba/project_blms/project_blms_Ale"

# Main input and output folders.
data_dir <- file.path(base_dir, "data")
task_dir <- file.path(base_dir, "Task_1")
results_dir <- file.path(task_dir, "results")

# Subfolders for Task 1 results.
fits_dir <- file.path(results_dir, "fits")
plots_dir <- file.path(results_dir, "plots")
diagnostics_dir <- file.path(results_dir, "diagnostics")

# Create result folders if they do not already exist.
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

# Path to the input dataset.
data_path <- file.path(data_dir, "hourly_load_factor.csv")

# ------------------------------------------------------------
# Progress helper
# This function prints clean progress messages during execution.
# ------------------------------------------------------------

# Print a timestamped progress message.
cat_step <- function(message) {
  cat(format(Sys.time(), "%H:%M:%S"), "-", message, "\n")
  flush.console()
}

# Round only the numeric columns of a data frame.
round_numeric_df <- function(df, digits = 4) {
  out <- df
  num_cols <- sapply(out, is.numeric)
  out[num_cols] <- lapply(out[num_cols], round, digits = digits)
  out
}

# Extract log-likelihood matrix: rows = posterior draws, columns = observations.
extract_log_lik <- function(samps) {
  do.call(
    rbind,
    lapply(samps, function(chain) {
      chain[, grep("^log_lik\\[", colnames(chain)), drop = FALSE]
    })
  )
}

# Extract useful WAIC quantities into one row.
waic_row <- function(w, H) {
  e <- w$estimates
  data.frame(
    H = H,
    ELPD = e["elpd_waic", "Estimate"],
    p_eff = e["p_waic", "Estimate"],
    WAIC = e["waic", "Estimate"],
    SE = e["waic", "SE"]
  )
}

# Print the run header.
cat("\n============================================================\n")
cat("HOURLY CITY LOAD FACTOR - RUN STARTED\n")
cat("============================================================\n")

# ------------------------------------------------------------
# Data loading and empirical analysis
# Notation: x = (x_1, ..., x_n), where x_i is the observed load factor.
# ------------------------------------------------------------

cat_step("[1/7] Loading data and computing empirical summaries.")

# Load the dataset.
data <- read.csv(data_path)

# Extract the response variable.
x <- data$load_factor
n <- length(x)

# Check that all observations are inside the Beta support.
stopifnot(all(x > 0), all(x < 1))

# Compute simple empirical summaries.
empirical_summary <- summary(x)
empirical_range <- range(x)

cat_step(paste0("Data loaded: n = ", n, " observations."))

# Print empirical summaries.
cat("\n1. EMPIRICAL ANALYSIS\n")
cat("------------------------------------------------------------\n")
cat("Number of observations:", n, "\n")
cat("Range:", round(empirical_range[1], 4), "-", round(empirical_range[2], 4), "\n")
cat("\nEmpirical summary:\n")
print(empirical_summary)

# Save the empirical distribution plot.
png(
  filename = file.path(plots_dir, "empirical_distribution.png"),
  width = 1000,
  height = 700
)

hist(
  x,
  breaks = 35,
  freq = FALSE,
  col = "tomato",
  border = "white",
  xlab = expression(x[i]),
  ylab = "Density",
  main = "Empirical distribution of the observed sample"
)

lines(
  density(x, from = 0, to = 1),
  col = "blue",
  lwd = 2
)

rug(x)

dev.off()

# ------------------------------------------------------------
# Model fitting or loading
# ------------------------------------------------------------

# Number of mixture components to test.
H_grid <- 2:5

if (RERUN_JAGS) {
  
  # ------------------------------------------------------------
  # Single Beta model
  # We fit x_i | theta ~ Beta(alpha, beta), theta = (alpha, beta).
  # ------------------------------------------------------------
  
  cat_step("[2/7] Fitting single Beta model, H = 1.")
  
  # JAGS model for a single Beta distribution.
  model_beta_single <- "
model {

  # Likelihood
  for (i in 1:n) {
    x[i] ~ dbeta(alpha, beta)
    log_lik[i] <- logdensity.beta(x[i], alpha, beta)
  }

  # Priors
  alpha ~ dgamma(a_alpha, b_alpha)
  beta  ~ dgamma(a_beta, b_beta)
}
"

# Data and hyperparameters passed to JAGS.
data_beta_single <- list(
  x = x,
  n = n,
  a_alpha = 2,
  b_alpha = 0.1,
  a_beta  = 2,
  b_beta  = 0.1
)

# Initial values for the single Beta model.
inits_beta_single <- function() {
  list(
    alpha = 2,
    beta = 2
  )
}

cat_step("Compiling single Beta model.")

# Compile the JAGS model.
jm_beta_single <- jags.model(
  textConnection(model_beta_single),
  data = data_beta_single,
  inits = inits_beta_single,
  n.chains = 3,
  n.adapt = 2000,
  quiet = TRUE
)

cat_step("Running burn-in for single Beta model.")

# Run burn-in iterations.
update(jm_beta_single, 5000, progress.bar = "none")

cat_step("Sampling posterior draws for single Beta model.")

# Draw posterior samples.
samps_beta_single <- coda.samples(
  jm_beta_single,
  variable.names = c("alpha", "beta", paste0("log_lik[", 1:n, "]")),
  n.iter = 10000,
  thin = 5,
  progress.bar = "none"
)

cat_step("Computing WAIC for single Beta model.")

# Compute WAIC for the single Beta model.
log_lik_beta_single <- extract_log_lik(samps_beta_single)
waic_beta_single <- loo::waic(log_lik_beta_single)

# Store the fitted model objects.
fit_beta_single <- list(
  samples = samps_beta_single,
  log_lik = log_lik_beta_single,
  waic = waic_beta_single
)

# Save the single Beta fit.
saveRDS(
  fit_beta_single,
  file = file.path(fits_dir, "fit_beta_single_H1.rds")
)

cat_step("Completed single Beta model, H = 1.")

# ------------------------------------------------------------
# Beta mixture models
# General model for H = 1, ..., 5 components.
# z_i | p ~ Categorical(p_1, ..., p_H)
# x_i | z_i = h, theta_h ~ Beta(alpha_h, beta_h)
# ------------------------------------------------------------

cat_step("[3/7] Fitting Beta mixture models, H = 2, 3, 4, 5.")

# JAGS model for a Beta mixture with H components.
model_beta_mixture <- "
model {

  # Likelihood
  for (i in 1:n) {

    z[i] ~ dcat(p[])
    x[i] ~ dbeta(alpha[z[i]], beta[z[i]])

    # Marginal log-likelihood for WAIC
    for (h in 1:H) {
      lcomp[i,h] <- log(p[h]) + logdensity.beta(x[i], alpha[h], beta[h])
    }

    mx[i] <- max(lcomp[i,])

    for (h in 1:H) {
      ecomp[i,h] <- exp(lcomp[i,h] - mx[i])
    }

    log_lik[i] <- mx[i] + log(sum(ecomp[i,]))
  }

  # Component priors
  for (h in 1:H) {
    alpha[h] ~ dgamma(a_alpha, b_alpha)
    beta[h]  ~ dgamma(a_beta, b_beta)
  }

  # Mixture weights prior
  p[1:H] ~ ddirich(a[])
}
"

# Initial values for the Beta mixture model.
inits_beta_mixture <- function(H) {
  function() {
    qs <- as.numeric(quantile(x, probs = seq(0.2, 0.8, length.out = H)))
    
    list(
      alpha = pmax(1, 5 * qs),
      beta  = pmax(1, 5 * (1 - qs)),
      p = rep(1 / H, H)
    )
  }
}

# Empty lists to store fits and WAIC values.
fits_beta_mixture <- vector("list", length(H_grid))
waic_beta_mixture <- vector("list", length(H_grid))

# Fit one mixture model for each value of H.
for (k in seq_along(H_grid)) {
  
  H <- H_grid[k]
  
  cat_step(paste0("Starting Beta mixture model, H = ", H, "."))
  
  # Data and hyperparameters passed to JAGS.
  data_beta_mixture <- list(
    x = x,
    n = n,
    H = H,
    a = rep(1, H),
    a_alpha = 2,
    b_alpha = 0.1,
    a_beta  = 2,
    b_beta  = 0.1
  )
  
  cat_step(paste0("Compiling Beta mixture model, H = ", H, "."))
  
  # Compile the mixture model.
  jm_beta_mixture <- jags.model(
    textConnection(model_beta_mixture),
    data = data_beta_mixture,
    inits = inits_beta_mixture(H),
    n.chains = 3,
    n.adapt = 2000,
    quiet = TRUE
  )
  
  cat_step(paste0("Running burn-in for Beta mixture model, H = ", H, "."))
  
  # Run burn-in iterations.
  update(jm_beta_mixture, 5000, progress.bar = "none")
  
  cat_step(paste0("Sampling posterior draws for Beta mixture model, H = ", H, "."))
  
  # Draw posterior samples.
  samps_beta_mixture <- coda.samples(
    jm_beta_mixture,
    variable.names = c(
      "alpha",
      "beta",
      "p",
      paste0("log_lik[", 1:n, "]")
    ),
    n.iter = 10000,
    thin = 5,
    progress.bar = "none"
  )
  
  cat_step(paste0("Computing WAIC for Beta mixture model, H = ", H, "."))
  
  # Compute WAIC for the current mixture model.
  log_lik_beta_mixture <- extract_log_lik(samps_beta_mixture)
  waic_H <- loo::waic(log_lik_beta_mixture)
  
  # Store the fitted model.
  fits_beta_mixture[[k]] <- list(
    H = H,
    samples = samps_beta_mixture,
    log_lik = log_lik_beta_mixture,
    waic = waic_H
  )
  
  waic_beta_mixture[[k]] <- waic_H
  
  # Save the current mixture fit.
  saveRDS(
    fits_beta_mixture[[k]],
    file = file.path(fits_dir, paste0("fit_beta_mixture_H", H, ".rds"))
  )
  
  cat_step(
    paste0(
      "Completed Beta mixture model, H = ",
      H,
      ". WAIC = ",
      round(waic_H$estimates["waic", "Estimate"], 3),
      "."
    )
  )
}

} else {
  
  cat_step("[2/7] Loading previously saved fitted models.")
  
  # Load the single Beta model.
  fit_beta_single <- readRDS(
    file = file.path(fits_dir, "fit_beta_single_H1.rds")
  )
  
  waic_beta_single <- fit_beta_single$waic
  
  # Load all mixture models.
  fits_beta_mixture <- vector("list", length(H_grid))
  waic_beta_mixture <- vector("list", length(H_grid))
  
  for (k in seq_along(H_grid)) {
    
    H <- H_grid[k]
    
    fits_beta_mixture[[k]] <- readRDS(
      file = file.path(fits_dir, paste0("fit_beta_mixture_H", H, ".rds"))
    )
    
    waic_beta_mixture[[k]] <- fits_beta_mixture[[k]]$waic
    
    cat_step(paste0("Loaded Beta mixture model, H = ", H, "."))
  }
  
  cat_step("Previously saved fitted models loaded.")
}

# ------------------------------------------------------------
# WAIC model comparison
# Lower WAIC indicates better estimated predictive performance.
# ------------------------------------------------------------

cat_step("[4/7] Comparing models using WAIC.")

# Build the WAIC comparison table.
waic_table <- rbind(
  waic_row(waic_beta_single, H = 1),
  do.call(
    rbind,
    Map(waic_row, waic_beta_mixture, H_grid)
  )
)

# Sort models by WAIC.
waic_table <- waic_table[order(waic_table$WAIC), ]
waic_curve <- waic_table[order(waic_table$H), ]

# Select the best model according to the lowest WAIC.
best_H <- waic_table$H[1]

cat_step(paste0("WAIC comparison completed. Best model: H = ", best_H, "."))

# ------------------------------------------------------------
# MCMC diagnostics
# This section stores trace plots, autocorrelation plots, and effective sample size.
# Burn-in and thinning have already been applied during sampling.
# ------------------------------------------------------------

cat_step("[5/7] Computing MCMC diagnostics.")

# Compute basic MCMC diagnostics for selected parameters.
mcmc_diagnostics <- function(samps, pars) {
  
  pars_available <- pars[pars %in% varnames(samps)]
  
  if (length(pars_available) == 0) {
    stop("None of the requested parameters are available in the MCMC sample.")
  }
  
  par_samps <- samps[, pars_available]
  
  list(
    summary = summary(par_samps),
    effective_size = effectiveSize(par_samps),
    samples = par_samps
  )
}

# Diagnostics for the single Beta model.
diagnostics_beta_single <- mcmc_diagnostics(
  samps = fit_beta_single$samples,
  pars = c("alpha", "beta")
)

# Diagnostics for the mixture models.
diagnostics_beta_mixture <- vector("list", length(fits_beta_mixture))

for (k in seq_along(fits_beta_mixture)) {
  
  H <- fits_beta_mixture[[k]]$H
  
  pars_H <- c(
    paste0("alpha[", 1:H, "]"),
    paste0("beta[", 1:H, "]"),
    paste0("p[", 1:H, "]")
  )
  
  diagnostics_beta_mixture[[k]] <- mcmc_diagnostics(
    samps = fits_beta_mixture[[k]]$samples,
    pars = pars_H
  )
  
  cat_step(paste0("Computed diagnostics for Beta mixture model, H = ", H, "."))
}

names(diagnostics_beta_mixture) <- paste0("H=", H_grid)

# Save all diagnostics.
saveRDS(
  list(
    single_beta = diagnostics_beta_single,
    beta_mixture = diagnostics_beta_mixture
  ),
  file = file.path(diagnostics_dir, "mcmc_diagnostics.rds")
)

cat_step("MCMC diagnostics completed.")

# ------------------------------------------------------------
# Parameter summaries
# This section creates clean posterior summaries for the fitted parameters.
# ------------------------------------------------------------

cat_step("[6/7] Computing posterior parameter summaries.")

# Summarize posterior samples using mean, SD, and credible intervals.
summarize_params <- function(samps, pars) {
  
  pars_available <- pars[pars %in% varnames(samps)]
  mat <- as.matrix(samps[, pars_available])
  
  q <- t(apply(
    mat,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975)
  ))
  
  data.frame(
    Parameter = colnames(mat),
    Mean = colMeans(mat),
    SD = apply(mat, 2, sd),
    Q2.5 = q[, 1],
    Median = q[, 2],
    Q97.5 = q[, 3],
    row.names = NULL
  )
}

# Summarize mixture parameters after ordering components by their mean.
summarize_mixture_params <- function(samps, H) {
  
  mat <- as.matrix(samps)
  
  alpha <- matrix(NA, nrow = nrow(mat), ncol = H)
  beta  <- matrix(NA, nrow = nrow(mat), ncol = H)
  p     <- matrix(NA, nrow = nrow(mat), ncol = H)
  
  # Extract parameters for each component.
  for (h in 1:H) {
    alpha[, h] <- mat[, paste0("alpha[", h, "]")]
    beta[, h]  <- mat[, paste0("beta[", h, "]")]
    p[, h]     <- mat[, paste0("p[", h, "]")]
  }
  
  # Relabel components by increasing component mean alpha / (alpha + beta).
  alpha_ord <- alpha
  beta_ord  <- beta
  p_ord     <- p
  
  component_mean <- alpha / (alpha + beta)
  
  for (s in 1:nrow(mat)) {
    ord <- order(component_mean[s, ])
    alpha_ord[s, ] <- alpha[s, ord]
    beta_ord[s, ]  <- beta[s, ord]
    p_ord[s, ]     <- p[s, ord]
  }
  
  out <- data.frame()
  
  # Build the summary table component by component.
  for (h in 1:H) {
    
    block <- data.frame(
      Component = h,
      Parameter = c("alpha", "beta", "p"),
      Mean = c(mean(alpha_ord[, h]), mean(beta_ord[, h]), mean(p_ord[, h])),
      SD = c(sd(alpha_ord[, h]), sd(beta_ord[, h]), sd(p_ord[, h])),
      Q2.5 = c(
        quantile(alpha_ord[, h], 0.025),
        quantile(beta_ord[, h], 0.025),
        quantile(p_ord[, h], 0.025)
      ),
      Median = c(
        quantile(alpha_ord[, h], 0.5),
        quantile(beta_ord[, h], 0.5),
        quantile(p_ord[, h], 0.5)
      ),
      Q97.5 = c(
        quantile(alpha_ord[, h], 0.975),
        quantile(beta_ord[, h], 0.975),
        quantile(p_ord[, h], 0.975)
      ),
      row.names = NULL
    )
    
    out <- rbind(out, block)
  }
  
  out
}

# Posterior summary for the single Beta model.
param_summary_single <- summarize_params(
  samps = fit_beta_single$samples,
  pars = c("alpha", "beta")
)

# Posterior summaries for the mixture models.
param_summary_mixture <- vector("list", length(fits_beta_mixture))

for (k in seq_along(fits_beta_mixture)) {
  
  H <- fits_beta_mixture[[k]]$H
  
  param_summary_mixture[[k]] <- summarize_mixture_params(
    samps = fits_beta_mixture[[k]]$samples,
    H = H
  )
}

names(param_summary_mixture) <- paste0("H=", H_grid)

cat_step("Posterior parameter summaries completed.")

# ------------------------------------------------------------
# Posterior predictive density
# This section computes the posterior mean fitted density for the best model.
# ------------------------------------------------------------

cat_step("[7/7] Computing posterior predictive density for selected model.")

# Compute the posterior mean density for the single Beta model.
posterior_density_single <- function(grid, samps) {
  
  mat <- as.matrix(samps)
  
  sapply(grid, function(g) {
    mean(dbeta(g, mat[, "alpha"], mat[, "beta"]))
  })
}

# Compute the posterior mean density for a Beta mixture model.
posterior_density_mixture <- function(grid, samps, H) {
  
  mat <- as.matrix(samps)
  
  alpha <- matrix(NA, nrow = nrow(mat), ncol = H)
  beta  <- matrix(NA, nrow = nrow(mat), ncol = H)
  p     <- matrix(NA, nrow = nrow(mat), ncol = H)
  
  # Extract mixture parameters.
  for (h in 1:H) {
    alpha[, h] <- mat[, paste0("alpha[", h, "]")]
    beta[, h]  <- mat[, paste0("beta[", h, "]")]
    p[, h]     <- mat[, paste0("p[", h, "]")]
  }
  
  # Average the mixture density over posterior draws.
  sapply(grid, function(g) {
    dens_s <- rep(0, nrow(mat))
    
    for (h in 1:H) {
      dens_s <- dens_s + p[, h] * dbeta(g, alpha[, h], beta[, h])
    }
    
    mean(dens_s)
  })
}

# Grid used to evaluate the fitted density.
grid_x <- seq(0.001, 0.999, length.out = 300)

# Compute the fitted density for the selected model.
if (best_H == 1) {
  
  best_density <- posterior_density_single(
    grid = grid_x,
    samps = fit_beta_single$samples
  )
  
} else {
  
  best_fit_index <- which(H_grid == best_H)
  
  best_density <- posterior_density_mixture(
    grid = grid_x,
    samps = fits_beta_mixture[[best_fit_index]]$samples,
    H = best_H
  )
}

cat_step("Posterior predictive density completed.")

# ------------------------------------------------------------
# Final report
# This section prints the main results in a clean order.
# ------------------------------------------------------------

# Print the final report header.
cat("\n============================================================\n")
cat("HOURLY CITY LOAD FACTOR - BAYESIAN BETA MODELS\n")
cat("============================================================\n")

# Print MCMC diagnostic summaries.
cat("\n2. MCMC DIAGNOSTICS\n")
cat("------------------------------------------------------------\n")
cat("Burn-in iterations discarded: 5000\n")
cat("Posterior sampling iterations: 10000\n")
cat("Thinning interval: 5\n")
cat("Number of chains: 3\n")

cat("\nEffective sample size - Single Beta model:\n")
print(round(diagnostics_beta_single$effective_size, 2))

for (k in seq_along(diagnostics_beta_mixture)) {
  cat("\nEffective sample size - Beta mixture", names(diagnostics_beta_mixture)[k], ":\n")
  print(round(diagnostics_beta_mixture[[k]]$effective_size, 2))
}

# Save diagnostic plots for the selected model.
cat("\nDiagnostic plots for the selected model have been saved.\n")

if (best_H == 1) {
  
  png(
    filename = file.path(plots_dir, "traceplot_selected_model.png"),
    width = 1200,
    height = 800
  )
  
  traceplot(diagnostics_beta_single$samples)
  
  dev.off()
  
  png(
    filename = file.path(plots_dir, "autocorrelation_selected_model.png"),
    width = 1200,
    height = 800
  )
  
  autocorr.plot(diagnostics_beta_single$samples)
  
  dev.off()
  
} else {
  
  best_fit_index <- which(H_grid == best_H)
  
  png(
    filename = file.path(plots_dir, "traceplot_selected_model.png"),
    width = 1200,
    height = 800
  )
  
  traceplot(diagnostics_beta_mixture[[best_fit_index]]$samples)
  
  dev.off()
  
  png(
    filename = file.path(plots_dir, "autocorrelation_selected_model.png"),
    width = 1200,
    height = 800
  )
  
  autocorr.plot(diagnostics_beta_mixture[[best_fit_index]]$samples)
  
  dev.off()
}

# Print posterior parameter summaries.
cat("\n3. POSTERIOR PARAMETER ESTIMATES\n")
cat("------------------------------------------------------------\n")

cat("\nSingle Beta model, H = 1:\n")
print(round_numeric_df(param_summary_single, 4))

for (k in seq_along(param_summary_mixture)) {
  cat("\nBeta mixture model,", names(param_summary_mixture)[k], ":\n")
  cat("Components are ordered by increasing posterior component mean.\n")
  print(round_numeric_df(param_summary_mixture[[k]], 4))
}

# Save the posterior predictive density plot of the selected model.
cat("\n4. POSTERIOR PREDICTIVE DENSITY\n")
cat("------------------------------------------------------------\n")
cat("Selected model according to WAIC: H =", best_H, "\n")

png(
  filename = file.path(plots_dir, "posterior_predictive_density_selected_model.png"),
  width = 1000,
  height = 700
)

hist(
  x,
  breaks = 35,
  freq = FALSE,
  col = "tomato",
  border = "white",
  xlab = expression(x[i]),
  ylab = "Density",
  main = paste0("Empirical density and fitted posterior density, H = ", best_H)
)

lines(
  density(x, from = 0, to = 1),
  col = "blue",
  lwd = 2
)

lines(
  grid_x,
  best_density,
  col = "black",
  lwd = 2
)

legend(
  "topleft",
  legend = c("Empirical kernel density", "Posterior fitted density"),
  col = c("blue", "black"),
  lwd = 2,
  bty = "n"
)

dev.off()

# Print the WAIC comparison.
cat("\n5. WAIC MODEL COMPARISON\n")
cat("------------------------------------------------------------\n")
cat("Lower WAIC indicates better estimated predictive performance.\n\n")

print(round_numeric_df(waic_table, 3))

cat("\nBest model according to WAIC: H =", best_H, "\n")

# Save WAIC plot as a function of the number of components.
png(
  filename = file.path(plots_dir, "waic_model_comparison.png"),
  width = 1000,
  height = 700
)

plot(
  waic_curve$H,
  waic_curve$WAIC,
  type = "b",
  pch = 19,
  xlab = "Number of components H",
  ylab = "WAIC",
  main = "WAIC comparison for Beta mixture models"
)

arrows(
  waic_curve$H,
  waic_curve$WAIC - waic_curve$SE,
  waic_curve$H,
  waic_curve$WAIC + waic_curve$SE,
  angle = 90,
  code = 3,
  length = 0.05
)

dev.off()

# Print saved file paths.
cat("\n6. SAVED OBJECTS\n")
cat("------------------------------------------------------------\n")
cat("Saved files:\n")
cat("- ", file.path(fits_dir, "fit_beta_single_H1.rds"), "\n", sep = "")

for (H in H_grid) {
  cat("- ", file.path(fits_dir, paste0("fit_beta_mixture_H", H, ".rds")), "\n", sep = "")
}

cat("- ", file.path(diagnostics_dir, "mcmc_diagnostics.rds"), "\n", sep = "")
cat("- ", file.path(plots_dir, "empirical_distribution.png"), "\n", sep = "")
cat("- ", file.path(plots_dir, "traceplot_selected_model.png"), "\n", sep = "")
cat("- ", file.path(plots_dir, "autocorrelation_selected_model.png"), "\n", sep = "")
cat("- ", file.path(plots_dir, "posterior_predictive_density_selected_model.png"), "\n", sep = "")
cat("- ", file.path(plots_dir, "waic_model_comparison.png"), "\n", sep = "")

cat("\nAnalysis completed.\n")


# ------------------------------------------------------------
# Selected model: H = 3
# Posterior summaries and posterior plots
# ------------------------------------------------------------

best_H <- 3
best_fit_index <- which(H_grid == best_H)

cat_step(paste0("Computing posterior summaries and posterior plots for H = ", best_H, "."))

# Extract posterior draws of the selected model.
posterior_mat_raw <- as.matrix(fits_beta_mixture[[best_fit_index]]$samples)

# Extract alpha, beta and p.
alpha <- matrix(NA, nrow = nrow(posterior_mat_raw), ncol = best_H)
beta  <- matrix(NA, nrow = nrow(posterior_mat_raw), ncol = best_H)
p     <- matrix(NA, nrow = nrow(posterior_mat_raw), ncol = best_H)

for (h in 1:best_H) {
  alpha[, h] <- posterior_mat_raw[, paste0("alpha[", h, "]")]
  beta[, h]  <- posterior_mat_raw[, paste0("beta[", h, "]")]
  p[, h]     <- posterior_mat_raw[, paste0("p[", h, "]")]
}

# Order components by increasing posterior component mean.
component_mean <- alpha / (alpha + beta)

alpha_ord <- alpha
beta_ord  <- beta
p_ord     <- p

for (s in 1:nrow(posterior_mat_raw)) {
  ord <- order(component_mean[s, ])
  alpha_ord[s, ] <- alpha[s, ord]
  beta_ord[s, ]  <- beta[s, ord]
  p_ord[s, ]     <- p[s, ord]
}

# Build posterior summary table.
posterior_summary_H3 <- data.frame()

for (h in 1:best_H) {
  
  alpha_h <- alpha_ord[, h]
  beta_h  <- beta_ord[, h]
  p_h     <- p_ord[, h]
  
  block <- data.frame(
    Component = h,
    Parameter = c("alpha", "beta", "p"),
    Mean = c(mean(alpha_h), mean(beta_h), mean(p_h)),
    Variance = c(var(alpha_h), var(beta_h), var(p_h)),
    SD = c(sd(alpha_h), sd(beta_h), sd(p_h)),
    Q2.5 = c(
      quantile(alpha_h, 0.025),
      quantile(beta_h, 0.025),
      quantile(p_h, 0.025)
    ),
    Median = c(
      quantile(alpha_h, 0.5),
      quantile(beta_h, 0.5),
      quantile(p_h, 0.5)
    ),
    Q97.5 = c(
      quantile(alpha_h, 0.975),
      quantile(beta_h, 0.975),
      quantile(p_h, 0.975)
    ),
    row.names = NULL
  )
  
  posterior_summary_H3 <- rbind(posterior_summary_H3, block)
}

# Save posterior summary table.
write.csv(
  posterior_summary_H3,
  file = file.path(diagnostics_dir, "posterior_summary_selected_model_H3.csv"),
  row.names = FALSE
)

# Create data frame for posterior plots.
posterior_draws_H3 <- data.frame(
  alpha_1 = alpha_ord[, 1],
  beta_1  = beta_ord[, 1],
  p_1     = p_ord[, 1],
  alpha_2 = alpha_ord[, 2],
  beta_2  = beta_ord[, 2],
  p_2     = p_ord[, 2],
  alpha_3 = alpha_ord[, 3],
  beta_3  = beta_ord[, 3],
  p_3     = p_ord[, 3]
)

# Save posterior draws.
write.csv(
  posterior_draws_H3,
  file = file.path(diagnostics_dir, "posterior_draws_selected_model_H3.csv"),
  row.names = FALSE
)

# Save one combined posterior plot.
png(
  filename = file.path(plots_dir, "posterior_distributions_selected_model_H3.png"),
  width = 1400,
  height = 1000
)

par(mfrow = c(3, 3))

for (par_name in colnames(posterior_draws_H3)) {
  
  hist(
    posterior_draws_H3[[par_name]],
    breaks = 40,
    freq = FALSE,
    col = "lightgray",
    border = "white",
    xlab = par_name,
    ylab = "Density",
    main = paste("Posterior of", par_name)
  )
  
  lines(
    density(posterior_draws_H3[[par_name]]),
    col = "blue",
    lwd = 2
  )
}

dev.off()

cat_step("Posterior summaries and posterior plots for H = 3 saved.")