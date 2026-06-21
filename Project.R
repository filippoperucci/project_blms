library(rjags)
library(coda)
library(loo)
library(dplyr)

rm(list = ls())

lf <- read.csv("hourly_load_factor.csv")

y <- pmin(pmax(lf$load_factor, 1e-3), 1 - 1e-3)
N <- length(y)

summary(y)
hist(lf$load_factor, breaks = 35, freq = FALSE,
     col = "tomato", border = "white", 
     xlab = "Load Factor (hourly/daily peak)",
     main = "Hourly city load factor - empirical density (n = 2000)")
lines(density(lf$load_factor,from = 0.01, to = 0.99),
      col = "blue", lwd = 2)


#helpers
extract_ll <- function(samps) {
  do.call(rbind,
          lapply(samps, function(ch)
            ch[, grep("^log_lik\\[", colnames(ch)), drop = FALSE]))
}

ic_row <- function(w, nm) {
  e <- w$estimates
  data.frame(
    Model = nm,
    ELPD  = e["elpd_waic", "Estimate"],
    p_eff = e["p_waic",    "Estimate"],
    WAIC  = e["waic",      "Estimate"],
    SE    = e["waic",      "SE"]
  )
}

# Post-hoc relabelling of mixture draws by increasing component "key" (e.g. the
# mean). Returns the posterior draws with columns reordered per iteration so
# that key[1] <= key[2] <= ... This tames label switching for reporting without
# touching the (permutation-invariant) fit.
relabel_by <- function(samp_mat, H, key_cols, move_groups) {
  out <- samp_mat
  for (i in seq_len(nrow(samp_mat))) {
    o <- order(samp_mat[i, key_cols])
    for (g in move_groups) out[i, g] <- samp_mat[i, g[o]]
  }
  out
} 

#general beta mixture
model_beta_mix <- "
model {

  # Likelihood with latent allocation + marginal log-likelihood for WAIC/LOO.
  for (i in 1:N) {
    z[i] ~ dcat(p[])
    y[i] ~ dnorm(mu[z[i]], tau[z[i]])

    # Per-component log mixture contribution, log(p_h) + log N(y_i | mu_h, tau_h).
    # We use logdensity.norm (exact, no manual exp/sqrt) and combine the H
    # components with a numerically safe log-sum-exp. exp() is a *scalar*
    # function in JAGS, so the exponential is taken element-by-element in the
    # loop and only then summed.
    for (h in 1:H) {
      lcomp[i,h] <- log(p[h]) + logdensity.norm(y[i], mu[h], tau[h])
    }
    mx[i] <- max(lcomp[i,])
    for (h in 1:H) {
      e[i,h] <- exp(lcomp[i,h] - mx[i])
    }
    log_lik[i] <- mx[i] + log(sum(e[i,]))
  }

  # Priors on the component parameters.
  # The data are standardised, so we put a weakly-informative prior directly on
  # the component standard deviations (avoids the spike-at-zero pathology of a
  # vague Gamma(0.01,0.01) on the precision, which is what produced
  # 'Invalid parent values for sigma[.]').
  for (h in 1:H) {
    mu[h]    ~ dnorm(0, 0.01)
    sigma[h] ~ dunif(0.05, 5)          # standardised data => sd in a sane range
    tau[h]   <- 1 / (sigma[h] * sigma[h])
  }

  # Dirichlet prior on the weights
  p[1:H] ~ ddirich(a[])
}
"

#initial values
beta_inits <- function(y, H){
  
  function(){
    
    qs <- as.numeric(
      quantile(y, probs = seq(0.1, 0.9, length.out = H))
    )
    
    list(
      alpha = pmax(1, 20*qs),
      beta  = pmax(1, 20*(1-qs)),
      p      = rep(1/H, H)
    )
  }
}

#fit H
H_grid <- 2:5

fits <- vector("list", length(H_grid))

for(k in seq_along(H_grid)) {
  
  H <- H_grid[k]
  Hk <- as.integer(H_grid[k])
  
  dataList <- list(
    y = y,
    N = length(y),
    H = Hk,
    a = rep(1, Hk)
  )
  stopifnot(length(dataList$a) == dataList$H)
  str(dataList$a)
  length(dataList$a)
  Hk
  jm <- jags.model(
    textConnection(model_beta_mix),
    data = dataList,
    inits = beta_inits(y, H),
    n.chains = 3,
    n.adapt = 2000
  )
  
  update(jm, 5000)
  
  fits[[k]] <- coda.samples(
    jm,
    variable.names =
      c("alpha",
        "beta",
        "p",
        paste0("log_lik[",1:N,"]")),
    n.iter = 10000,
    thin = 5
  )
}

#WAIC model selection
names(fits) <- paste0("H=", H_grid)
waic_list <- lapply(fits, function(s) {
  ll <- extract_ll(s)
  ll <- t(ll)   # <-- transposition
  loo::waic(ll)
})

apply(ll, 2, sd)[1:20]
plot(apply(ll, 1, function(x) sd(x)))

waic_table <-
  bind_rows(Map(ic_row,
                waic_list,
                names(waic_list))) |>
  arrange(WAIC)

waic_table

#compare models
loo::loo_compare(waic_list)

waic_table
#WAIC curve
waic_curve <-
  bind_rows(
    Map(
      function(w, nm)
        data.frame(
          H = as.integer(sub("H=","",nm)),
          WAIC = w$estimates["waic","Estimate"],
          SE   = w$estimates["waic","SE"]
        ),
      waic_list,
      names(waic_list)
    )
  )

plot(
  waic_curve$H,
  waic_curve$WAIC,
  type = "b",
  pch = 19,
  xlab = "Number of components (H)",
  ylab = "WAIC"
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

