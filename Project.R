rm(list = ls())
lf = read.csv("hourly_load_factor.csv")
dim(lf); head(lf); summary(lf$load_factor)
hist(lf$load_factor, breaks = 35, freq = FALSE,
     col = "tomato", border = "white", 
     xlab = "Load Factor (hourly/daily peak)",
     main = "Hourly city load factor - empirical density (n = 2000)")
lines(density(lf$load_factor,from = 0.01, to = 0.99),
      col = "blue", lwd = 2)


