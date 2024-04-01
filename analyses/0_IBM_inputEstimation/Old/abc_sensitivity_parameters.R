# creates parameters for ABC sensitivity analysis

library(tidyverse)

beta_mu <- seq(0.003, 0.009, 0.002)
beta_sigma <- c(0.01, 0.05)

params_df <- expand_grid(beta_mu, beta_sigma)
saveRDS(params_df, "data/processed/abc_sensitivity_parameters.rds")
