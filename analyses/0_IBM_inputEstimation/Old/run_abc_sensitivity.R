library(tidyverse)
library(yfmonkeyabc)
source("src/R/helper.R")
args <- commandArgs(trailingOnly = TRUE)

params_df <- readRDS("data/processed/abc_sensitivity_parameters.rds")
id <- args[1]
params_df <- params_df[id, ]
beta_mu <- params_df$beta_mu
beta_sigma <- params_df$beta_sigma

abc_data <- readRDS("data/processed/abc_data.rds")

# get data from monkey_death and observation_time to act as
# priors
days <- readRDS("data/processed/monkey_death_days.rds")
mu <- 1 / days
fit <- readRDS("data/processed/observation_time_fit.rds")
days_obs <- rstan::extract(fit, "days_simulated")[[1]]
p_obs <- 1 / days_obs

prior_sample <- function() {
  a_beta <- rnorm_truncated(beta_mu, beta_sigma)
  a_mu <- sample(mu, 1)
  a_p_obs <- sample(p_obs, 1)
  list(beta=a_beta, mu=a_mu, p_obs=a_p_obs)
}

ndraws <- 1000
max_tries <- 250000
rmse_threshold <- 6
draws <- abc(ndraws, rmse_threshold,
             abc_data$initial_states,
             abc_data$df_real,
             abc_data$D_final_real,
             prior_sample,
             max_tries, print_to_screen = T)
draws <- draws %>%
  mutate(R0=beta / mu * abc_data$max_monkeys)
filename <- paste0("data/processed/abc_draws_sensitivity_",
                   id, ".rds")
saveRDS(draws, filename)
