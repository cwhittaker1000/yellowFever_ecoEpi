library(tidyverse)
library(yfmonkeyabc)
source("src/R/helper.R")
source("src/R/abc_helpers.R")

abc_data <- readRDS("data/processed/abc_data.rds")

# get data from monkey_death and observation_time to act as
# priors
days <- readRDS("data/processed/monkey_death_days.rds")
mu <- 1 / days
fit <- readRDS("data/processed/observation_time_fit.rds")
days_obs <- rstan::extract(fit, "days_simulated")[[1]]
p_obs <- 1 / days_obs

prior_sample <- function() {
  a_beta <- rnorm_truncated(0.005, 0.01)
  a_mu <- sample(mu, 1)
  a_p_obs <- sample(p_obs, 1)
  list(beta=a_beta, mu=a_mu, gamma=0.1, p_obs=a_p_obs)
}

ndraws <- 1000
max_tries <- 1000
rmse_threshold <- 20
data=c(Dobs=abc_data$df_real$D_obs, D_final_obs=abc_data$D_final_real)
initial_states <- c(S=79, E=1, I=0, D_unobs=0, D_obs=0)
draws <- abc(ndraws, rmse_threshold,
             simulator,
             prior_sample,
             initial_states,
             data,
             max_tries, parameter_names=c("beta", "gamma", "mu", "p_obs"),
             print_to_screen = T,
             tries_between_prints = 200,
             model_seid)
draws <- draws %>%
  mutate(R0=beta / mu * abc_data$max_monkeys)
saveRDS(draws, "data/processed/abc_draws.rds")
