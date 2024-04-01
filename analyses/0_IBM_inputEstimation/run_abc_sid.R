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
  list(beta=a_beta, mu=a_mu, p_obs=a_p_obs)
}

f_R0 <- function(parameters) {
  parameters$beta / parameters$mu * abc_data$max_monkeys
}

ndraws <- 1000
par_names <- c("beta", "mu", "p_obs")
prior_draws <- prior_predictive_R0(ndraws, prior_sample, f_R0,
                                   par_names)
qs <- round(quantile(prior_draws$R0, c(0.25, 0.75)), 1)
p <- qplot(prior_draws$R0) +
  xlab("R0") +
  ylab("Count") +
  annotate("text", x=5, y=30, label=paste0("R0: ", paste(as.list(qs), collapse="-")))
ggsave("outputs/abc_prior_predictive_R0_sid.pdf", p)
saveRDS(prior_draws, "data/processed/prior_draws_sid.rds")

max_tries <- 200000
rmse_threshold <- 5
data=c(Dobs=abc_data$df_real$D_obs, D_final_obs=abc_data$D_final_real)
initial_states <- c(S=79, I=1, D_unobs=0, D_obs=0)
draws <- abc(ndraws, rmse_threshold,
             simulator,
             prior_sample,
             initial_states,
             data,
             max_tries, parameter_names=par_names,
             print_to_screen = T,
             tries_between_prints = 200,
             model_sid)
draws <- draws %>%
  mutate(R0=beta / mu * abc_data$max_monkeys)
saveRDS(draws, "data/processed/abc_draws_sid.rds")

qs <- round(quantile(draws$R0, c(0.25, 0.75)), 1)
p <- qplot(draws$R0) +
  xlab("R0") +
  ylab("Count") +
  annotate("text", x=5, y=30, label=paste0("R0: ", paste(as.list(qs), collapse="-")))
ggsave("outputs/abc_posterior_predictive_R0_sid.pdf", p)

# posterior predictive
draws <- draws %>% 
  select(-rmse, -R0)
nreps <- 400
compartments <- c("S", "I", "D_unobs", "D_obs")
big_df <- posterior_predictive_sims(nreps, model_sid, draws, initial_states,
                          compartments)

g <- plot_posterior_predictive_check(big_df, abc_data)

ggsave("outputs/abc_posterior_predictive_sid.pdf", g)
