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
  a_beta <- rnorm_truncated(0.1, 0.01)
  a_mu <- sample(mu, 1) * 3
  a_p_obs <- sample(p_obs, 1)
  list(beta=a_beta, gamma=0.249, mu=a_mu, p_obs=a_p_obs)
}

f_R0 <- function(parameters) {
  parameters$beta * 79 * 3 / parameters$mu * (parameters$gamma / (parameters$gamma + parameters$mu))^3
}

ndraws <- 1000
par_names <- c("beta", "gamma", "mu", "p_obs")
prior_draws <- prior_predictive_R0(ndraws, prior_sample, f_R0,
                                   par_names)
qs <- round(quantile(prior_draws$R0, c(0.25, 0.75)), 1)
p <- qplot(prior_draws$R0) +
  xlab("R0") +
  ylab("Count") +
  annotate("text", x=10, y=100, label=paste0("R0: ", paste(as.list(qs), collapse="-")))
ggsave("outputs/abc_prior_predictive_R0_many_states.pdf", p)
saveRDS(prior_draws, "data/processed/prior_draws_many_states.rds")


max_tries <- 100000
rmse_threshold <- 5
data=c(Dobs=abc_data$df_real$D_obs, D_final_obs=abc_data$D_final_real)
initial_states <- c(S=79, E0_0=1, E0_1=0, E0_2=0,
                    E1_0=0, E1_1=0, E1_2=0,
                    E2_0=0, E2_1=0, E2_2=0,
                    I_0=0, I_1=0, I_2=0,
                    D_unobs=0, D_obs=0)
draws <- abc(ndraws, rmse_threshold,
             simulator,
             prior_sample,
             initial_states,
             data,
             max_tries, parameter_names=c("beta", "gamma", "mu", "p_obs"),
             print_to_screen = T,
             tries_between_prints = 200,
             model_manystates)
draws <- draws %>%
  mutate(R0=beta * 79 * 3 / mu * (gamma / (gamma + mu))^3)
saveRDS(draws, "data/processed/abc_draws_many_states.rds")

qs <- round(quantile(draws$R0, c(0.25, 0.75)), 1)
p <- qplot(draws$R0) +
  xlab("R0") +
  ylab("Count") +
  annotate("text", x=10, y=100, label=paste0("R0: ", paste(as.list(qs), collapse="-")))
ggsave("outputs/abc_posterior_predictive_R0_many_states.pdf", p)

# posterior predictive
draws <- draws %>% 
  select(-rmse, -R0)
nreps <- 400
compartments <- c("S",
                  "E0_0", "E0_1", "E0_2",
                  "E1_0", "E1_1", "E1_2",
                  "E2_0", "E2_1", "E2_2",
                  "I_0", "I_1", "I_2",
                  "D_unobs", "D_obs")
big_df <- posterior_predictive_sims(nreps, model_manystates, draws, initial_states,
                                    compartments)

g <- plot_posterior_predictive_check(big_df, abc_data)

ggsave("outputs/abc_posterior_predictive_many_states.pdf", g)
