library(tidyverse)
library(mvtnorm)
source("src/R/helper.R")

# generate draws
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

nreplicates <- 1000
m_results <- matrix(nrow = nreplicates, ncol = 4)
for(i in 1:nreplicates) {
  vals <- prior_sample()
  R0 <- vals$beta / vals$mu * 84
  m_results[i, ] <- c(vals$beta, vals$mu, vals$p_obs, R0)
}

colnames(m_results) <- c("beta", "mu", "p_obs", "R0")
m_results <- as.data.frame(m_results)
saveRDS(m_results, "data/processed/abc_prior_predictive.rds")

g <- qplot(m_results$R0) +
  scale_x_log10() +
  xlab("R0") +
  ylab("Count")
ggsave("outputs/abc_prior_predictive.pdf", g,
       width = 10, height = 6)