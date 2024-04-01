library(tidyverse)
library(rstan)
options(mc.cores=4)
rstan_options(auto_write=TRUE)

df <- readRDS("data/processed/epidemic_series.rds") %>% 
  mutate(day=seq_along(date))
model <- stan_model("src/stan/white_pagano.stan")

# assume first case of outbreak is on day 65
df_short <- df %>% 
  filter(day >= 56)

stan_data <- list(T = nrow(df_short),
                  N = df_short$cases,
                  k = 30,
                  alpha_R0=2,
                  beta_R0=1,
                  a_1=5,
                  b_1=1,
                  a_2=1,
                  b_2=2)
fit <- sampling(model, stan_data, iter=2000, chains=4)

R0 <- extract(fit, "R0")[[1]]
mean(R0)
quantile(R0)
qplot(R0)
