library(tidyverse)
library(rstan)
options(mc.cores=4)
rstan_options(auto_write=TRUE)

df <- readRDS("data/processed/laemmert.rds") %>% 
  mutate(truncated=if_else(is.na(days_post_infection), 1, 0)) %>% 
  mutate(days_post_infection=if_else(is.na(days_post_infection), -99, days_post_infection))

model <- stan_model("src/stan/monkey_death.stan")


data_stan <- list(N=nrow(df),
                  days=df$days_post_infection,
                  truncated=df$truncated,
                  days_truncated=10,
                  a_1=1,
                  a_2=0.5,
                  b_1=1,
                  b_2=2)

fit <- sampling(model, data=data_stan, iter=2000, chains=4)
saveRDS(fit, "data/processed/monkey_death_fit.rds")
