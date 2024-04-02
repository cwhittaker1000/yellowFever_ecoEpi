# fits time from death to corpse observed

library(tidyverse)
library(rstan)
options(mc.cores=4)
rstan_options(auto_write=TRUE)

df <- readRDS("data/processed/cleaned_data.rds")

# fit only uninterpolated data
df <- df %>% 
  filter(days_since_death_interpolated != 1)

model <- stan_model("src/stan/observation_time.stan")

data_stan <- list(N=nrow(df),
                  days=df$days_since_death,
                  a_1=2, a_2=1,
                  b_1=3, b_2=3)

fit <- sampling(model, data=data_stan, iter=2000, chains=4)

saveRDS(fit, "data/processed/observation_time_fit.rds")