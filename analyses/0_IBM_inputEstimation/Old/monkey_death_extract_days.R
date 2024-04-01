library(rstan)
library(tidyverse)

fit <- readRDS("data/processed/monkey_death_fit.rds")

days <- rstan::extract(fit, "days_simulated")[[1]]
saveRDS(days, "data/processed/monkey_death_days.rds")
