# recreates table 1 from "The susceptibility of Howler monkeys to yellow fever virus",
# Laemmert & Klumm, 1950.

## Loading required libraries
library(tidyverse); library(rstan)

## Interpret (e.g.) 6th day after infection as being 5 days post-infection
days_death_post_infection <- c(5, 8, 2, 10, 6, 4, 3, 4, 3, 9)
## Interpret - "trace" doesn't count; if detectable on 1st day, then 0.5; count 6th day after infection as being 5 days post-infection
days_virus_post_infection <- c(2, 2, 1, 2, 2, 1, 1, 1, 0.5, 2)

df <- tibble(id=seq_along(days_virus_post_infection),
             days_virus_post_infection=days_virus_post_infection,
             days_death_post_infection=days_death_post_infection)

saveRDS(df, "data/processed/laemmert.rds")