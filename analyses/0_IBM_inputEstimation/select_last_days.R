# selects days beginning towards the end of the series (where there is
# an outbreak). Here we assume that the monkey which died on the
# 2017-12-01 was infectious a time before according to the time
# taken for a monkey to die from the monkey_death analysis

library(tidyverse)

days <- readRDS("data/processed/monkey_death_days.rds")
median_death_days <- round(median(days))

df_real <- readRDS("data/processed/epidemic_series.rds")
first_day_reported <- 56 - median_death_days

df_real <- df_real %>% 
  filter(day >= first_day_reported) %>% 
  mutate(day=seq_along(date) - 1) %>% 
  mutate(D_obs=cumsum(corpses_reported))

saveRDS(df_real, "data/processed/epidemic_series_last_part.rds")