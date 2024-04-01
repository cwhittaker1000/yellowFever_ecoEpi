library(tidyverse)

df_real <- readRDS("data/processed/epidemic_series_last_part.rds")

days <- readRDS("data/processed/monkey_death_days.rds")
median_death_days <- round(median(days))

# calculate corpses reported before truncation
df_full <- readRDS("data/processed/epidemic_series.rds")
first_day_reported <- 56 - median_death_days
n_corpses_before <- df_full %>% 
  filter(day < first_day_reported) %>% 
  summarise(corpses_reported=sum(corpses_reported)) %>% 
  pull(corpses_reported)

# assume that only n_corpses before have died and only one monkey
# infected at start
original_monkeys <- 84
max_monkeys <- original_monkeys - n_corpses_before + 1
initial_states <- list(S=original_monkeys-n_corpses_before, I=1, D=0, D_obs=0)
D_final_real <- max_monkeys

abc_data <- list(df_real=df_real,
                 initial_states=initial_states,
                 D_final_real=D_final_real,
                 max_monkeys=max_monkeys)
saveRDS(abc_data, "data/processed/abc_data.rds")