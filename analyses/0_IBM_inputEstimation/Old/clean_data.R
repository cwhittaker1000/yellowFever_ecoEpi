library(tidyverse)

df <- readxl::read_xlsx("data/raw/EPIZOOTIAS_HORTO _SP_2017.xlsx") %>% 
  mutate(date_corpse_reported=as.Date(Death_repoorting)) %>% 
  select(-Death_repoorting)

# remove alive monkey
df <- df %>% 
  filter(decomposition != "Alive")

# convert unknown days_since_death to medians
df <- df %>% 
  mutate(days_since_death=as.numeric(days_since_death))
median_days <- median(df$days_since_death, na.rm = T)
df <- df %>% 
  mutate(days_since_death_interpolated=if_else(is.na(days_since_death), 1, 0)) %>% 
  mutate(days_since_death=if_else(is.na(days_since_death), median_days, days_since_death)) %>% 
  mutate(date_death_estimated=date_corpse_reported-days_since_death)

saveRDS(df, "data/processed/cleaned_data.rds")