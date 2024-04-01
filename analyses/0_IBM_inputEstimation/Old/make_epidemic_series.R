library(tidyverse)

df <- readRDS("data/processed/cleaned_data.rds")

# generate a series that runs over both reported corpses and
# estimated death dates
first_date_death_estimated <- min(df$date_death_estimated)
last_date_reported <- max(df$date_corpse_reported)
dates <- seq(first_date_death_estimated, last_date_reported, 1)

# determine number of estimated deaths and corpses reported
# on each day
deaths <- vector(length = length(dates))
corpses <- vector(length = length(dates))
for(i in seq_along(dates)) {
  df_deaths <- df %>% 
    filter(as.character(date_death_estimated)==dates[i])
  deaths[i] <- nrow(df_deaths)
  df_corpses <- df %>% 
    filter(as.character(date_corpse_reported)==dates[i])
  corpses[i] <- nrow(df_corpses)
}

stopifnot(sum(deaths)==sum(corpses))
stopifnot(sum(deaths)==nrow(df))

time_series <- tibble(date=as.Date(dates),
                      deaths=deaths,
                      corpses_reported=corpses) %>% 
  mutate(day=seq_along(date))
saveRDS(time_series, "data/processed/epidemic_series.rds")