
library(tidyverse)

days <- readRDS("data/processed/monkey_death_days.rds")
median_death_days <- round(median(days))

df_real <- readRDS("data/processed/epidemic_series.rds")
first_day_reported <- 56 - median_death_days
first_date_reported <- df_real$date[first_day_reported]

g <- df_real %>% 
  ggplot(aes(x=date, y=corpses_reported)) +
  geom_col() +
  geom_vline(xintercept = first_date_reported,
             linetype=2) +
  xlab("Date") +
  ylab("Number of corpses reported")
ggsave("outputs/cutoff_dates.pdf", g, width = 10, height = 6)