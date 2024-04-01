library(tidyverse)

# consider only alouatta in PEAL outbreak
df <- read.csv("data/raw/PEAL PEC 2017-2018 FINAL.csv") %>% 
  rename_all(tolower) %>% 
  filter(peal_outbreak %in% c("yes_lab_positive", "yes_lab_negative")) %>% 
  filter(popular_name == "howler_monkey")

tmp <- df %>% 
  filter(date_death != "unknown") %>% 
  select(date, date_death, date_notification)
