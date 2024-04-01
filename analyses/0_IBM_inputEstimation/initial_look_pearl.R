library(tidyverse)

df <- read.csv("data/raw/PEAL PEC 2017-2018 ALL.csv") %>% 
  rename_all(tolower) %>% 
  mutate(date=as.Date(date, "%d/%m/%Y")) %>% 
  mutate(date_death=as.Date(date_death, "%d/%m/%Y")) %>% 
  mutate(date_notification=as.Date(date_notification, "%d/%m/%Y")) %>% 
  mutate(location=if_else(location != "Norte", "Mairipora", "Norte"))

table(df$scientific_name)

df %>% 
  group_by(location, scientific_name, date) %>% 
  count() %>% 
  ggplot(aes(x=date, y=n, fill=scientific_name)) +
  geom_col() +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~location)

df %>% 
  ggplot(aes(x=longitude, y=latitude, colour=location)) +
  geom_point() +
  facet_wrap(~scientific_name)
