library(tidyverse)
library(yfmonkeyabc)

abc_data <- readRDS("data/processed/abc_data.rds")
draws <- readRDS("data/processed/abc_draws.rds")

nreps <- 400

for(i in 1:nreps) {
  id <- sample(nrow(draws), 1)
  parameters <- draws[id, ] %>% 
    select(-R0)
  df_sim <- gillespie_daily(83, abc_data$initial_states, parameters) %>% 
    mutate(replicate=i)
  if(i == 1)
    big_df <- df_sim
  else
    big_df <- big_df %>% bind_rows(df_sim)
}

big_df1 <- big_df %>% 
  mutate(type="simulated") %>% 
  bind_rows(abc_data$df_real %>% mutate(time=day, type="real"))

g <- ggplot(big_df1 %>% filter(type=="simulated"), aes(x=time, y=D_obs)) +
  geom_line(aes(group=as.factor(replicate)), alpha=0.5) +
  geom_point(data=big_df1 %>% filter(type=="real"),
             colour="orange") +
  ylab("Observed corpses") +
  xlab("Time, days after 2017-11-26")

ggsave("outputs/abc_posterior_predictive.pdf", g, width = 12, height = 8)