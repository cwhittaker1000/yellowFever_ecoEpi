library(rstan)
library(tidyverse)

diags <- readRDS("data/processed/monkey_death_fit_diagnostics.rds")
if(Reduce(`+`, diags) != 0)
  stop("Monkey death Stan model not converged.")

fit <- readRDS("data/processed/monkey_death_fit.rds")

# check p values
p_min <- mean(rstan::extract(fit, "is_min")[[1]])
p_truncated <- mean(rstan::extract(fit, "is_truncated")[[1]])
saveRDS(list(min=p_min, truncated=p_truncated),
        "data/processed/monkey_death_posterior_predictive_p.rds")

# plot actual and simulated data
df <- readRDS("data/processed/laemmert.rds") %>% 
  mutate(type="actual") %>% 
  filter(!is.na(days_post_infection))

days_sim <- rstan::extract(fit, "days_simulated")[[1]]
df_sim <- tibble(days_post_infection=days_sim, iteration=seq_along(days_sim)) %>% 
  mutate(type="simulated")

df_both <- df %>% 
  bind_rows(df_sim)

g <- ggplot(df_both, aes(x=days_post_infection, colour=type)) +
  geom_density() +
  geom_vline(xintercept = 10, linetype=2) +
  scale_x_continuous(limits=c(0, 15)) +
  scale_color_brewer("Data", palette="Dark2") +
  xlab("Days from infection to death") +
  ylab("Density")

ggsave("outputs/monkey_death_posterior_predictive.pdf", g,
       width=10, height = 6)
