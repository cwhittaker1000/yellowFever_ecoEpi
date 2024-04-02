library(rstan)
library(tidyverse)

diags <- readRDS("data/processed/observation_time_fit_diagnostics.rds")
if(Reduce(`+`, diags) != 0)
  stop("Observation time Stan model not converged.")

fit <- readRDS("data/processed/observation_time_fit.rds")

# check p values -- not sure how useful these are: i) 0.5 days is
# probably the least time an observer would say and ii) lots of
# corpses weren't observed so the max in the data is likely a massive
# understatement. This is probably why the p values are near 1
p_min <- mean(rstan::extract(fit, "is_min")[[1]])
p_max <- mean(rstan::extract(fit, "is_max")[[1]])
saveRDS(list(min=p_min, max=p_max),
        "data/processed/observation_time_posterior_predictive_p.rds")

# plot actual and simulated data
df <- readRDS("data/processed/cleaned_data.rds") %>% 
  filter(days_since_death_interpolated != 1) %>% 
  mutate(type="actual")

days_sim <- rstan::extract(fit, "days_simulated")[[1]]
df_sim <- tibble(days_since_death=days_sim, iteration=seq_along(days_sim)) %>% 
  mutate(type="simulated")

df_both <- df %>% 
  bind_rows(df_sim)

g <- ggplot(df_both %>% filter(type=="simulated"),
            aes(x=days_since_death)) +
  geom_histogram(data=df_both %>% filter(type=="actual"),
                 aes(y = ..density..), fill="blue") +
  geom_density() +
  scale_x_continuous(limits=c(0, 15)) +
  scale_color_brewer("Data", palette="Dark2") +
  xlab("Days from death to corpse observation") +
  ylab("Density")

ggsave("outputs/observation_time_posterior_predictive.pdf", g,
       width=10, height = 6)
