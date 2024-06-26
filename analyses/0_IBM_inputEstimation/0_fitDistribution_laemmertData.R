## Loading required libraries
library(tidyverse); library(rstan)

# Recreating table 1 from "The susceptibility of Howler monkeys to yellow fever virus" Laemmert & Klumm, 1950.

## Interpret (e.g.) 6th day after infection as being 5 days post-infection; assume monkey 4 dies on day 10
days_death_post_infection <- c(5, 8, 2, 10, 6, 4, 3, 4, 3, 9)
## Interpret - "trace" doesn't count; if detectable on 1st day, then 0.5; count 6th day after infection as being 5 days post-infection
days_virus_post_infection <- c(2, 2, 1, 2, 2, 1, 1, 1, 0.5, 2)

df <- tibble(id = seq_along(days_virus_post_infection),
             days_virus_post_infection = days_virus_post_infection,
             days_death_post_infection = days_death_post_infection,
             days_death_post_virus = days_death_post_infection - days_virus_post_infection)  %>% 
  mutate(death_truncated = if_else(is.na(days_death_post_infection), 1, 0)) %>% 
  mutate(days_death_post_infection = if_else(is.na(days_death_post_infection), -99, days_death_post_infection))

## Fitting the time between infection and death
model <- stan_model("analyses/0_IBM_inputEstimation/models/laemmert_monkey_death.stan")
data_stan <- list(N = nrow(df),
                  days = df$days_death_post_infection,
                  truncated = df$death_truncated,
                  days_truncated = 10,
                  a_1 = 1,
                  a_2 = 0.5,
                  b_1 = 1,
                  b_2 = 2)
fit <- sampling(model, data=data_stan, iter=2000, chains=4)
hist(rstan::extract(fit, "days_simulated")[[1]], breaks = 50)
summary(fit)
saveRDS(fit, "analyses/0_IBM_inputEstimation/outputs/infection_deathDist_stanFit.rds")

## Comparing fit to empirical
days_sim <- rstan::extract(fit, "days_simulated")[[1]]
df_sim <- tibble(days_death_post_infection = days_sim, iteration = seq_along(days_sim)) %>% 
  mutate(type = "simulated")
df_real <- tibble(days_death_post_infection = df$days_death_post_infection[!is.na(days_death_post_infection)], 
                  iteration = 1:10,
                  type = "actual")
df_both <- df_sim %>% 
  bind_rows(df_real)
a <- ggplot(df_both, aes(x=days_death_post_infection, colour=type)) +
  geom_density() +
  geom_vline(xintercept = 10, linetype=2) +
  scale_x_continuous(limits = c(0, 15)) +
  scale_color_brewer("Data", palette = "Dark2") +
  xlab("Days from exposure to death") +
  ylab("Density")
mean(df$days_death_post_infection)
median(df$days_death_post_infection)
mean(days_sim)
median(days_sim)

## Fitting the time between infection and death
data_stan <- list(N = nrow(df),
                  days = df$days_virus_post_infection,
                  truncated = df$death_truncated,
                  days_truncated = 10,
                  a_1 = 1,
                  a_2 = 0.5,
                  b_1 = 1,
                  b_2 = 2)
fit <- sampling(model, data=data_stan, iter=2000, chains=4)
hist(rstan::extract(fit, "days_simulated")[[1]], breaks = 50)
summary(fit)
saveRDS(fit, "analyses/0_IBM_inputEstimation/outputs/exposure_infectiousDist_stanFit.rds")

## Comparing fit to empirical
days_sim <- rstan::extract(fit, "days_simulated")[[1]]
df_sim <- tibble(days_virus_post_infection = days_sim, iteration = seq_along(days_sim)) %>% 
  mutate(type = "simulated")
df_real <- tibble(days_virus_post_infection = df$days_virus_post_infection[!is.na(days_virus_post_infection)], 
                  iteration = 1:10,
                  type = "actual")
df_both <- df_sim %>% 
  bind_rows(df_real)
b <- ggplot(df_both, aes(x=days_virus_post_infection, colour=type)) +
  geom_density() +
  geom_vline(xintercept = 10, linetype=2) +
  scale_x_continuous(limits = c(0, 15)) +
  scale_color_brewer("Data", palette = "Dark2") +
  xlab("Days from exposure to infectious") +
  ylab("Density")
mean(df$days_virus_post_infection)
median(df$days_virus_post_infection)
mean(days_sim)
median(days_sim)

## Fitting the time between infectiousness (virus detectable) and death
data_stan <- list(N = nrow(df),
                  days = df$days_death_post_virus,
                  truncated = df$death_truncated,
                  days_truncated = 10,
                  a_1 = 1,
                  a_2 = 0.5,
                  b_1 = 1,
                  b_2 = 2)
fit <- sampling(model, data=data_stan, iter=2000, chains=4)
hist(rstan::extract(fit, "days_simulated")[[1]], breaks = 50)
summary(fit)
saveRDS(fit, "analyses/0_IBM_inputEstimation/outputs/infectious_deathDist_stanFit.rds")

## Comparing fit to empirical
days_sim <- rstan::extract(fit, "days_simulated")[[1]]
df_sim <- tibble(days_death_post_virus = days_sim, iteration = seq_along(days_sim)) %>% 
  mutate(type = "simulated")
df_real <- tibble(days_death_post_virus = df$days_death_post_virus, 
                  iteration = 1:10,
                  type = "actual")
df_both <- df_sim %>% 
  bind_rows(df_real)
c <- ggplot(df_both, aes(x=days_death_post_virus, colour=type)) +
  geom_density() +
  geom_vline(xintercept = 10, linetype=2) +
  scale_x_continuous(limits = c(0, 15)) +
  scale_color_brewer("Data", palette = "Dark2") +
  xlab("Days from infectious to death") +
  ylab("Density")
mean(df$days_death_post_virus)
median(df$days_death_post_virus)
mean(days_sim)
median(days_sim)

cowplot::plot_grid(a, b, c, nrow = 1)