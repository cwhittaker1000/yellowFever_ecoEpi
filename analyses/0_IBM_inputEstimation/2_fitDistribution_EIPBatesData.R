## Loading required libraries
library(tidyverse); library(rstan); library(loo)

## Sorucing helper DIC function
calculate_DIC <- function(stan_fit) {
  deviance <- -2 * rstan::extract(stan_fit, "lp__")[[1]]
  mean_deviance <- mean(deviance)
  pD_gelman <- 0.5 * var(deviance)
  DIC <- mean_deviance + pD_gelman
  return(DIC)
}

## Recreating data from "The Development of the Virus of Yellow Fever in Haemagogus Mosquitoes"
##                       by Bates and Roca-Garcia 1946
## Table 2 @ 30 degrees - ignoring the two lots at Day 14 where nothing was infected
## Bit confused about the bit of the "first transmission at Day 13".
degrees30_df_adults <- tibble(days = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20),
                              total= c(28,26,27,29,27,26,30,30,29,30,30, 29, 17, 30, 18, 12),    
                              died = c(15,6, 1, 6 , 6, 9, 8, 2,10,15,14, 14, 9, 19, 15, 7),
                              mod =  c(1,1, 1, 6 , 6, 9, 8, 2,10,15,14, 14, 9, 19, 15, 7))
degrees30_df_adults$perc <- degrees30_df_adults$mod/degrees30_df_adults$total

degrees30_df_babies <- tibble(days = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20),
                              total= c(23,25,15,25,24,25,23,25,25,25,22, 25, 13, 24, 14, 10),
                              died = c(23,14,5, 18,18,24,17,18,25,25,21, 25, 13, 23, 14, 9),
                              mod  = c(8, 8, 5, 18,18,24,17,18,25,25,21, 25, 13, 23, 14, 9))
degrees30_df_babies$perc <- degrees30_df_babies$mod/degrees30_df_babies$total

### Fitting gamma to babies and adult data
model_gamma <- stan_model("analyses/0_IBM_inputEstimation/models/EIP_gamma.stan")

## Babies fit
data_stan_babies <- list(N = length(degrees30_df_babies$days),
                         day = degrees30_df_babies$days,
                         infected = degrees30_df_babies$total,
                         died = degrees30_df_babies$mod,
                         a_1 = 0.1,
                         a_2 = 10,
                         b_1 = 0.1,
                         b_2 = 10,
                         min_p_death_prior_mean = 0.33,
                         min_p_death_prior_sd = 0.02,
                         max_p_death_prior_mean = 1,
                         max_p_death_prior_sd = 0.1)
fit_babies <- sampling(model_gamma, data = data_stan_babies, iter = 2000, chains = 4)
summary(fit_babies)

## Adults fit
data_stan_adults <- list(N = length(degrees30_df_adults$days),
                         day = degrees30_df_adults$days,
                         infected = degrees30_df_adults$total,
                         died = degrees30_df_adults$mod,
                         a_1 = 0.1,
                         a_2 = 10,
                         b_1 = 0.1,
                         b_2 = 10,
                         min_p_death_prior_mean = 0.04,
                         min_p_death_prior_sd = 0.02,
                         max_p_death_prior_mean = 0.6,
                         max_p_death_prior_sd = 0.1)
fit_adults <- sampling(model_gamma, data = data_stan_adults, iter = 2000, chains = 4)
summary(fit_adults)

## Plotting
par(mfrow = c(2, 2))
hist(rstan::extract(fit_babies, "days_simulated")[[1]], breaks = 50)
mean(rstan::extract(fit_babies, "days_simulated")[[1]])
pmin_babies <- mean(rstan::extract(fit_babies, "p_death_min")[[1]])
pmax_babies <- mean(rstan::extract(fit_babies, "p_death_max")[[1]])
degrees30_df_babies$perc_mod <- (degrees30_df_babies$perc - pmin_babies) / (pmax_babies - pmin_babies)
plot(degrees30_df_babies$days, degrees30_df_babies$perc_mod, ylim = c(0, 1), type = "l")
x <- seq(0, 20, by = 0.01)
modelled_a <- mean(rstan::extract(fit_babies, "a")[[1]])
modelled_b <- mean(rstan::extract(fit_babies, "b")[[1]])
cdf_values <- pgamma(x, shape = modelled_a, rate = modelled_b)
lines(x, cdf_values, type = "l", col = "red")

hist(rstan::extract(fit_adults, "days_simulated")[[1]], breaks = 50)
mean(rstan::extract(fit_adults, "days_simulated")[[1]])
pmin_adults <- mean(rstan::extract(fit_adults, "p_death_min")[[1]])
pmax_adults <- mean(rstan::extract(fit_adults, "p_death_max")[[1]])
degrees30_df_adults$perc_mod <- (degrees30_df_adults$perc - pmin_adults) / (pmax_adults - pmin_adults)
plot(degrees30_df_adults$days, degrees30_df_adults$perc_mod, ylim = c(0, 1.5), type = "l")
x <- seq(0, 20, by = 0.01)
modelled_a <- mean(rstan::extract(fit_adults, "a")[[1]])
modelled_b <- mean(rstan::extract(fit_adults, "b")[[1]])
cdf_values <- pgamma(x, shape = modelled_a, rate = modelled_b)
lines(x, cdf_values, type = "l", col = "blue")




#############

# Stan model running
model <- stan_model("analyses/0_IBM_inputEstimation/models/EIP_Offsetgamma.stan")
data_stan <- list(N = length(degrees30_df_babies$days),
                  day = degrees30_df_babies$days,
                  infected = degrees30_df_babies$total,
                  died = degrees30_df_babies$mod,
                  a_1 = 0.01,
                  a_2 = 10,
                  b_1 = 0.01,
                  b_2 = 10,
                  offset_mean = 2)
fit <- sampling(model, data=data_stan, iter=2000, chains=4)
summary(fit)
print(fit)
pairs(fit)
hist(rstan::extract(fit, "a")[[1]], breaks = 50)
hist(rstan::extract(fit, "b")[[1]], breaks = 50)
hist(rstan::extract(fit, "days_simulated")[[1]], breaks = 100, xlim = c(0, 20))

plot(degrees30_df_babies$days, degrees30_df_babies$perc_mod, ylim = c(0, 1), type = "l")
x <- seq(0, 20, by = 0.01)
modelled_a <- mean(rstan::extract(fit, "a")[[1]])
modelled_b <- mean(rstan::extract(fit, "b")[[1]])
modelled_offset <- mean(rstan::extract(fit, "offset")[[1]])
cdf_values <- pgamma(x, shape = modelled_a, rate = modelled_b)
lines(x + modelled_offset, cdf_values, type = "l", col = "red")

# Stan model running
model2 <- stan_model("analyses/0_IBM_inputEstimation/models/EIP_weibull.stan")
data_stan <- list(N = length(degrees30_df_babies$days),
                  day = degrees30_df_babies$days,
                  infected = degrees30_df_babies$total,
                  died = degrees30_df_babies$mod,
                  a_1 = 1,
                  a_2 = 10,
                  b_1 = 1,
                  b_2 = 10)
fit2 <- sampling(model2, data=data_stan, iter=2000, chains=1)
summary(fit2)
hist(rstan::extract(fit2, "a")[[1]], breaks = 50)
hist(rstan::extract(fit2, "b")[[1]], breaks = 50)
hist(rstan::extract(fit2, "days_simulated")[[1]], breaks = 50)

plot(degrees30_df_babies$days, degrees30_df_babies$perc_mod, ylim = c(0, 1), type = "l")
x <- seq(0, 20, by = 0.01)
modelled_a <- mean(rstan::extract(fit2, "a")[[1]])
modelled_b <- mean(rstan::extract(fit2, "b")[[1]])
cdf_values <- pweibull(x, shape = modelled_a, scale = modelled_b)
lines(x, cdf_values, type = "l", col = "red")

