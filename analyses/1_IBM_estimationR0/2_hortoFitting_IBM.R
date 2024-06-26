# Load required libraries
library(individual); library(dplyr); library(EasyABC); library(tidyverse); library(parallel)

# Sourcing required functions
### -> Do I want to be able to model importations explicitly? (Or at least spread out initial infections over multiple days)
###    Might require model alterations/elaborations
### -> Need tMRCA and importation frequency / initial infections from the genomic analyses
### -> Do I want to be conditioning my simulations on there being an outbreak? Worried that we're biasing
###    estimates towards higher R0 values because of extinctions @ low R0.
### -> Is frequency or density dependent transmission more appropriate here do we think? Surely frequency
###    given population is entirely wiped out.
### -> Try outputing d_obs and d_unobs - seeing some v slight weirdness where the ABC results tend to end
###    at a lower cumulative incidence than the empirical data (only slightly) - wondering whether
###    that's indicative of an issue somewhere?
source("analyses/1_IBM_estimationR0/functions/run_simulation.R")
source("analyses/1_IBM_estimationR0/functions/abc_functions.R")

# Loading in and processing Horto/PEAL data for model fitting
horto_df <- readRDS("analyses/1_IBM_estimationR0/data/processed_HortoData.rds") %>%
  filter(!is.na(zone_peal)) %>%
  filter(final_yfv_result != "negative")
epi_curve <- incidence::incidence(horto_df$date_collection)
plot(epi_curve)

# Generating cumulative incidence data and cutting off first 4 infections
start_date <- as.Date("2017-12-01")
horto_df_fitting <- horto_df %>%
  filter(date_collection > start_date) %>%
  group_by(date_collection) %>%
  summarise(count = n()) %>%
  complete(date_collection = seq.Date(start_date, 
                                      as.Date("2018-01-08"), 
                                      by = "days"),
           fill = list(count = 0)) %>%
  mutate(cumsum = cumsum(count))
plot(horto_df_fitting$date_collection, horto_df_fitting$count)
plot(horto_df_fitting$date_collection, horto_df_fitting$cumsum)
horto_df_fitting$time <- 1:nrow(horto_df_fitting)

# Setting up parameters list

## Loading in fitted parameters for model
death_observation_fit <- readRDS("analyses/0_IBM_inputEstimation/outputs/deathObservation_distExp_stanFit.rds")
death_observation_gamma_rate <- mean(rstan::extract(death_observation_fit, "a")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
latent_period_fit <- readRDS("analyses/0_IBM_inputEstimation/outputs/exposure_infectiousDist_stanFit.rds")
latent_period_gamma_shape <- mean(rstan::extract(latent_period_fit, "a")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
latent_period_gamma_rate <- mean(rstan::extract(latent_period_fit, "b")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
infectious_period_fit <- readRDS("analyses/0_IBM_inputEstimation/outputs/infectious_deathDist_stanFit.rds")
infectious_period_gamma_shape <- mean(rstan::extract(infectious_period_fit, "a")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
infectious_period_gamma_rate <- mean(rstan::extract(infectious_period_fit, "b")[[1]]) # note exp used here and gamma below, but if shape set to 1, then is an exponential
EIP_gamma_fit <- readRDS("analyses/0_IBM_inputEstimation/outputs/EIP_adultMice_gammaParams.rds")
EIP_gamma_shape <- EIP_gamma_fit$gamma_a
EIP_gamma_rate <- EIP_gamma_fit$gamma_b

N <- 86 - 3 # (3 negative monkeys - assume the rest not in database killed by yellow fever)
past_length_max <- 15
parameters_list <- list(seed = 15,
                        steps = nrow(horto_df_fitting) / 0.2, 
                        dt = 0.2,
                        overall_run_length = nrow(horto_df_fitting) / 0.2, 
                        N = N, 
                        beta = 0.1,
                        initial_infections = 3,
                        death_obs_prop = sum(horto_df_fitting$count) / N, 
                        latent_period_gamma_shape = latent_period_gamma_shape, 
                        latent_period_gamma_rate = latent_period_gamma_rate,
                        infectious_period_gamma_shape = infectious_period_gamma_shape, 
                        infectious_period_gamma_rate = infectious_period_gamma_rate,
                        death_observation_gamma_shape = 1, 
                        death_observation_gamma_rate = death_observation_gamma_rate,
                        past_length = past_length_max, 
                        past_weightings_vector = rev(dgamma(1:past_length_max, shape = EIP_gamma_shape, rate = EIP_gamma_rate)), # from re-analysis of Bates and Roca-Garcia 1946
                        ### does this need to be 0:past_length_max and do I need to be changing the structure of lagged vector in model to account for this?
                        lagged_I_input = NULL)

# Model fitting 
set.seed(123)
observed_data <- horto_df_fitting$cumsum
tolerance <- 10
n <- 100
accepted_parameters <- abc_parallel_rmse(
  observed_data = observed_data,
  prior_sampler = prior_sampler,
  model_simulate = model_simulate,
  parameters_list = parameters_list, ## note that the argument input has to be called "parameters_list" so it can be exported...
  rmse = rmse,
  tolerance = tolerance,  # Set a tolerance level, adjusted for RMSE scale
  n = n  # Number of accepted samples we want
)
gamma <- 1 / (parameters_list$infectious_period_gamma_shape / parameters_list$infectious_period_gamma_rate)

hist(unlist(accepted_parameters$accepted_parameters) * N / gamma)
mean(unlist(accepted_parameters$accepted_parameters) * N / gamma)
median(unlist(accepted_parameters$accepted_parameters) * N / gamma)

x <- matrix(data = unlist(accepted_parameters$model_output), byrow = TRUE,
            nrow = n, ncol = length(accepted_parameters$model_output[[1]]))
colnames(x) <- seq(0:(length(accepted_parameters$model_output[[1]]) - 1))
y <- data.frame(x, iteration = 1:n) %>%
  pivot_longer(cols = -iteration, names_to = "time", values_to = "cumsum") %>%
  mutate(time = as.numeric(str_remove(time, "X")))
ggplot() +
  geom_line(data = y, aes(x = time, y = cumsum, group = iteration)) +
  geom_line(data = horto_df_fitting, aes(x = time, y = cumsum), colour = "red")

## Fits with positives only 
start_date <- as.Date("2017-12-01")
horto_df_fitting_pos_only <- horto_df %>%
  filter(date_collection > start_date) %>%
  filter(final_yfv_result == "positive") %>%
  group_by(date_collection) %>%
  summarise(count = n()) %>%
  complete(date_collection = seq.Date(start_date, 
                                      as.Date("2018-01-08"), 
                                      by = "days"),
           fill = list(count = 0)) %>%
  mutate(cumsum = cumsum(count))
plot(horto_df_fitting_pos_only$date_collection, horto_df_fitting_pos_only$count)
plot(horto_df_fitting_pos_only$date_collection, horto_df_fitting_pos_only$cumsum)
horto_df_fitting_pos_only$time <- 1:nrow(horto_df_fitting_pos_only)

## Parameters list
N <- 86 - 3 # (3 negative monkeys - assume the rest not in database killed by yellow fever)
past_length_max <- 15
parameters_list <- list(seed = 15,
                        steps = nrow(horto_df_fitting_pos_only) / 0.2, 
                        dt = 0.2,
                        overall_run_length = nrow(horto_df_fitting_pos_only) / 0.2, 
                        N = N, 
                        beta = 0.1,
                        initial_infections = 3,
                        death_obs_prop = sum(horto_df_fitting_pos_only$count) / N, 
                        latent_period_gamma_shape = latent_period_gamma_shape, 
                        latent_period_gamma_rate = latent_period_gamma_rate,
                        infectious_period_gamma_shape = infectious_period_gamma_shape, 
                        infectious_period_gamma_rate = infectious_period_gamma_rate,
                        death_observation_gamma_shape = 1, 
                        death_observation_gamma_rate = death_observation_gamma_rate,
                        past_length = past_length_max, 
                        past_weightings_vector = rev(dgamma(1:past_length_max, shape = 6.9 * 0.5, rate = 0.5)), # from Kraemer et al (Aedes)
                        # past_weightings_vector = rev(dgamma(1:past_length_max, shape = 10 * 0.5, rate = 0.5)),
                        lagged_I_input = NULL)

# Model fitting 
set.seed(123)
observed_data <- horto_df_fitting_pos_only$cumsum
tolerance <- 6
n <- 100
accepted_parameters_pos_only <- abc_parallel_rmse(
  observed_data = observed_data,
  prior_sampler = prior_sampler,
  model_simulate = model_simulate,
  parameters_list = parameters_list, ## note that the argument input has to be called "parameters_list" so it can be exported...
  rmse = rmse,
  tolerance = tolerance,  # Set a tolerance level, adjusted for RMSE scale
  n = n  # Number of accepted samples we want
)
gamma <- 1 / (parameters_list$infectious_period_gamma_shape / parameters_list$infectious_period_gamma_rate)

hist(unlist(accepted_parameters_pos_only$accepted_parameters) * N / gamma)
mean(unlist(accepted_parameters_pos_only$accepted_parameters) * N / gamma)
median(unlist(accepted_parameters_pos_only$accepted_parameters) * N / gamma)

x <- matrix(data = unlist(accepted_parameters_pos_only$model_output), byrow = TRUE,
            nrow = n, ncol = length(accepted_parameters_pos_only$model_output[[1]]))
colnames(x) <- seq(0:(length(accepted_parameters_pos_only$model_output[[1]]) - 1))
y <- data.frame(x, iteration = 1:n) %>%
  pivot_longer(cols = -iteration, names_to = "time", values_to = "cumsum") %>%
  mutate(time = as.numeric(str_remove(time, "X")))
ggplot() +
  geom_line(data = y, aes(x = time, y = cumsum, group = iteration)) +
  geom_line(data = horto_df_fitting_pos_only, aes(x = time, y = cumsum), colour = "red")
