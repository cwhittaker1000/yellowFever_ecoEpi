# Load required libraries
library(individual); library(dplyr); library(EasyABC); library(tidyverse); library(parallel)

# Sourcing required functions
source("analyses/1_IBM_estimationR0/functions/run_simulation.R")
source("analyses/1_IBM_estimationR0/functions/abc_functions.R")

# Generating synthetic data

## Defining parameters
R0 <- 3
N <- 86
infectious_period_gamma_shape <- 8
infectious_period_gamma_rate <- 2
gamma <- 1 / (infectious_period_gamma_shape / infectious_period_gamma_rate)
beta_synth <- R0 * gamma / N
past_length_max <- 20
parameters_list <- list(seed = 100268,
                        steps = 400, 
                        dt = 0.2,
                        overall_run_length = 400, 
                        N = N, 
                        beta = beta_synth,
                        initial_infections = 1,
                        death_obs_prop = 1, 
                        infectious_period_gamma_shape = infectious_period_gamma_shape, 
                        infectious_period_gamma_rate = infectious_period_gamma_rate,
                        latent_period_gamma_shape = 8, 
                        latent_period_gamma_rate = 2,
                        death_observation_gamma_shape = 8, 
                        death_observation_gamma_rate = 2,
                        past_length = past_length_max, 
                        past_weightings_vector = rev(dgamma(1:past_length_max, shape = 12 * 0.5, rate = 0.5)),
                        lagged_I_input = NULL)

## Running model and summarising output
synthetic_data <- run_simulation(parameters_list)
synthetic_data$result$time <- floor(synthetic_data$result$timestep * parameters_list$dt)
observed_incidence <- synthetic_data$result %>%
  mutate(incidence = c(0, diff(Dobs_count))) %>%
  group_by(time) %>%
  summarise(daily_incidence = sum(incidence)) %>%
  mutate(cumulative_incidence = cumsum(daily_incidence))
observed_data <- observed_incidence$daily_incidence
plot(observed_data)
plot(observed_incidence$cumulative_incidence)

# Test running ABC with parallel processing and RMSE
set.seed(123)
observed_data <- observed_incidence$cumulative_incidence
tolerance <- 15
n <- 250
accepted_parameters <- abc_parallel_rmse(
  observed_data = observed_data,
  prior_sampler = prior_sampler,
  model_simulate = model_simulate,
  parameters_list = parameters_list, ## note that the argument input has to be called "parameters_list" so it can be exported...
  rmse = rmse,
  tolerance = tolerance,  # Set a tolerance level, adjusted for RMSE scale
  n = n  # Number of accepted samples we want
)
hist(unlist(accepted_parameters$accepted_parameters) * N / gamma)
mean(unlist(accepted_parameters$accepted_parameters) * N / gamma)
median(unlist(accepted_parameters$accepted_parameters) * N / gamma)
quantile(unlist(accepted_parameters$accepted_parameters) * N / gamma, probs = c(0.05, 0.5, 0.95))

x <- matrix(data = unlist(accepted_parameters$model_output), byrow = TRUE,
            nrow = n, ncol = length(accepted_parameters$model_output[[1]]))
colnames(x) <- seq(0:(length(accepted_parameters$model_output[[1]]) - 1))
y <- data.frame(x, iteration = 1:n) %>%
  pivot_longer(cols = -iteration, names_to = "time", values_to = "cumsum") %>%
  mutate(time = as.numeric(str_remove(time, "X")))

ggplot() +
  geom_line(data = y, aes(x = time, y = cumsum, group = iteration)) +
  geom_line(data = observed_incidence, aes(x = time, y = cumulative_incidence), colour = "red")

