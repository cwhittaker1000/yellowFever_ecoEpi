## calculate rmse between observed and simulated
rmse <- function(observed, simulated) {
  sqrt(mean((observed - simulated)^2))
}

## draw from the prior distribution
prior_sampler <- function(parameters_list) {
  R0 <- runif(n = 1, min = 1, max = 10) 
  gamma <- 1 / (parameters_list$infectious_period_gamma_shape / parameters_list$infectious_period_gamma_rate)
  beta_sim <- R0 * gamma / parameters_list$N
  return(beta_sim)
}

## runs model with parameter inputs and returns quantity to be compared (cumulative incidence)
model_simulate <- function(parameters_list, seed) {
  
  ## Running model
  parameters_list$seed <- seed
  model_output <- run_simulation(parameters_list = parameters_list)
  
  ## Calculating cumulative incidence
  model_output$result$time <- floor(model_output$result$timestep * 0.2)
  model_output2 <- model_output$result %>%
    mutate(incidence = c(0, diff(Dobs_count))) %>%
    group_by(time) %>%
    summarise(daily_incidence = sum(incidence)) %>%
    mutate(cumulative_incidence = cumsum(daily_incidence))
  return(model_output2$cumulative_incidence)
}

# Parallelised ABC function
abc_parallel_rmse <- function(observed_data, prior_sampler, model_simulate, parameters_list, rmse, tolerance, n) {
  
  parameters_list <- parameters_list
  
  # Function to try a set of parameters and return NA if not accepted, param if accepted
  try_parameter <- function(param, parameters_list) {
    temp_parameters_list <- parameters_list
    temp_parameters_list$beta <- param
    simulated_data <- model_simulate(parameters_list = temp_parameters_list, 
                                     seed = rnbinom(n = 1, mu = 100000, size = 2))
    if(rmse(observed_data, simulated_data) < tolerance) {
      return(list(param = param, model_output = simulated_data))
    } else {
      return(list(param = NA, model_output = NA))  # Return NA if the parameter is not accepted
    }
  }
  
  # Detect the number of available cores and set up cluster with required inputs
  no_cores <- detectCores() - 3
  cl <- makeCluster(no_cores)
  clusterExport(cl, c("prior_sampler", "rmse", "observed_data", "tolerance", "run_simulation", "model_simulate", "parameters_list"))
  clusterEvalQ(cl, library(individual))
  clusterEvalQ(cl, library(dplyr))
  
  # Store accepted parameters, loop until enough parameters are accepted
  accepted_parameters <- numeric(0)
  model_output_list <- list()
  while(length(accepted_parameters) < n) {

    # Run 100 different parameter combinations    
    params <- replicate(100, prior_sampler(parameters_list = parameters_list))
    trial_results <- parLapply(cl, as.list(params), try_parameter, parameters_list = parameters_list)
    parameter_outputs <- unlist(lapply(trial_results, "[[", 1))

    # Filter out NAs and append to accepted parameters
    accepted_index <- which(!is.na(parameter_outputs))
    accepted_parameters <- c(accepted_parameters, unlist(Filter(Negate(is.na), parameter_outputs)))
    model_output_list <- c(model_output_list, lapply(trial_results[accepted_index], "[[", 2))
    
    # If more parameters are accepted than needed, truncate the list
    if(length(accepted_parameters) > n) {
      accepted_parameters <- accepted_parameters[1:n]
      model_output_list <- model_output_list[1:n]
    }
    print(length(accepted_parameters))
  }
  stopCluster(cl)
  
  return(list(accepted_parameters = accepted_parameters,
              model_output = model_output_list))
}