library(SimBIID)
library(tidyverse)
source("src/R/helper.R")

transitions <- c(
  "S -> beta * S * (I_0 + I_1 + I_2) -> E0_0",
  "E0_0 -> gamma * E0_0 -> E1_0",
  "E1_0 -> gamma * E1_0 -> E2_0",
  "E2_0 -> gamma * E2_0 -> I_0",
  "E0_1 -> gamma * E0_1 -> E1_1",
  "E1_1 -> gamma * E1_1 -> E2_1",
  "E2_1 -> gamma * E2_1 -> I_1",
  "E0_2 -> gamma * E0_2 -> E1_2",
  "E1_2 -> gamma * E1_2 -> E2_2",
  "E2_2 -> gamma * E2_2 -> I_2",
  "E0_0 -> mu * E0_0 -> E0_1",
  "E0_1 -> mu * E0_1 -> E0_2",
  "E0_2 -> mu * E0_2 -> D_unobs",
  "E1_0 -> mu * E1_0 -> E1_1",
  "E1_1 -> mu * E1_1 -> E1_2",
  "E1_2 -> mu * E1_2 -> D_unobs",
  "E2_0 -> mu * E2_0 -> E2_1",
  "E2_1 -> mu * E2_1 -> E2_2",
  "E2_2 -> mu * E2_2 -> D_unobs",
  "I_0 -> mu * I_0 -> I_1",
  "I_1 -> mu * I_1 -> I_2",
  "I_2 -> mu * I_2 -> D_unobs",
  "D_unobs -> p_obs * D_unobs -> D_obs"
)

compartments <- c("S",
                  "E0_0", "E0_1", "E0_2",
                  "E1_0", "E1_1", "E1_2",
                  "E2_0", "E2_1", "E2_2",
                  "I_0", "I_1", "I_2",
                  "D_unobs", "D_obs")
pars <- c("beta", "gamma", "mu", "p_obs")
model <- mparseRcpp(
  transitions = transitions, 
  compartments = compartments,
  pars = pars,
  tspan=TRUE
)
model_manystates <- compileRcpp(model)

simulator <- function(pars, initial_conditions, data, model) {
  pars <- unlist(pars)
  len_d <- length(data)
  D_obs <- data[1:(len_d - 1)]
  sims <- model(pars, initial_conditions, tspan = seq(1, 50, 1),
                tstart=0, tstop=100)[[2]]
  diff_1 <- sims[1:36, ncol(sims)] - D_obs
  D <- last(sims[, (ncol(sims) - 1)] + sims[, ncol(sims)])
  diff_2 <- D - last(data)
  diff <- c(diff_1, diff_2)
  total_rmse <- sqrt(mean(diff^2))
  total_rmse
}

f_R0_many_states <- function(beta, mu, gamma) {
  beta * 79 * 3 / mu * (gamma / (gamma + mu))^3
}

transitions <- c(
  "S -> beta * S * I -> E",
  "E -> gamma * E -> I",
  "E -> mu * E -> D_unobs",
  "I -> mu * I -> D_unobs",
  "D_unobs -> p_obs * D_unobs -> D_obs"
)

compartments <- c("S", "E", "I", "D_unobs", "D_obs")
pars <- c("beta", "gamma", "mu", "p_obs")
model <- mparseRcpp(
  transitions = transitions, 
  compartments = compartments,
  pars = pars,
  tspan=TRUE
)
model_seid <- compileRcpp(model)


transitions <- c(
  "S -> beta * S * I -> I",
  "I -> mu * I -> D_unobs",
  "D_unobs -> p_obs * D_unobs -> D_obs"
)

compartments <- c("S", "I", "D_unobs", "D_obs")
pars <- c("beta", "mu", "p_obs")
model <- mparseRcpp(
  transitions = transitions, 
  compartments = compartments,
  pars = pars,
  tspan=TRUE
)
model_sid <- compileRcpp(model)

prior_predictive_R0 <- function(ndraws, prior_sample, f_R0, parameter_names) {
  m_results <- matrix(nrow = ndraws, ncol = (1 + length(parameter_names)))
  for(i in 1:ndraws) {
    params <- prior_sample()
    R0 <- f_R0(params)
    m_results[i, ] <- c(as.vector(unlist(params)), R0)
  }
  colnames(m_results) <- c(parameter_names, "R0")
  m_results <- m_results %>%
    as.data.frame()
}

posterior_predictive_sims <- function(nreps, model, draws, initial_states,
                                      compartments) {
  cnames <- c("time", compartments)
  for(i in 1:nreps) {
    id <- sample(nrow(draws), 1)
    sims <- model(as.vector(unlist(draws[id, ])), initial_states, tspan = seq(1, 50, 1),
                      tstart=0, tstop=100)[[2]] %>% 
      as.data.frame()
    colnames(sims) <- cnames
    sims <- sims %>% mutate(iteration=i)
    if(i == 1)
      big_df <- sims
    else
      big_df <- big_df %>% bind_rows(sims)
  }
  big_df
}

plot_posterior_predictive_check <- function(big_df, abc_data) {
  big_df1 <- big_df %>% 
    mutate(type="simulated") %>% 
    bind_rows(abc_data$df_real %>% mutate(time=day, type="real"))
  
  ggplot(big_df1 %>% filter(type=="simulated"), aes(x=time, y=D_obs)) +
    geom_line(aes(group=as.factor(iteration)), alpha=0.5) +
    geom_point(data=big_df1 %>% filter(type=="real"),
               colour="orange") +
    ylab("Observed corpses") +
    xlab("Time, days after 2017-11-26")
}
