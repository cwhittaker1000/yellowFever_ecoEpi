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
model <- compileRcpp(model)
cnames <- c("time", compartments)
inits <- c(S=79, E0_0=1, E0_1=0, E0_2=0,
           E1_0=0, E1_1=0, E1_2=0,
           E2_0=0, E2_1=0, E2_2=0,
           I_0=0, I_1=0, I_2=0,
           D_unobs=0, D_obs=0)

nreplicates <- 200
for(i in 1:nreplicates) {
  sims <- model(pars = c(beta=0.05, gamma=0.1 * 3, mu=0.2 * 3, p_obs=0.2),
              u = inits,
              tstart=0, tstop=100, tspan=seq(1, 36, 1))[[2]] %>% 
    as.data.frame()
  colnames(sims) <- cnames
  sims <- sims %>% mutate(iteration=i)
  if(i == 1)
    big_df <- sims
  else
    big_df <- big_df %>% bind_rows(sims)
}

big_df %>% 
  ggplot(aes(x=time, y=D_obs)) +
  geom_line(aes(group=as.factor(iteration)))

# try simple ABC
df_real <- readRDS("data/processed/epidemic_series_last_part.rds") %>% 
  mutate(day=seq_along(day))

## have stacked data like this because this is how ABCSMC needs it!
simSEIDDobs <- function(pars, data, tols, u) {
  len_d <- length(data)
  D_obs <- data[1:(len_d - 1)]
  sims <- model(pars, u, tspan = seq(1, 50, 1),
                tstart=0, tstop=100)[[2]]
  diff_1 <- sims[1:36, ncol(sims)] - D_obs
  D <- last(sims[, (ncol(sims) - 1)] + sims[, ncol(sims)])
  diff_2 <- D - last(data)
  diff <- c(diff_1, diff_2)
  total_rmse <- sqrt(mean(diff^2))
  if(total_rmse < tols[1])
    return(total_rmse)
  else
    return(NA)
}

f_R0 <- function(beta, mu, gamma) {
  beta * 79 * 3 / mu * (gamma / (gamma + mu))^3
}

days <- readRDS("data/processed/monkey_death_days.rds")
mu <- 1 / days
fit <- readRDS("data/processed/observation_time_fit.rds")
days_obs <- rstan::extract(fit, "days_simulated")[[1]]
p_obs <- 1 / days_obs
prior_sample <- function() {
  a_beta <- rnorm_truncated(0.03, 0.01)
  a_gamma <- rnorm_truncated(0.3, 0.05)
  a_mu <- sample(mu, 1) * 3
  a_p_obs <- sample(p_obs, 1)
  list(beta=a_beta, gamma=a_gamma, mu=a_mu, p_obs=a_p_obs)
}
threshold <- 3
ntries <- 10000
vbeta <- c()
vmu <- c()
vgamma <- c()
vpobs <- c()
vrmse <- c()
vR0_prior <- c()
for(i in 1:ntries) {
  params <- prior_sample()
  R0 <- f_R0(params$beta, params$mu, params$gamma)
  vR0_prior <- c(vR0_prior, R0)
  vals <- simSEIDDobs(pars = unlist(params),
              data=c(Dobs=df_real$D_obs, D_final_obs=80),
              tols=threshold,
              u = inits)
  if(!is.na(vals)) {
    vbeta <- c(vbeta, params$beta)
    vgamma <- c(vgamma, params$gamma)
    vmu <- c(vmu, params$mu)
    vpobs <- c(vpobs, params$p_obs)
    vrmse <- c(vrmse, vals)
  }
}
quantile(vR0_prior)

# uses slide 36 here: https://indico.ictp.it/event/7960/session/3/contribution/19/material/slides/0.pdf
draws <- tibble(beta=vbeta, gamma=vgamma, mu=vmu, p_obs=vpobs, rmse=vrmse) %>% 
  mutate(R0=f_R0(beta, mu, gamma))
qplot(draws$rmse, draws$R0)

nreps <- 400

for(i in 1:nreps) {
  id <- sample(nrow(draws), 1)
  parameters <- draws[id, ]
  df_sim <- model(pars = c(beta=parameters$beta,
                 gamma=parameters$gamma,
                 mu=parameters$mu,
                 p_obs=parameters$p_obs),
        u = inits,
        tstart=0, tstop=100, tspan=seq(1, 36, 1))[[2]] %>% 
    as.data.frame()
  colnames(df_sim) <- cnames
  df_sim <- df_sim %>% mutate(replicate=i)
  if(i == 1)
    big_df <- df_sim
  else
    big_df <- big_df %>% bind_rows(df_sim)
}

big_df1 <- big_df %>% 
  mutate(type="simulated") %>% 
  bind_rows(df_real %>% mutate(time=day, type="real"))

ggplot(big_df1 %>% filter(type=="simulated"), aes(x=time, y=D_obs)) +
  geom_line(aes(group=as.factor(replicate)), alpha=0.5) +
  geom_point(data=big_df1 %>% filter(type=="real"),
             colour="orange") +
  ylab("Observed corpses") +
  xlab("Time, days after 2017-11-26")
quantile(draws$R0)

# try ABC SMC -- can't get to work because I can't see how to make tols
# match data

simSEIDDobs(pars = c(beta=0.02, gamma=0.2, mu=0.2, p_obs=0.2),
            data=c(Dobs=df_real$D_obs, D_final_obs=80),
            tols=200,
            u = c(S=79, E=0, I=1, D_unobs=0, D_obs=0))

priors <- data.frame(parnames=c("beta", "gamma", "mu", "p_obs"),
                     dist=c("norm", "norm", "norm", "norm"),
                     stringsAsFactors = FALSE)
priors$p1 <- c(0.01, 0.1, 0.2, 0.2)
priors$p2 <- c(0.001, 0.01, 0.02, 0.02)
data=c(Dobs=df_real$D_obs, D_final_obs=80)
inits = c(S=79, E=0, I=1, D_unobs=0, D_obs=0)
tols <- 50
post <- ABCSMC(
  x = data, 
  priors = priors, 
  func = simSEIDDobs, 
  u = iniStates, 
  tols = c(100, rep(0, 36)), 
  ptol = 0.2, 
  ngen = 2, 
  npart = 50,
  model = model
)
