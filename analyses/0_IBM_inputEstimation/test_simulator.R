library(tidyverse)
library(yfmonkeyabc)
source("src/R/helper.R")
source("src/R/abc_helpers.R")


# test length of time it typically takes death to occur
compartments <- c("S",
                  "E0_0", "E0_1", "E0_2",
                  "E1_0", "E1_1", "E1_2",
                  "E2_0", "E2_1", "E2_2",
                  "I_0", "I_1", "I_2",
                  "D_unobs", "D_obs")
initial_states <- c(S=0, E0_0=1, E0_1=0, E0_2=0,
                    E1_0=0, E1_1=0, E1_2=0,
                    E2_0=0, E2_1=0, E2_2=0,
                    I_0=0, I_1=0, I_2=0,
                    D_unobs=0, D_obs=0)
# make p_obs really small
params <- c(beta=0.1, gamma=0.3, mu=0.6, p_obs=0.000001)
cnames <- c("time", compartments)

nreps <- 400
times_death <- vector(length = nreps)
max_stage_reacted <-  vector(length = nreps)
for(i in 1:nreps) {
  sims <- model_manystates(params, initial_states, tspan = seq(1, 50, 1),
                           tstart=0, tstop=100)[[2]] %>% 
    as.data.frame()
  colnames(sims) <- cnames
  times_death[i] <- sims$time[min(which(sims$D_unobs==1))]
  
  df <- sims %>% 
    pivot_longer(-time) %>% 
    mutate(exposed=str_detect(name, "E")) %>% 
    mutate()
  E1 <- sum(sims$E1_0+sims$E1_1+sims$E1_2) > 0
  E2 <- sum(sims$E2_0+sims$E2_1+sims$E2_2) > 0
  I <- sum(sims$I_0+sims$I_1+sims$I_2) > 0
  max_stage_reacted[i] <- if_else(I, 3, if_else(E2, 2, 1))
}
qplot(times_death)
qplot(max_stage_reacted)
