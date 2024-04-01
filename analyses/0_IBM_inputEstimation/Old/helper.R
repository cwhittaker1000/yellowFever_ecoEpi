library(mvtnorm)

rnorm_truncated <- function(mu, sigma) {
  val <- rnorm(1, mu, sigma)
  while(val < 0)
    val <- rnorm(1, mu, sigma)
  val
}

check_diagnostics <- function(filename) {
  fit <- readRDS(filename)
  fit_diagnostics <- summarise_draws(fit)
  rhat_above_1.01 <- sum(fit_diagnostics$rhat > 1.01, na.rm = T)
  ess_bulk_below_400 <- sum(fit_diagnostics$ess_bulk < 400, na.rm = T)
  ess_tail_below_400 <- sum(fit_diagnostics$ess_tail < 400, na.rm = T)
  list(rhat=rhat_above_1.01,
       ess_bulk=ess_bulk_below_400,
       ess_tail=ess_tail_below_400)
}

rmvrnorm2D <- function(n, mux, muy, sigmax, sigmay, rho){
  rmvnorm(n, c(mux, muy),
          matrix(c(sigmax^2, sigmax * sigmay * rho,
                   sigmax * sigmay * rho, sigmay^2),
                 ncol = 2))
}
rmvrnorm2D_truncated <- function(n, mux, muy, sigmax, sigmay, rho) {
  ndraws <- 0
  beta <- c()
  mu <- c()
  while(ndraws < n) {
    draw <- rmvrnorm2D(1, mux, muy, sigmax, sigmay, rho)
    while(draw[1, 1] < 0 || draw[1, 2] < 0) {
      draw <- rmvrnorm2D(1, mux, muy, sigmax, sigmay, rho)
    }
    beta <- c(beta, draw[1, 1])
    mu <- c(mu, draw[1, 2])
    ndraws <- ndraws + 1
  }
  tibble(beta, mu)
}

mvt_prior_sample <- function(n, beta_bar, mu_bar, beta_sigma, mu_sigma, rho,
                             p_obs_lower, p_obs_upper, max_monkeys=84) {
  vals <- rmvrnorm2D_truncated(n, beta_bar, mu_bar, beta_sigma, mu_sigma, rho)
  p_obs <- runif(n, p_obs_lower, p_obs_upper)
  tibble(beta=vals$beta, mu=vals$mu, p_obs=p_obs) %>% 
    mutate(R0=beta / mu * max_monkeys)
}

abc_prior_predictive <- function(params,
                                 days,
                                 filename_end,
                                 ndraws=1000) {
  beta_bar <- params$beta_bar
  mu_bar <- params$mu_bar
  beta_sigma <- params$beta_sigma
  mu_sigma <- params$mu_sigma
  rho <- params$rho
  p_obs_lower <- params$p_obs_lower
  p_obs_upper <- params$p_obs_upper
  
  draws <- mvt_prior_sample(ndraws, beta_bar, mu_bar, beta_sigma, mu_sigma, rho,
                            p_obs_lower, p_obs_upper)
  draws
}

abc_run <- function(abc_params,
                    prior_params,
                    abc_data,
                    filename_end,
                    ndraws=1000) {
  
  prior_sample <- function() {
    vals <- mvt_prior_sample(1, beta_bar, mu_bar, beta_sigma, mu_sigma, rho,
                     p_obs_lower, p_obs_upper)
    list(beta=vals$beta, mu=vals$mu, p_obs=p_obs)
  }
  
  ndraws <- abc_params$ndraws
  max_tries <- abc_params$max_tries
  rmse_threshold <- abc_params$rmse_threshold
  draws <- abc(ndraws, rmse_threshold,
               abc_data$initial_states,
               abc_data$df_real,
               abc_data$D_final_real,
               prior_sample,
               max_tries, print_to_screen = T)
  draws %>%
    mutate(R0=beta / mu * abc_data$max_monkeys)
}