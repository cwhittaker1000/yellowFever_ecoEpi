# checks that custom functions in Stan code work ok

library(rstan)

expose_stan_functions("src/stan/white_pagano.stan")

# check mu function ------
N <- rpois(10, 3)
R0 <- 2
k <- 4
p <- c(0.1, 0.3, 0.2, 0.4)

mu_1 <- mu(R0, N, p, 1, k)
if(mu_1 != 0)
  stop("mu_1 should equal 0")
mu_2 <- mu(R0, N, p, 2, k)
mu_2_expected <- N[1] * p[1] * R0
if(mu_2 != mu_2_expected)
  stop("mu_2 wrong")

t <- 8
mu_8 <- mu(R0, N, p, t, k)
a_sum <- 0
for(i in 1:k) {
  a_sum <- a_sum + N[t - i] * p[i]
}
a_sum <- a_sum * R0
if(mu_8 != a_sum)
  stop("mu_8 wrong")

k <- 1
p <- c(1)
mu_2 <- mu(R0, N, p, 2, k)
if(mu_2 != N[1] * R0)
  stop("another mu_2 wrong")

# check discretised gamma
k <- 5
alpha <- 8
beta <- 3
gamma_pdfs <- vector(length = k)
for(i in 1:5)
  gamma_pdfs[i] <- discrete_gamma_pdf(i, alpha, beta, k)
if(sum(gamma_pdfs) != 1)
  stop("gamma not normalised")

file_conn <- file("data/processed/test_stan_results.txt")
writeLines(c("All ok with Stan functions"), file_conn)
close(file_conn)