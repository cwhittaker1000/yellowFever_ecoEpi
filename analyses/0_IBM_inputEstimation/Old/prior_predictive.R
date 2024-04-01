library(tidyverse)

a_1 <- 5
b_1 <- 1
a_2 <- 1
b_2 <- 2
n <- 10000
alpha <- rgamma(n, a_1, b_1)
beta <- rgamma(n, a_2, b_2)

delay <- map2_dbl(alpha, beta, ~rgamma(1, .x, .y))
quantile(delay)
