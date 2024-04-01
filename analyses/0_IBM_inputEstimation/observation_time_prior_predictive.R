library(tidyverse)

nreplicates <- 10000
a <- rgamma(nreplicates, 2, 1)
b <- rgamma(nreplicates, 3, 3)
days <- map2_dbl(a, b, ~rgamma(1, .x, .y))
g <- qplot(days) +
  scale_x_log10(limits=c(0.1, NA)) +
  xlab("Days from death to corpse observation") +
  ylab("Count")
ggsave("outputs/observation_time_prior_predictive.pdf", g,
       width = 10, height = 6)