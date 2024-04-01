library(tidyverse)

nreplicates <- 10000
a <- rgamma(nreplicates, 1, 0.5)
b <- rgamma(nreplicates, 1, 2)
days <- map2_dbl(a, b, ~rgamma(1, .x, .y))
g <- qplot(days) +
  scale_x_log10(limits=c(0.1, NA)) +
  xlab("Days from infection to death") +
  ylab("Count")
ggsave("outputs/monkey_death_prior_predictive.pdf")