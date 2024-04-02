data {
  int N;
  int day[N];
  int infected[N];
  int died[N];
  
  // priors for the gamma distribution
  real a_1;
  real a_2;
  real b_1;
  real b_2;
  real offset_mean;
}

parameters {
  real<lower=0> a; // gamma dist parameter
  real<lower=0> b; // gamma dist parameter
  real<lower=0,upper=1> p_death_min; // baseline mortality parameter
  real<lower=0.5,upper=5> offset; // mean of the offset prior
  // real<lower=0> a_offset; // gamma dist parameter (but offset)
  // real<lower=0> b_offset; // gamma dist parameter (but offset)
}

model {
  p_death_min ~ normal(0.33, 0.02);
  a ~ uniform(a_1, a_2);
  b ~ uniform(b_1, b_2);
  offset ~ normal(offset_mean, 0.5);
  for (i in 1:N) {
    real temp_day = 0.0; 
    if (day[i] - offset > 0) {
      temp_day = day[i] - offset;
    }
    real mortality_p = p_death_min + (1 - p_death_min) * gamma_cdf(temp_day, a, b);
    died[i] ~ binomial(infected[i], mortality_p);
  }
}

generated quantities {
  real days_simulated = normal_rng(offset, 0.1) + gamma_rng(a, b);
  // real days_simulated_offset = normal_rng(offset, 0.1) + gamma_rng(a_offset, b_offset);
  // {
  //   vector[N] days_sim;
  //   vector[N] days_sim_offset;
  //   for(i in 1:N)
  //     days_sim[i] = gamma_rng(a, b);
  //     days_sim_offset[i] = normal_rng(offset, 0.1) + gamma_rng(a_offset, b_offset);
  // }
}
