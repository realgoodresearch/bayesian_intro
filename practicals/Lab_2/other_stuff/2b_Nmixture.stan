data {
  int<lower=0> J; // number of survey rounds
  array[J] int<lower=0> y; // counts of individuals in each survey round
}

parameters {
  int<lower=0> N; // total population size
  real<lower=max(y)> lambda; // expected value of total population size
  real<lower=0, upper=1> p; // detection probability
}

model {
  
  // likelihood
  y ~ binomial(N, p);
  
  // priors
  p ~ beta(1, 1);
  N ~ poisson(lambda);
  lambda ~ uniform(max(y), 1e6);
}
