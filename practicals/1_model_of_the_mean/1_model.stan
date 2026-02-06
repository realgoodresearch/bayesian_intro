data {
  int<lower=0> n; // sample size
  vector[n] y; // data
}
parameters {
  real mu; // mean
  real<lower=0> sigma; // standard deviation
}
model {
  // model of the mean
  y ~ normal(mu, sigma);
  
  // priors
  mu ~ normal(10, 0.1); // default N(0,100)
  sigma ~ uniform(0, 1e4);
}
