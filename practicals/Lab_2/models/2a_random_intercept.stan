data {
  int<lower=0> n; // sample size
  int<lower=0> K; // number of predictors
  int<lower=0> G; // number of groups for random intercept
  array[n] int<lower=0, upper=G> group; // group
  vector[n] y; // response variable
  matrix[n, K] x; // predictor variables
}

parameters {
  vector[G] alpha; // random intercept
  vector[K] beta; // slope
  real<lower=0> sigma; // residual variation
  real mu_alpha; // mean alpha among groups
  real<lower=0> sigma_alpha; // standard deviation of alpha among groups
}

transformed parameters {
  vector[n] mu; // regression expected value
  
  mu = alpha[group] + x * beta;
}

model {
  
  // likelihood
  y ~ normal(mu, sigma);
  
  // random intercept
  alpha ~ normal(mu_alpha, sigma_alpha);
  
  // priors
  sigma ~ cauchy(0, 5);
  beta ~ normal(0, 10);
  mu_alpha ~ normal(0, 10);
  sigma_alpha ~ cauchy(0, 5);
}
