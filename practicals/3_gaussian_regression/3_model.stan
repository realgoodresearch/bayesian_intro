data {
  int<lower=0> n; // sample size
  int<lower=0> K; // number of predictors
  vector[n] y; // response variable
  matrix[n, K] x; // predictor variables
}

parameters {
  real alpha; // intercept
  vector[K] beta; // slope
  real<lower=0> sigma; // residual variation
}

transformed parameters {
  vector[n] mu; // regression expected value
  
  mu = alpha + x * beta; // x * beta is a matrix multiplication
}

model {
  
  // likelihood
  y ~ normal(mu, sigma);
  
  // priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ cauchy(0, 10);
}
