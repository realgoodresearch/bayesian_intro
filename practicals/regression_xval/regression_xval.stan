data {
  int<lower=0> n_train; // sample size
  int<lower=0> n_test; // sample size for xval
  int<lower=0> K; // number of predictors
  vector[n_train] y; // response variable
  matrix[n_train, K] x; // predictor variables
  matrix[n_test, K] x_test; // predictor variables
}

parameters {
  real alpha; // intercept
  vector[K] beta; // slope
  real<lower=0> sigma; // residual variation
}

transformed parameters {
  vector[n_train] mu; // regression expected value
  
  mu = alpha + x * beta;
}

model {
  
  // likelihood
  y ~ normal(mu, sigma);
  
  // priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ cauchy(0, 10);
}

generated quantities {
  array[n_test] real y_hat; // response variable for test data

  y_hat = normal_rng(alpha + x_test * beta, sigma);
}
