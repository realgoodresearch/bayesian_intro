data {
  int<lower=0> N; // sample size
  int<lower=0> K; // number of predictors
  int<lower=0> J; // number of survey rounds at each site
  matrix<lower=0>[N, J] m; // marked individuals (round 1)
  array[N] int<lower=0> c; // captured individuals (round 2)
  array[N] int<lower=0> r; // recaptured individuals (round 2)
  matrix[N, K] x; // predictor variables
}

parameters {
  vector<lower=rho_min>[N] rho; // population
  real alpha; // intercept
  vector[K] beta; // slopes
  real<lower=0> sigma; // residual variation
}

transformed parameters {
  vector[N] mu; // regression expected value
  vector<lower=0, upper=1>[N] theta; // detection probability
  
  mu = alpha + x * beta;
  
  theta = m ./ rho;
}

model {
  
  // process model
  rho ~ lognormal(mu, sigma);
  
  // observation model
  r ~ binomial(c, theta); // note: ./ is element-wise division of vectors

  // priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ cauchy(0, 10);
}
