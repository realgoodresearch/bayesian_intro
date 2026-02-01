data {
  int<lower=0> T; // number of time steps
  int<lower=0> K; // number of covariates                   
  matrix[T,K] x; // covariates for population growth rates
  
  vector<lower=0>[T] M; // marked individuals (round 1)
  array[T] int<lower=0> C; // captured individuals (round 2)
  array[T] int<lower=0> R; // recaptured individuals (round 2)
  array[T] real<lower=0> U; // unique individuals (across both rounds)
}

parameters {
  vector<lower=0>[T] N; // population sizes
  real alpha; // intercept for regression on population growth rate (lambda)
  vector[K] beta; // effect sizes for regression on lambda
  real<lower=0> sigma; // demographic stochastisticity
}

transformed parameters {
  vector<lower=0>[T] mu; // expected value of population size
  vector<lower=0>[T] lambda; // population growth rate

  lambda = exp(alpha + x * beta);

  mu[1] = N[1];
  mu[2:T] = N[1:(T-1)] .* lambda[2:T];
}

model {
  
  // process model
  N ~ lognormal(log(mu), sigma);  

  // observation model
  R ~ binomial(C, M ./ N);

  // priors
  N[1] ~ uniform(U[1], 1e7);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 5);
  // sigma ~ uniform(0, 1e5);
}
