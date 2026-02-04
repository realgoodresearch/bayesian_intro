data {
  int<lower=0> T; // number of time steps
  vector<lower=0>[T] M; // marked individuals (round 1)
  array[T] int<lower=0> C; // captured individuals (round 2)
  array[T] int<lower=0> R; // recaptured individuals (round 2)
  array[T] real<lower=0> U; // unique individuals (across both rounds)
}

parameters {
  vector<lower=0>[T] N; // population sizes
  vector<lower=0>[T] lambda; // rates of population change
  real<lower=0> sigma; // demographic stochastisticity
}

transformed parameters {
  vector<lower=0>[T] mu;
  
  mu[1] = N[1];
  mu[2:T] = N[1:(T-1)] .* lambda[2:T];
}

model {
  
  // process model
  N ~ lognormal(log(mu), sigma);  

  // observation model
  R ~ binomial(C, M ./ N);

  // priors
  N[1] ~ uniform(U[1], U[1]*100);
  lambda ~ lognormal(0, 0.4);
  sigma ~ uniform(0, 1.5);
}
