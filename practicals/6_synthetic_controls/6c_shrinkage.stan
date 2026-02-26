// bayesian synthetic controls
// negative binomial likelihood
// horseshoe shrinkage
data {
  int<lower=1> T; // Total number of months
  int<lower=1> T0; // End of pre-treatment period
  int<lower=1> K; // Number of donor countries
  array[T] int y; // Treated country migration counts (Ukraine)
  matrix[T, K] x; // Donor country migration counts/rates
}
parameters {
  real alpha;
  real phi;
  
  vector<lower=0>[K] weights_raw;
  vector<lower=0>[K] local_scales;
  real<lower=0> global_scale;
}
transformed parameters {
  vector[T] mu;
  vector<lower=0>[K] weights_shrunk;
  simplex[K] weights;
  
  weights_shrunk = weights_raw .* local_scales * global_scale;
  weights = weights_shrunk / sum(weights_shrunk);
  
  mu = exp(alpha) .* (x * weights);
}
model {
  // likelihood
  y[1 : T0] ~ neg_binomial_2(mu[1 : T0], phi);
  
  // priors
  alpha ~ normal(0, 10);
  phi ~ exponential(0.1);
  
  weights_raw ~ std_normal();
  local_scales ~ cauchy(0, 1);
  global_scale ~ cauchy(0, 1);
}
generated quantities {
  array[T] int y_synthetic;
  array[T] int effect;
  
  for (t in 1 : T) {
    y_synthetic[t] = neg_binomial_2_rng(mu[t], phi);
    effect[t] = y[t] - y_synthetic[t];
  }
}
