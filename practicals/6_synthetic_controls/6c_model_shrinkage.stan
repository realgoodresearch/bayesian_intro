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
  vector[T] lambda;
  vector<lower=0>[K] weights_shrunk;
  simplex[K] weights;
  
  weights_shrunk = weights_raw .* local_scales * global_scale;
  weights = weights_shrunk / sum(weights_shrunk);
  
  lambda = exp(alpha) .* (x * weights);
}
model {
  // likelihood
  y[1 : T0] ~ neg_binomial_2(lambda[1 : T0], phi);
  
  // priors
  alpha ~ normal(0, 10);
  phi ~ exponential(0.1);
  
  weights_raw ~ std_normal();
  local_scales ~ cauchy(0, 1);
  global_scale ~ cauchy(0, 1);
}
generated quantities {
  vector[T] y_synthetic = lambda;
  vector[T] effect = to_vector(y) - y_synthetic;
  
  // Posterior predictive check: Generate new counts based on the model
  array[T] int y_rep;
  for (t in 1 : T) {
    y_rep[t] = neg_binomial_2_rng(lambda[t], phi);
  }
}
