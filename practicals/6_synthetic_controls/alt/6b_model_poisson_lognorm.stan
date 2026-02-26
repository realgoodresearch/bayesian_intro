data {
  int<lower=1> T; // Total number of months
  int<lower=1> T0; // End of pre-treatment period
  int<lower=1> K; // Number of donor countries
  array[T] int y; // Treated country migration counts (Ukraine)
  matrix[T, K] x; // Donor country migration counts/rates
}
parameters {
  real alpha;
  simplex[K] weights;
  real<lower=0> sigma_eps;
  vector[T] epsilon_raw;
}
transformed parameters {
  vector[T] lambda;
  vector[T] epsilon = epsilon_raw * (sigma_eps + 1e-4);
  
  lambda = exp(alpha + epsilon) .* (x * weights);
}
model {
  // likelihood
  y[1 : T0] ~ poisson(lambda[1 : T0]);
  
  // priors
  weights ~ dirichlet(rep_vector(1, K));
  
  alpha ~ normal(0, 10);
  sigma_eps ~ exponential(1);
  epsilon_raw ~ std_normal();
}
generated quantities {
  vector[T] y_synthetic = lambda;
  vector[T] effect = to_vector(y) - y_synthetic;
  
  // Posterior predictive check: Generate new counts based on the model
  array[T] int y_rep;
  for (t in 1 : T) {
    y_rep[t] = poisson_rng(lambda[t]);
  }
}
