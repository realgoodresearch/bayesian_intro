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
}
transformed parameters {
  vector[T] lambda;
  lambda = exp(alpha) * (x * weights);
}
model {
  // likelihood
  y[1 : T0] ~ poisson(lambda[1 : T0]);
  
  // priors
  alpha ~ normal(0, 10);
  weights ~ dirichlet(rep_vector(1, K));
}
generated quantities {
  vector[T] y_synthetic = lambda;
  vector[T] effect = to_vector(y) - y_synthetic;
  
  // posterior predictive check
  array[T] int y_rep;
  for (t in 1 : T) {
    y_rep[t] = poisson_rng(lambda[t]);
  }
}
