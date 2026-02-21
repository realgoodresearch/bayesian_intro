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
  real phi;
}
transformed parameters {
  vector[T] lambda;
  
  lambda = exp(alpha) .* (x * weights);
}
model {
  // likelihood
  y[1 : T0] ~ neg_binomial_2(lambda[1 : T0], phi);
  
  // priors
  weights ~ dirichlet(rep_vector(1, K));
  
  alpha ~ normal(0, 10);
  phi ~ exponential(0.1);
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
