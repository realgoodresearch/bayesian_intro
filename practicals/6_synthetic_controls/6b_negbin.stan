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
  vector[T] mu;
  
  mu = exp(alpha) .* (x * weights);
}
model {
  // likelihood
  y[1 : T0] ~ neg_binomial_2(mu[1 : T0], phi);
  
  // priors
  weights ~ dirichlet(rep_vector(1, K));
  
  alpha ~ normal(0, 10);
  phi ~ exponential(0.1);
}
generated quantities {
  array[T] int y_synthetic;
  array[T] int effect;
  
  for (t in 1 : T) {
    y_synthetic[t] = neg_binomial_2_rng(mu[t], phi);
    effect[t] = y[t] - y_synthetic[t];
  }
}
