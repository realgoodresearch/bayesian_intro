data {
  int<lower=0> M; // marked individuals (round 1)
  int<lower=0> C; // captured individuals (round 2)
  int<lower=0> R; // recaptured individuals (round 2)
}

parameters {
  real<lower=(C-R+M)> N; // population size
}

model {
  
  // mark-recapture model
  R ~ binomial(C, M / N);

  // priors
  N ~ uniform(C-R+M, 1e6);
}
