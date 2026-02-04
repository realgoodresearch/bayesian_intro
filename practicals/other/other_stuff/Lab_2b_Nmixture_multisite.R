# cleanup
rm(list=ls()); gc(); cat("\014"); try(dev.off(), silent=T)

# load libraries
library(cmdstanr)
library(bayesplot)
library(posterior)

# check working directory
getwd()

# [optional] set different working directory
# setwd('mypath')

# directories
outdir <- file.path(getwd(), 'out')
dir.create(outdir, showWarnings=F, recursive=T)

# set seed for random number generators (important for reproducibility)
seed <- round(runif(1, 1, 1e6))
set.seed(seed)



#---- simulate data ----#

# define the sample size (number of locations surveyed)
N <- 100

# number of covariates to predict population
K <- 2

# we are going to keep rerunning the simulation until we have populations with at least 100 individuals and no more than 10M
rho <- rep(0, N)
while(min(rho) < 100 | max(rho) > 1e6) {

  # create covariates
  x <- matrix(NA, nrow=N, ncol=K)
  for(k in 1:K){
    x[,k] <- rnorm(N, 0, 1)
  }
  
  # define regression parameters
  alpha <- rnorm(1, 0, 10)
  beta <- rnorm(K, 0, 5)
  sigma <- runif(1, 0, 1)
  
  # generate response variable
  # NOTE:  We exponentiate the right side so that the log of the response variable has a linear relationship with the covariates
  rho <- exp(alpha + as.vector(x %*% beta) + rnorm(N, 0, sigma))

}
alpha
beta
sigma

# visualise covariate relationships with the latent parameter (avg population size)
plot(rho~x[,1])
plot(rho~x[,2])

# visualise covariate relationships with log(avg population size)
plot(log(rho)~x[,1])
plot(log(rho)~x[,2])

# summary statistics of latent parameter (avg population size)
summary(rho)
summary(log(rho))

# density plot of latent parameter (avg population size)
plot(density(rho))
rug(rho)

plot(density(log(rho)))
rug(log(rho))



### Simulate repeat survey data

# number of repeat surveys conducted at each site
J <- 3

# setup count data
y <- matrix(NA, nrow=N, ncol=J)

# probability of an individual being detected in a survey
detection <- runif(1, 0.2, 0.8)
detection

# individuals captured in each survey round at each location
for(j in 1:J){
  y[,j] <- rbinom(N, round(rho), detection)  
}



#---- run Bayesian model (CmdStandR) ----#

# define location of the stan model file
model_file <- file.path('models', '2b_mark_recapture.stan')

# compile the stan model
mod <- cmdstan_model(model_file)



# format data as a list (required by stan)
md <- list(N = N,
           K = K,
           x = x,
           m = m,
           c = c,
           r = r,
           rho_min = c - r + m,
           rho_true = rho,
           alpha_true = alpha,
           beta_true = beta,
           sigma_true = sigma,
           detection_true = detection,
           seed = seed) # important to set and save seed for reproducibility

# save data to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_regression.rds'))



# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()
  
  result[['alpha']] <- runif(1, 0, 10)
  result[['beta']] <- runif(md$K, -3, 3)
  result[['sigma']] <- runif(1, 0, 1)

  result[['rho']] <- result[['mu']] <- runif(md$N, md$m/0.8, md$m/0.2)

  return(result)
}


# MCMC config
chains <- 4
warmup <- 1e3
samples <- 2e3
inits <- lapply(1:chains, function(id) init_generator(md=md, chain_id=id))

# run MCMC
fit <- mod$sample(data = md,
                  parallel_chains = chains,
                  init = inits,
                  iter_sampling = samples,
                  iter_warmup = warmup,
                  save_warmup = TRUE,
                  seed = md$seed)

# save model to disk
fit$save_object(file=file.path(outdir, 'fit_mean.rds'))




#---- model diagnostics ----#

# traceplots to visually inspect MCMC chains
bayesplot::mcmc_trace(fit$draws(), 
                      pars = c(paste0('alpha[',1:md$G,']'),
                               paste0('beta[',1:md$K,']'), 
                               'sigma'))

# MCMC convergence diagnostics; rhat=1 means the chains converged
# (note: stan will return a warning automatically if there are MCMC problems)
posterior::rhat(fit$draws('alpha[1]'))
posterior::rhat(fit$draws('alpha[2]'))
posterior::rhat(fit$draws('alpha[3]'))
posterior::rhat(fit$draws('beta[1]'))
posterior::rhat(fit$draws('beta[2]'))
posterior::rhat(fit$draws('sigma'))

# MCMC sampler diagnostics
# (note: stan will return a warning automatically if there are MCMC problems)
fit$diagnostic_summary()



#---- explore results and compare to simulated data ----#


# parameters from the true unobserved population (i.e. from simulation)
alpha
beta
sigma

# summary of model results
pars <- c('alpha', 'beta', 'sigma')
fit$summary(variables=pars)

# visualise marginal posterior distributions for parameters as density plots
pars <- c(paste0('alpha[',1:md$G,']'), 
          paste0('beta[',1:md$K,']'), 
          'sigma')
bayesplot::mcmc_areas(fit$draws(), pars=pars, prob=0.95)

# ... as histograms
bayesplot::mcmc_hist(fit$draws(), pars=pars)



