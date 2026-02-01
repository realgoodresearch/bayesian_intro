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

# total population
N <- sample(1e3:1e6, 1)
N

# number of repeat surveys conducted
J <- sample(2:5, 1)
J

# probability of an individual being detected in a survey
p <- runif(1, 0.2, 0.8)
p

# individuals captured in each survey round
y <- rbinom(J, N, p)  



#---- run Bayesian model (CmdStandR) ----#

# define location of the stan model file
model_file <- file.path('models', '2b_Nmixture.stan')

# compile the stan model
mod <- cmdstan_model(model_file)



# format data as a list (required by stan)
md <- list(J = J,
           y = y,
           N_true = N,
           p_true = p,
           seed = seed) # important to set and save seed for reproducibility

# save data to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_Nmixture.rds'))



# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()
  
  result[['lambda']] <- runif(1, max(md$y), 1e6)
  result[['p']] <- runif(1, 0.2, 0.8)

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
fit$save_object(file=file.path(outdir, 'fit_Nmixture.rds'))




#---- model diagnostics ----#

# traceplots to visually inspect MCMC chains
bayesplot::mcmc_trace(fit$draws(), pars = c('lambda', 'p'))

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



