# cleanup
rm(list=ls()); gc(); cat("\014"); try(dev.off(), silent=T)

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



#---- SECTION 1: Simulate data ----#

# define the sample size
n <- 1e3

# define the mean value
mu <- 5

# define the standard deviation
sigma <- 3

# simulate data by drawing random samples from a normal distribution
# note:  this distribution should be the same as the likelihood that you will use in your model
y <- rnorm(n = n,
           mean = mu,
           sd = sigma)

# summary statistics of data
mean(y)
sd(y)
summary(y)

# histogram of data
hist(y)

# density plot of data
plot(density(y))
rug(y)

# format data as a list (required by stan)
md <- list(n = n,
           y = y, 
           mu_true = mu,
           sigma_true = sigma,
           seed = seed) # important to set and save seed for reproducibility

# save data to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_mean.rds'))





#---- SECTION 2: MCMC for Bayesian model ----#

# load libraries
library(cmdstanr)
library(posterior)
library(bayesplot)

# define location of the stan model file
model_file <- file.path('models', '1a_mean.stan')

# compile the stan model
mod <- cmdstan_model(model_file)

# check the path of this compiled executable file
mod$exe_file()

# print the model to the R console
mod$print()



# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()
  
  result[['mu']] <- runif(1, -1e3, 1e3)
  result[['sigma']] <- runif(1, 0, 1e3)
  
  return(result)
}



# MCMC configuration
chains <- 4
warmup <- 1e3
samples <- 2e3
inits <- lapply(1:chains, function(id) init_generator(md=md, chain_id=id))

# run MCMC to sample from the posterior distribution of our model, given our data
fit <- mod$sample(data = md,
                  parallel_chains = chains,
                  init = inits,
                  iter_sampling = samples,
                  iter_warmup = warmup,
                  save_warmup = TRUE,
                  seed = md$seed)

# save fitted model to disk
fit$save_object(file=file.path(outdir, 'fit_mean.rds'))



#---- SECTION 3: Model diagnostics ----#

# load model and data
md <- readRDS(file.path(outdir, 'md_mean.rds'))
fit <- readRDS(file.path(outdir, 'fit_mean.rds'))

## MCMC diagnostics

# extract MCMC draws as data frame
fit$draws(format='df')

# extract MCMC draws as an array (dimensions = iterations x chains x variables)
fit$draws()

# traceplots to visually inspect MCMC chains
bayesplot::mcmc_trace(fit$draws(), 
                      pars = c('mu', 'sigma'))

# check to make sure the warmup period was long enough
bayesplot::mcmc_trace(fit$draws(inc_warmup=TRUE), 
                      pars = c('mu', 'sigma'), 
                      n_warmup = warmup)


# MCMC convergence diagnostics; rhat=1 means the chains converged
# (note: stan will return a warning automatically if there are MCMC problems)
posterior::rhat(fit$draws('mu'))
posterior::rhat(fit$draws('sigma'))

# MCMC sampler diagnostics
# (note: stan will return a warning automatically if there are MCMC problems)
fit$diagnostic_summary()


## Assess parameter estimates (compare to simulated true parameters)

# mu and sigma from the true unobserved population (i.e. from simulation)
mu
sigma

# mu and sigma from the sample data (i.e. from simulation)
mean(md$y)
sd(md$y)


# summary of model results
fit$summary()

# manually calculate point estimate (i.e. mean) of parameter estimates from MCMC draws
mean(fit$draws('mu'))
mean(fit$draws('sigma'))

# manually calculate 95% credible intervals of parameter estimates from MCMC draws
quantile(fit$draws('mu'), probs=c(0.025, 0.975))
quantile(fit$draws('sigma'), probs=c(0.025, 0.975))


# visualise marginal posterior distributions for parameters as density plots
bayesplot::mcmc_areas(fit$draws('mu'), prob=0.95)
bayesplot::mcmc_areas(fit$draws('sigma'), prob=0.95)

bayesplot::mcmc_areas(fit$draws(), pars=c('mu', 'sigma'), prob=0.95)

# ... as histograms
bayesplot::mcmc_hist(fit$draws('mu'))
bayesplot::mcmc_hist(fit$draws('sigma'))

bayesplot::mcmc_hist(fit$draws(), pars=c('mu', 'sigma'))



