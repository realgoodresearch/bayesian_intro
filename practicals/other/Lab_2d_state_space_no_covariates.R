# cleanup
rm(list=ls()); gc(); cat("\014"); try(dev.off(), silent=T)

# load libraries
library(cmdstanr)
library(bayesplot)
library(posterior)
library(ggplot2)

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


## population time series

# number of time steps (years)
Y <- sample(x=10:50, size=1)
Y

# population growth rates
lambda <- rlnorm(n=Y, meanlog=0, sdlog=0.1)
lambda

# population size at time 1 (and expected value)
N <- mu <- sample(x=1e3:1e5, size=1)
N

# create deterministic population time series
for(t in 2:Y){
  mu[t] <- mu[t-1] * lambda[t]
}
mu

# define residual variation in population sizes
sigma <- runif(n = 1, 
               min = 0, 
               max = 0.5)
sigma

# add residual variation to our population sizes
N[2:Y] <- rlnorm(n = Y-1, 
                 meanlog = log(mu[2:Y]), 
                 sdlog = sigma)

# round N so they represent counts 
# (integers will be needed below when N is used in rbinom() function)
N <- round(N)
N

# plot population
plot(y=N, x=1:Y, type='l')


## mark-recapture data

# define the annual probabilities that an individual from the population 
# will be detected in the survey
p <- runif(Y, 0.1, 0.5)
p

# number of individuals "marked" in the first survey round
# (ensure that at least one individual was marked)
m <- rbinom(Y, N, p)
m[m==0] <- 1
m

# number of individuals "captured" in the second survey round
c <- rbinom(Y, N, p)
c

# number of individuals "recaptured in the second survey round
r <- rbinom(Y, c, p)
r

## save data

# format data as a list (required by stan)
md <- list(T = Y, 
           M = m,
           C = c,
           R = r,
           U = c-r+m,
           N_true = N,
           lambda_true = lambda,
           sigma_true = sigma,
           p_true = p,
           seed = seed) # important to set and save seed for reproducibility

# save to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_state_space_no_covariates.rds'))




#---- run MCMC (CmdStandR) ----#

# compile the stan model
mod <- cmdstan_model(file.path('models', '2d_state_space_no_covariates.stan'))


# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()
  
  # Lincoln-Peterson deterministic estimator of population
  pop <- md$M / (md$R / md$C)
  lam <- c(1, pop[2:md$T] / pop[1:(md$T-1)])
  
  result[['mu']] <- rlnorm(md$T, log(pop), 0.1)
  result[['lambda']] <- rlnorm(md$T, log(lam), 0.1)
  result[['sigma']] <- runif(1, 0, 0.1)
  
  result[['N']] <- rlnorm(md$T, log(pop), result[['sigma']])

  return(result)
}


# MCMC config
chains <- 4
warmup <- 2e3
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
fit$save_object(file=file.path(outdir, 'fit_state_space_no_covariates.rds'))




#---- results ----#

# load data
md <- readRDS(file.path(outdir, 'md_state_space_no_covariates.rds'))
fit <- readRDS(file.path(outdir, 'fit_state_space_no_covariates.rds'))


# summary of model results
f <- fit$summary(variables=c('N', 'lambda'))
f
summary(f$rhat)


# traceplots to visually inspect MCMC chains
t <- sample(1:md$T, 6)
bayesplot::mcmc_trace(fit$draws(), pars=paste0('N[', t, ']'))
bayesplot::mcmc_trace(fit$draws(), pars=paste0('lambda[', t, ']'))
bayesplot::mcmc_trace(fit$draws(), pars='sigma')


# true versus predicted population
plot(y = fit$summary(variables='N')$mean, 
     x = md$N_true)
abline(0,1)

# true versus predicted lambda
plot(y=fit$summary(variables='lambda')$mean, x=md$lambda_true)
abline(0,1)

# sigma fit
plot(density(fit$draws('sigma')),
     xlab = 'sigma',
     main = 'Marginal posterior distribution for sigma')
abline(v=md$sigma_true, col='red')
legend('topright',
       legend = c('Posterior', 'True value'),
       lty = c(1, 1),
       col = c('black', 'red'))

# visualise marginal posterior distributions for parameters as density plots
# colour-coded by rhat convergence statistics
pars <- c('sigma')
mcmc_areas(fit$draws(pars), 
           pars = pars, 
           prob = 0.95,
           rhat = bayesplot::rhat(fit)[pars])

pars <- paste0('N[', 1:md$T, ']')
mcmc_areas(fit$draws(), 
           pars = pars, 
           prob = 0.95,
           rhat = bayesplot::rhat(fit)[pars])

pars <- paste0('lambda[', 1:md$T, ']')
mcmc_areas(log(fit$draws(pars)), 
           pars = pars, 
           prob = 0.95,
           rhat = bayesplot::rhat(fit)[pars]) + 
  ggplot2::xlim(-2, 2)


