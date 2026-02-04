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


## covariates to predict annual population growth rates

# number of time steps (years)
Y <- sample(x=10:50, size=1)
Y

# number of covariates
K <- 2

# create covariates
x <- matrix(NA, nrow=Y, ncol=K)
for(k in 1:K){
  x[,k] <- rnorm(n=Y, mean=0, sd=1)
}
head(x)

## population growth rates

# this "while loop" will re-sample simulation parameters until we get 
# population growth rates (i.e. lambdas on the natural scale) that are within 
# our specified criteria
lambda <- rnorm(n=Y, mean=10, sd=10) # initiate lambda with really bad values
while(abs(mean(lambda)) > 0.1 | sd(lambda) > 0.2){
  
  # draw random regression coefficients
  alpha <- rnorm(n=1, mean=0, sd=1)
  beta <- rnorm(n=K, mean=0, sd=1)
  
  # calculate annual rates of population change from covariates
  lambda <- alpha + as.vector(x %*% beta)

}
alpha
beta

hist(lambda)
plot(y=lambda, x=1:Y)

# exponentiate lambda
lambda <- exp(lambda)

hist(lambda)
plot(y=lambda, x=1:Y)


## population time series

# population size at time 1 (and expected value)
N <- mu <- sample(x=1e3:1e5, size=1)
N

# create deterministic population time series
for(t in 2:Y){
  mu[t] <- mu[t-1] * lambda[t]
}

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
           K = K,
           M = m,
           C = c,
           R = r,
           U = c-r+m,
           x = x,
           N_true = N,
           lambda_true = lambda,
           alpha_true = alpha,
           beta_true = beta,
           sigma_true = sigma,
           p_true = p,
           seed = seed) # important to set and save seed for reproducibility

# save to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_state_space.rds'))




#---- run Bayesian model (CmdStandR) ----#

# load data
md <- readRDS(file.path(outdir, 'md_state_space.rds'))

# compile the stan model
mod <- cmdstan_model(file.path('models', '2c_state_space.stan'))


# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()
  
  # Lincoln-Peterson deterministic estimator of population
  pop <- md$M / (md$R / md$C)

  result[['alpha']] <- rnorm(1, alpha, 0.1)
  result[['beta']] <- rnorm(md$K, beta, 0.1)
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
fit$save_object(file=file.path(outdir, 'fit_state_space.rds'))



#---- results ----#

# load data
md <- readRDS(file.path(outdir, 'md_state_space.rds'))
fit <- readRDS(file.path(outdir, 'fit_state_space.rds'))


# summary of model results
f <- fit$summary(variables=c('N', 'lambda'))
f

# check convergence (rhat=1 is converged)
summary(f$rhat)


# traceplots to visually inspect MCMC chains
mcmc_trace(fit$draws(), pars=c('sigma',
                               'alpha',
                               paste0('beta[', 1:md$K, ']')))

t <- sample(1:md$T, 6)
mcmc_trace(fit$draws(), pars=paste0('N[', t, ']'))
mcmc_trace(fit$draws(), pars=paste0('lambda[', t, ']'))


# visualise marginal posterior distributions for parameters as density plots
# colour-coded by rhat convergence statistics
pars <- c('alpha', paste0('beta[', 1:md$K, ']'), 'sigma')
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
mcmc_areas(fit$draws(pars), 
           transformations = log,
           pars = pars, 
           prob = 0.95,
           rhat = bayesplot::rhat(fit)[pars]) + 
  ggplot2::xlim(-2, 2)



# true versus predicted population
plot(y = fit$summary(variables='N')$mean, 
     x = md$N_true,
     ylab = 'Predicted (mean of posterior)',
     xlab = 'True value',
     main = 'Population size (N)')
abline(0,1) # 1:1 line

# true versus predicted lambda
plot(y = fit$summary(variables='lambda')$mean, 
     x = md$lambda_true,
     ylab = 'Predicted (mean of posterior)',
     xlab = 'True value',
     main = 'Population growth rate (lambda)')
abline(0,1) # 1:1 line

# sigma fit
plot(density(fit$draws('sigma')),
     xlab = 'sigma',
     main = 'Marginal posterior distribution for sigma')
abline(v=md$sigma_true, col='red')
legend('topright',
       legend = c('Posterior', 'True value'),
       lty = c(1, 1),
       col = c('black', 'red'))




