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



#---- simulate data (Gaussian) ----#

# define the sample size
n_train <- 1000
n_test <- 500

# number of covariates
K <- 2

# create covariates
x <- matrix(NA, nrow=n_train, ncol=K)
for(k in 1:K){
  x[,k] <- rnorm(n_train, 0, 1)
}

x_test <- matrix(NA, nrow=n_test, ncol=K)
for(k in 1:K){
  x_test[,k] <- rnorm(n_test, 0, 1)
}


# define regression parameters
alpha <- rnorm(1, 0, 1)
beta <- rnorm(K, 0, 3)
sigma <- runif(1, 0, 2)

# generate response variable for model fitting
y <- alpha + as.vector(x %*% beta) + rnorm(n_train, 0, sigma)

# generate out-of-sample data for cross-validation
y_test <-  alpha + as.vector(x_test %*% beta) + rnorm(n_test, 0, sigma)


# visualise covariate relationships with response variable
plot(y~x[,1])
plot(y_test~x_test[,1])

plot(y~x[,2])
plot(y_test~x_test[,2])

# summary statistics of response variable
summary(y)
summary(y_test)

# density plot of response variable
plot(density(y))
rug(y)

plot(density(y_test))
rug(y_test)


# format data as a list (required by stan)
md <- list(n_train = n_train,
           n_test = n_test,
           K = K,
           y = as.vector(y), 
           x = x,
           x_test = x_test,
           alpha_true = alpha,
           beta_true = beta,
           sigma_true = sigma,
           seed = seed) # important to set and save seed for reproducibility

# save data to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_regression.rds'))





#---- run Bayesian model (CmdStandR) ----#

# compile the stan model
mod <- cmdstan_model('regression_xval.stan')


# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()
  
  result[['alpha']] <- runif(1, -10, 10)
  result[['beta']] <- runif(md$K, -10, 10)
  result[['sigma']] <- runif(1, 0, 1e2)
  
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
fit$save_object(file=file.path(outdir, 'fit_regression.rds'))




#---- model diagnostics ----#

# traceplots to visually inspect MCMC chains
bayesplot::mcmc_trace(fit$draws(), 
                      pars = c('alpha', 
                               paste0('beta[',1:md$K,']'), 
                               'sigma'))

# MCMC convergence diagnostics; rhat=1 means the chains converged
# (note: stan will return a warning automatically if there are MCMC problems)
fit$summary(variables=c('alpha', 'beta[1]', 'beta[2]', 'sigma'))

# MCMC sampler diagnostics
# (note: stan will return a warning automatically if there are MCMC problems)
fit$diagnostic_summary()



#---- explore results and compare to simulated data ----#


# parameters from the true unobserved population (i.e. from simulation)
alpha
beta
sigma

# summary of model results
pars <- c('alpha', 'beta[1]', 'beta[2]', 'sigma')
fit$summary(variables=pars)

# visualise marginal posterior distributions for parameters as density plots
bayesplot::mcmc_areas(fit$draws(), pars=pars, prob=0.95)

# ... as histograms
bayesplot::mcmc_hist(fit$draws(), pars=pars)


#---- out-of-sample cross validation ----#

# get yhat values
yhat <- as.data.frame(fit$draws(format='df', variables='y_hat'))

# residuals
resid <- matrix(NA, nrow=nrow(yhat), ncol=md$n_test)
for(i in 1:n_test){
  resid[,i] <- (yhat[,paste0('y_hat[',i,']')] - y[i]) / y[i]
}

# xval metrics
bias <- apply(resid, 1, mean)
summary(bias)
hist(bias)

imprecision <- apply(resid, 1, sd)
summary(imprecision)
hist(imprecision)

inaccuracy <- apply(resid, 1, function(x) {mean(abs(x))})
summary(inaccuracy)
hist(inaccuracy)


# axis limits
ylim <- range(yhat)
xlim <- range(y)

# empty plot
plot(ylim=range(yhat),
     xlim=range(y),
     xlab='observed (y)',
     ylab='predicted (y_hat)',
     main='cross-validation')




