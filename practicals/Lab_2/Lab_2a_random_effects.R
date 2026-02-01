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

# define the sample size
n <- 1e3

# number of covariates
K <- 2

# number of groups for random intercept
G <- 3

# assign samples to groups
group <- sample(1:G, size=n, replace=TRUE, prob=rep(1/G,G))

# create covariates
x <- matrix(NA, nrow=n, ncol=K)
for(k in 1:K){
  x[,k] <- rnorm(n, 0, 1)
}

# define regression parameters
alpha <- rnorm(G, 0, 3)
beta <- rnorm(K, 0, 3)
sigma <- runif(1, 0, 2)

# generate response variable
y <- alpha[group] + as.vector(x %*% beta) + rnorm(n, 0, sigma)

# visualise covariate relationships with response variable
# assign colours based on group membership (i.e. group for the random intercept)
plot(y~x[,1], col=as.factor(group))
plot(y~x[,2], col=as.factor(group))

# summary statistics of response variable
summary(y)

# density plot of response variable
plot(density(y))
rug(y)


# format data as a list (required by stan)
md <- list(n = n,
           K = K,
           G = G,
           y = y, 
           x = x,
           group = group,
           alpha_true = alpha,
           beta_true = beta,
           sigma_true = sigma,
           seed = seed) # important to set and save seed for reproducibility

# save data to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_random_intercept.rds'))



#---- run Bayesian model (CmdStandR) ----#

# load data
md <- readRDS(file.path(outdir, 'md_random_intercept.rds'))

# compile the stan model
mod <- cmdstan_model(file.path('models', '2a_random_intercept.stan'))

# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()
  
  result[['beta']] <- runif(md$K, -5, 5)
  result[['sigma']] <- runif(1, 0, 10)
  
  result[['alpha']] <- runif(md$G, -10, 10)
  result[['mu_alpha']] <- runif(1, -5, 5)
  result[['sigma_alpha']] <- runif(1, 0, 5)
  
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
fit$save_object(file=file.path(outdir, 'fit_random_intercept.rds'))



#---- results ----#

# load data
md <- readRDS(file.path(outdir, 'md_random_intercept.rds'))
fit <- readRDS(file.path(outdir, 'fit_random_intercept.rds'))


# traceplots to visually inspect MCMC chains
bayesplot::mcmc_trace(fit$draws(), 
                      pars = c(paste0('alpha[',1:md$G,']'),
                               paste0('beta[',1:md$K,']'), 
                               'sigma'))


# summary of model results
f <- fit$summary(variables=c('alpha', 'beta', 'sigma'))
f$true <- c(md$alpha, md$beta, md$sigma)
f$predicted <- f$mean
f$error <- f$predicted - f$true
f

# visualise marginal posterior distributions for parameters as density plots
# colour-coded by rhat convergence statistics
pars <- c(paste0('alpha[',1:md$G, ']'), 
          paste0('beta[',1:md$K, ']'), 
          'sigma')
mcmc_areas(fit$draws(pars), 
           pars = pars, 
           prob = 0.95,
           rhat = bayesplot::rhat(fit)[pars])



