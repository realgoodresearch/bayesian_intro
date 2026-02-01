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



#---- simulate data (Gaussian) ----#

# define the sample size
n <- 1e3

# define the mean value
mu <- 5

# define the standard deviation
sigma <- 3

# simulate data by drawing random samples from the distribution
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



#---- run Bayesian model (rstan) ----#

# load libraries
library(rstan)

# rstan recommended config
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)


# define stan model file
model_file <- file.path('models', '1_mean.stan')


# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()

  result[['mu']] <- runif(1, -1e3, 1e3)
  result[['sigma']] <- runif(1, 0, 1e3)

  return(result)
}

# identify parameters to save MCMC results
pars <- c('mu', 'sigma')


# format data as a list (required by stan)
md <- list(n = length(y),
           y = y,
           mu_true = mu,
           sigma_true = sigma,
           seed = seed) # important to set and save seed for reproducibility

# save data to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_mean.rds'))


# mcmc config
chains <- 4
warmup <- 1e3
samples <- 2e3
inits <- lapply(1:chains, function(id) init_generator(md=md, chain_id=id))

# run model
fit <- rstan::stan(file = model_file,
                   model_name = 'Our first Bayesian model (model of the mean)',
                   data = md,
                   pars = pars,
                   chains = chains,
                   init = inits,
                   iter = warmup + samples,
                   warmup = warmup,
                   seed = md$seed)

# save model to disk
saveRDS(fit, file.path(outdir, 'fit_mean.rds'))



#---- model diagnostics ----#

# view MCMC samples as data frame
fit_df <- as.data.frame(fit)
View(fit_df)

# traceplots to inspect MCMC chains
traceplot(object = fit,
          pars = c('mu', 'sigma'))

# check the warmup period
traceplot(object = fit,
          pars = c('mu', 'sigma'),
          inc_warmup = TRUE)

# summary of results
print(fit)

# manually calculate 95% credible intervals of mu estimate from MCMC samples
quantile(fit_df$mu, probs=c(0.025, 0.975))

# manually calculate point estimate (i.e. mean) of mu estimate from MCMC samples
mean(fit_df$mu)



#---- compare results to simulated population ----#

# compare to mu and sigma from the sample
mean(md$y)
sd(md$y)

# compare to mu and sigma from the true unobserved population
mu
sigma

# visualise marginal posterior distributions for parameters as boxplots
stan_plot(fit)

# marginal posterior density plots (probabilistic parameter estimates)
plot(density(fit_df$mu))
rug(fit_df$mu)

plot(density(fit_df$sigma))
rug(fit_df$sigma)


