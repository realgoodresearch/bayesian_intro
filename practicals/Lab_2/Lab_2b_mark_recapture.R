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

# population size
N <- sample(1e3:1e5, 1)
N

# probability of an individual being "captured"
p <- runif(1, 0.1, 0.5)
p

# number of individuals "marked" in the first survey round
# ensure that at least one individual was marked
m <- rbinom(1, N, p)
m[m==0] <- 1

# number of individuals "captured" in the second survey round
c <- rbinom(1, N, p)

# number of individuals "recaptured in the second survey round
r <- rbinom(1, c, p)


# format data as a list (required by stan)
md <- list(M = m,
           C = c,
           R = r,
           N_true = N,
           p_true = p,
           seed = seed) # important to set and save seed for reproducibility

# save data to disk
saveRDS(object = md,
        file = file.path(outdir, 'md_mark_recapture.rds'))


#---- run Bayesian model (CmdStandR) ----#

# load data
md <- readRDS(file.path(outdir, 'md_mark_recapture.rds'))

# define location of the stan model file
model_file <- file.path('models', '2b_mark_recapture.stan')

# compile the stan model
mod <- cmdstan_model(model_file)

# function to generate initial values for each parameter in the model
init_generator <- function(md=md, chain_id=1){
  result <- list()
  
  result[['N']] <- runif(1, md$C - md$R + md$M, 1e6)

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
fit$save_object(file=file.path(outdir, 'fit_mark_recapture.rds'))




#---- results ----#

# load data
md <- readRDS(file.path(outdir, 'md_mark_recapture.rds'))
fit <- readRDS(file.path(outdir, 'fit_mark_recapture.rds'))


# traceplots to visually inspect MCMC chains
mcmc_trace(fit$draws(), pars = c('N'))


# summary of model results
fit$summary()


# visualise marginal posterior distributions for parameters as density plots
# colour-coded by rhat convergence statistics
pars <- c('N')
mcmc_areas(fit$draws(pars), 
           pars = pars, 
           prob = 0.95,
           rhat = bayesplot::rhat(fit)[pars])


# assess fit of N
plot(density(fit$draws('N')),
     xlab = 'N',
     main = 'Marginal posterior distribution for N')
abline(v=md$N_true, col='red')
legend('topright',
       legend = c('Posterior', 'True value'),
       lty = c(1, 1),
       col = c('black', 'red'))





