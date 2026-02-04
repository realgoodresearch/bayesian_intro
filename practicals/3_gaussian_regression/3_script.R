# cleanup
rm(list = ls())
gc()
cat("\014")
try(dev.off(), silent = T)

# load libraries
library(cmdstanr)
library(bayesplot)
library(posterior)

# check working directory
# (it should be the root of the "bayesian_intro" github repository)
getwd()

# [optional] update/fix working directory
# setwd(file.path('my_github_repos/realgoodresearch/bayesian_intro'))

# directories
basedir <- file.path(getwd(), "practicals", "3_gaussian_regression")
outdir <- file.path(basedir, "results")
dir.create(outdir, showWarnings = F, recursive = T)

# set seed for random number generators (important for reproducibility)
seed <- round(runif(1, 1, 1e6))
set.seed(seed)



#---- simulate data (Gaussian) ----#

# define the sample size
n <- 1e3

# number of covariates (use K > 1)
K <- 2
if (K < 2) stop("Please define K as 2 or more.")

# create covariates
x <- matrix(NA, nrow = n, ncol = K)
for (k in 1:K) {
  x[, k] <- rnorm(n, 0, 1)
}

# define regression parameters (randomly)
alpha <- rnorm(1, 0, 1)
beta <- rnorm(K, 0, 3)
sigma <- runif(1, 0, 2)

# generate response variable
y <- alpha + as.vector(x %*% beta) + rnorm(n, 0, sigma)

# note: %*% is a matrix multiplication of a matrix of x values times a vector of beta values
# convince yourself that (x %*% beta) is equal to (beta1 * x1 + beta2 * x2)
head(x[, 1:2] %*% beta[1:2])
head(beta[1] * x[, 1] + beta[2] * x[, 2])

# visualise covariate relationships with response variable
plot(y ~ x[, 1])
plot(y ~ x[, 2])

# summary statistics of response variable
summary(y)

# density plot of response variable
plot(density(y))
rug(y)



#---- run Bayesian model (CmdStandR) ----#

# format data as a list (required by stan)
md <- list(
  n = n,
  K = K,
  y = as.vector(y),
  x = x,
  alpha_true = alpha,
  beta_true = beta,
  sigma_true = sigma,
  seed = seed # important to set and save seed for reproducibility
)

# save data to disk
saveRDS(
  object = md,
  file = file.path(outdir, "md.rds")
)


# define location of the stan model file
# NOTE: you need to revise the model provided to make it into regression!!
model_file <- file.path(basedir, "3_model.stan")

# compile the stan model
mod <- cmdstan_model(model_file)


# function to generate initial values for each parameter in the model
init_generator <- function(md = md, chain_id = 1) {
  result <- list()

  result[["alpha"]] <- runif(1, -10, 10)
  result[["beta"]] <- runif(md$K, -10, 10)
  result[["sigma"]] <- runif(1, 0, 1e2)

  return(result)
}


# MCMC config
chains <- 4
warmup <- 1e3
samples <- 2e3
inits <- lapply(1:chains, function(id) init_generator(md = md, chain_id = id))

# run MCMC
fit <- mod$sample(
  data = md,
  parallel_chains = chains,
  init = inits,
  iter_sampling = samples,
  iter_warmup = warmup,
  save_warmup = TRUE,
  seed = md$seed
)

# save model to disk
fit$save_object(file = file.path(outdir, "fit.rds"))


#---- model diagnostics ----#

# traceplots to visually inspect MCMC chains
bayesplot::mcmc_trace(
  fit$draws(),
  pars = c(
    "alpha",
    paste0("beta[", 1:md$K, "]"),
    "sigma"
  )
)

# MCMC convergence diagnostics; rhat=1 means the chains converged
# (note: stan will return a warning automatically if there are MCMC problems)
posterior::rhat(fit$draws("alpha"))
posterior::rhat(fit$draws("beta[1]"))
posterior::rhat(fit$draws("beta[2]"))
posterior::rhat(fit$draws("sigma"))

# MCMC sampler diagnostics
# (note: stan will return a warning automatically if there are MCMC problems)
fit$diagnostic_summary()



#---- explore results and compare to simulated data ----#


# parameters from the true unobserved population (i.e. from simulation)
alpha
beta
sigma

# summary of model results
pars <- c("alpha", "beta", "sigma")
fit$summary(variables = pars)

# visualise marginal posterior distributions for parameters as density plots
pars <- c("alpha", paste0("beta[", 1:md$K, "]"), "sigma")
bayesplot::mcmc_areas(fit$draws(), pars = pars, prob = 0.95)

# ... as histograms
bayesplot::mcmc_hist(fit$draws(), pars = pars)
