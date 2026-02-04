# cleanup
rm(list = ls())
gc()
cat("\014")
try(dev.off(), silent = T)
options(scipen = 999)

# libraries
library(cmdstanr)

# check working directory
# (it should be the root of the "bayesian_intro" github repository)
getwd()

# [optional] update/fix working directory
# setwd(file.path('my_github_repos/realgoodresearch/bayesian_intro'))

# directories
basedir <- file.path(getwd(), 'practicals', '2_likelihoods_priors_posteriors')
outdir <- file.path(basedir, 'results')
dir.create(outdir, showWarnings = F, recursive = T)

# load model and data from activity 1
fit <- readRDS(file.path(basedir, "..", "1_model_of_the_mean", "results", "fit.rds"))
md <- readRDS(file.path(basedir, "..", "1_model_of_the_mean", "results", "md.rds"))


#---- SECTION 1: A couple of functions to help ----#

# log-likelihood function
# (the density function--e.g. dnorm()-- should match the distribution of the likelihood from the model being assessed)
Lfun <- function(mu, sigma, dat = md$y) {
  sum(dnorm(x = dat, mean = mu, sd = sigma, log = TRUE))
}

# rescale a vector of log-likelihoods to have the min and max provided
rescale_log_likelihoods <- function(log_likelihoods, new_min, new_max) {
  # log-sum-exp trick (thanks, ChatGPT 4o)
  max_log_likelihoods <- max(log_likelihoods)
  likelihoods <- exp(log_likelihoods - max_log_likelihoods)

  # rescale (thanks again, ChatGPT 4o)
  old_min <- min(likelihoods, na.rm = TRUE)
  old_max <- max(likelihoods, na.rm = TRUE)
  scaled_likelihoods <- ((likelihoods - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min

  return(scaled_likelihoods)
}


#---- SECTION 2: Calculate values to plot ----#

# posterior estimate for mu
#** Note:  Update "from" and "to" define the range of values you want to plot (below) *
density_mu <- density(
  fit$draws("mu"),
  from = 4,
  to = 6
)

# take a look at the result from the density function
# y = probability densities from our posterior estimate for mu
# x = values of mu
density_mu
plot(density_mu)

# extract the probability density values (i.e. the values plotted on the y-axis)
posterior_mu <- density_mu$y


# probability density of the prior defined in the model
#** Note: This should be updated to match the prior for mu that was defined in the model being assessed *
prior_mu <- dnorm(
  x = density_mu$x,
  mean = 0,
  sd = 100
)


# likelihood estimate: i.e. the likelihood of the data given mu, i.e. L(x|mu)
log_likelihood_mu <- sapply(
  density_mu$x,
  function(m) {
    Lfun(
      mu = m,
      sigma = md$sigma_true
    )
  }
)

# rescale the log-likelihood values to match scale of posterior probabilities
# Note: This is necessary just to make a nice plot with posteriors, priors, and likelihoods together

likelihood_mu <- rescale_log_likelihoods(
  log_likelihood_mu,
  new_min = min(posterior_mu),
  new_max = max(posterior_mu)
)


#----  Section 3:  Plot posterior, prior, and likelihood ----#

# posterior
plot(
  x = density_mu$x,
  y = posterior_mu,
  type = "l",
  lwd = 0,
  col = "grey",
  xlab = "Parameter value",
  ylab = "Probability density or scaled Likelihood"
)

polygon(
  x = c(density_mu$x, rev(density_mu$x)),
  y = c(posterior_mu, rep(0, length(posterior_mu))),
  col = "grey",
  border = NA
)

# prior
lines(x = density_mu$x, y = prior_mu, lty = 3, lwd = 2)

# likelihood
lines(x = density_mu$x, y = likelihood_mu, lty = 2, lwd = 2)

# true value
abline(v = md$mu_true, col = "red")

# legend
legend(
  "topright",
  legend = c("Posterior", "Prior", "Likelihood", "True value"),
  fill = c("grey", NA, NA, NA),
  border = c(NA, NA, NA, NA),
  lty = c(NA, 3, 2, 1),
  lwd = c(NA, 2, 2, 2),
  col = c(NA, "black", "black", "red")
)
