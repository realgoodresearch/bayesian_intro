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
fit <- readRDS(file.path(
  basedir,
  "..",
  "1_model_of_the_mean",
  "results",
  "fit.rds"
))
md <- readRDS(file.path(
  basedir,
  "..",
  "1_model_of_the_mean",
  "results",
  "md.rds"
))


# load function
source(file.path(basedir, "2_script_fun.R"))

# make plot
myplot(
  fit = fit,
  md = md,
  par = "mu",
  xlim = c(4, 12),
  prior_mean = 10,
  prior_sd = 1
)
