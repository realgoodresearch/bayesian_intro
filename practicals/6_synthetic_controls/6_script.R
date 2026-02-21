# cleanup
rm(list = ls())
gc()
cat("\014")
try(dev.off(), silent = T)

# load libraries
library(cmdstanr)
library(bayesplot)
library(posterior)
library(tidyr)
library(tibble)
library(dplyr)
library(tidybayes)
library(ggplot2)

# check working directory
# (it should be the root of the "bayesian_intro" github repository)
getwd()

# [optional] update/fix working directory
# setwd(file.path('my_github_repos/realgoodresearch/bayesian_intro'))

# directories
basedir <- file.path(getwd(), "practicals", "6_synthetic_controls")
outdir <- file.path(basedir, "results")
dir.create(outdir, showWarnings = F, recursive = T)


# set seed for random number generators (important for reproducibility)
seed <- round(runif(1, 1, 1e6))
set.seed(seed)


#---- get data ----#

# download Meta global migration data from the Humanitarian Data Exchange
download.file(
  "https://data.humdata.org/dataset/e09595bd-f4c5-4a66-8130-5a05f14d5e64/resource/67b2d6fc-ff4f-4d04-a75e-80e3de81b072/download/international_migration_flow.csv",
  destfile = file.path(outdir, "data_meta_migration.csv")
)

# load data
dat <- read.csv(file.path(outdir, "data_meta_migration.csv"))


#---- USER OPTIONS ----#

# select countries
origin <- 'UA'
destination <- 'GB'

# identify last pre-treatment month
treatment_month <- '2022-01'

# select model file name
model_name <- '6c_model_shrinkage.stan'


#---- setup data ----#

# identify "donor" countries for synthetic controls (ordered as numeric index)
donors <- dat %>%
  filter(
    country_to == destination,
    country_from != origin
  ) %>%
  distinct(country_from) %>%
  arrange() %>%
  pull()

# identify months (ordered as numeric index)
months <- dat %>%
  distinct(migration_month) %>%
  arrange() %>%
  pull()

# get treated (y) data: flows from origin to destination
y_treated <- dat %>%
  filter(country_from == origin, country_to == destination) %>%
  slice(match(months, migration_month)) %>%
  pull(num_migrants) %>%
  replace_na(0)

# get donor (x) data: flows from other countries to destination
x_donors <- dat %>%
  filter(country_to == destination, country_from %in% donors) %>%
  pivot_wider(
    id_cols = migration_month,
    names_from = country_from,
    values_from = num_migrants,
    values_fill = 0
  ) %>%
  slice(match(months, migration_month)) %>%
  column_to_rownames("migration_month") %>%
  select(all_of(donors)) %>%
  as.matrix()

# checks
nrow(x_donors) == length(y_treated)
nrow(x_donors) == length(months)
ncol(x_donors) == length(donors)

# stan model data object
md <- list(
  T = length(months),
  T0 = which(row.names(x_donors) == treatment_month),
  K = length(donors),
  y = y_treated,
  x = x_donors,
  seed = seed
) # important to set and save seed for reproducibility

# save data to disk
saveRDS(
  object = md,
  file = file.path(outdir, "md.rds")
)


#---- run Bayesian model (CmdStandR) ----#

# function to generate initial values for each parameter in the model
init_generator <- function(md = md, model_name = model_name, chain_id = 1) {
  result <- list(
    alpha = rnorm(1, log(mean(md$y[1:md$T0])), 0.1),
    weights = rep(1 / md$K, md$K)
  )
  if (grepl("6b", model_name)) {
    result[['phi']] <- runif(1, 5, 15)
  }
  if (grepl("6c", model_name)) {
    result[['weights']] <- NULL
    result[['weights_raw']] <- rep(1, md$K)
    result[['local_scales']] <- rep(1, md$K)
    result[['global_scale']] <- 1
  }

  return(result)
}


# MCMC config
chains <- 4
warmup <- 1e3
samples <- 2e3
inits <- lapply(1:chains, function(id) {
  init_generator(md = md, model_name = model_name, chain_id = id)
})

# compile the stan model
mod <- cmdstan_model(file.path(basedir, model_name))

# run MCMC
fit <- mod$sample(
  data = md,
  parallel_chains = chains,
  init = inits,
  iter_sampling = samples,
  iter_warmup = warmup,
  adapt_delta = 0.99, # default = 0.80
  max_treedepth = 12, # deafult = 10
  seed = md$seed
)

# save model to disk
fit$save_object(file = file.path(outdir, paste0("fit_", model_name, ".rds")))


#---- model diagnostics ----#

# list of vars to inspect
potential_pars <- c("alpha", "phi")
weight_pars_sample <- paste0("weights[", sample(1:md$K, min(12, md$K)), "]")

# 2. Check which ones actually exist in the fit metadata
existing_vars <- fit$metadata()$variables
target_pars <- c(intersect(potential_pars, existing_vars), weight_pars_sample)

# traceplots to visually inspect MCMC chains
bayesplot::mcmc_trace(
  fit$draws(),
  pars = target_pars[1:12]
)

# Results and MCMC convergence diagnostic (r-hat); rhat < 1.01 means the chains converged (< 1.1 is okay for testing purposes)
# (note: stan will return a warning automatically if there are MCMC problems)
fit$summary(variables = target_pars)

#---- assess coverage ----#

# get prediction intervals
predictive_intervals <- fit$draws("y_rep", format = "df") %>%
  spread_draws(y_rep[t]) %>%
  median_qi(y_rep, .width = 0.95) %>% # Calculates 2.5% and 97.5% quantiles
  filter(t <= md$T0) # Only look at the pre-treatment period

# join with observed data
coverage_df <- predictive_intervals %>%
  mutate(observed = md$y[t]) %>%
  mutate(is_covered = observed >= .lower & observed <= .upper)

# calculate coverage
coverage_rate <- mean(coverage_df$is_covered) * 100

print(paste0("Pre-treatment Coverage: ", round(coverage_rate, 2), "%"))

# visualise coverage
ggplot(coverage_df, aes(x = t)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y_rep), color = "blue", linetype = "dashed") +
  geom_point(aes(y = observed, color = is_covered)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) +
  labs(
    title = paste(
      "Interval Coverage Diagnostic (",
      round(coverage_rate, 1),
      "%)"
    ),
    subtitle = "Red points fall outside the 95% prediction interval",
    x = "Month Index",
    y = "Migration Count"
  ) +
  theme_minimal()


#---- synthetic controls plot ----#

## prepare data for plotting
gen_quant_to_plot <- "y_rep"

# extract the posterior means for y_synthetic or y_rep
draws_df <- fit$draws(
  variables = gen_quant_to_plot,
  format = "df"
)

synth_draws <- draws_df %>%
  pivot_longer(
    cols = starts_with(gen_quant_to_plot),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(t = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", parameter))) %>%
  group_by(t) %>%
  summarise(
    y_hat = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(date = months[t])

# observed (real) data
real_data <- data.frame(
  date = row.names(md$x),
  y_obs = y_treated
)

# data for plot
plot_df <- left_join(synth_draws, real_data, by = "date") %>%
  mutate(date = as.Date(paste0(date, "-01"))) %>%
  mutate(diff = y_obs - y_hat)

treatment_date <- as.Date(paste0(months[md$T0], "-01"))


## plot reality and synthetic control

ggplot(plot_df, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "skyblue", alpha = 0.4) +
  geom_line(aes(y = y_hat, color = "Synthetic Control"), linetype = "dashed") +
  geom_line(aes(y = y_obs, color = "Actual Data")) +
  geom_vline(xintercept = treatment_date, color = "black") +
  scale_y_log10(labels = scales::comma) + # This makes the ribbon visible!
  labs(
    title = "Migration reality and synthetic control timeseries",
    subtitle = paste(origin, "to", destination),
  ) +
  ylab('Migration flows (log scale)') +
  theme_minimal()


## plot difference

ggplot(plot_df, aes(x = date)) +
  geom_ribbon(
    aes(ymin = y_obs - upper, ymax = y_obs - lower),
    fill = "firebrick",
    alpha = 0.4
  ) +
  geom_line(aes(y = diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = treatment_date, color = "black") +
  labs(
    title = "Estimated Treatment Effect (Actual - Synthetic)",
    subtitle = paste(origin, "to", destination),
  ) +
  theme_minimal()


#---- check shrinkage (i.e. dominant donors) ----#

if (grepl("shrink", model_name)) {
  # Set your threshold
  top_x <- 10

  # Summarize weights and map names
  weights_summary <- fit$summary("weights") %>%
    mutate(
      donor_index = 1:n(),
      # Map the index to the name in the 'donors' object
      donor_name = donors[donor_index]
    ) %>%
    arrange(desc(mean)) %>%
    mutate(cumulative_weight = cumsum(mean)) %>%
    slice_head(n = top_x)

  # Plot the top donors with names
  ggplot(weights_summary, aes(x = reorder(donor_name, -mean), y = mean)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_errorbar(
      aes(ymin = q5, ymax = q95),
      width = 0.2,
      color = "firebrick"
    ) +
    geom_text(aes(label = round(mean, 3)), vjust = -1.2, size = 3.5) +
    labs(
      title = paste(
        "Top",
        top_x,
        "donors for synthetic",
        origin,
        "->",
        destination
      ),
      subtitle = paste0(
        "Horseshoe Shrinkage | Total weight: ",
        round(max(weights_summary$cumulative_weight) * 100, 1),
        "%"
      ),
      x = NULL, # Remove axis label as names are self-explanatory
      y = "Posterior Weight (Mean)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      panel.grid.major.x = element_blank()
    )
}
