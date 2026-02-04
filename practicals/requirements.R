
## install R dependencies
# (see https://mc-stan.org/cmdstanr/articles/cmdstanr.html)


# function to check for packages and install them if missing
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("\n\n---------- Installing package:", pkg, " ----------\n\n"))
      if (pkg == "cmdstanr") {
        install.packages(
          "cmdstanr",
          repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
        )
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    }
  }
}

# list of required R packages
pkgs <- c(
  "lme4",
  "cmdstanr",
  "posterior",
  "bayesplot"
)

# install R packages
install_if_missing(pkgs)

# load libraries
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

# check CmdStan toolchain is available
check_cmdstan_toolchain()

# if needed, download and run the Rtools45 installer from here: 
# https://cran.rstudio.com/bin/windows/Rtools/rtools45/rtools.html

# install CmdStan
install_cmdstan(cores = 2)

# confirm
cmdstan_version()
