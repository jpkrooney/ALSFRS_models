library(tidyverse)
library(tableone)
library(brms)
library(bayesplot)
library(rstan)

cache_dir <- here::here("local_temp_data")

# Load models
fit_dims1 <- readRDS( paste0(cache_dir, "/ALSFRSdimensions1.rds") )
fit_multv <- readRDS( paste0(cache_dir, "/ALSFRSmultivariate.rds") )
fit_cratio <- readRDS( paste0(cache_dir, "/ALSFRScratio.rds") )
fit_custom <- readRDS( paste0(cache_dir, "/ALSFRScustom.rds") )

#
get_elapsed_time(fit_dims1)
get_elapsed_time(fit_multv)
get_elapsed_time(fit_cratio)
get_elapsed_time(fit_custom)


summary(fit_dims1)
summary(fit_multv)
summary(fit_cratio)
summary(fit_custom)


# Get fixef
fixef(fit_dims1)
fixef(fit_multv)
fixef(fit_cratio)
fixef(fit_custom)

# Get VarCov matrices
vcov_dims1 <- VarCorr(fit_dims1)
vcov_multv <- VarCorr(fit_multv)
vcov_cratio <- VarCorr(fit_cratio)
vcov_custom <- VarCorr(fit_custom)



vcov_dims1[[1]][[2]][,,12]
vcov_multv[[1]][[2]]
vcov_cratio[[1]][[2]][,,12]
vcov_custom[[1]][[2]][,,12]
