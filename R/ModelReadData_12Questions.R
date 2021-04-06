
# Load required packages
library(tidyverse)
library(tableone)
library(brms)
library(bayesplot)
source(here::here("R", "mv_probit.R"))

options(mc.cores = parallel::detectCores(), brms.backend = "cmdstanr")
cache_dir <- here::here("local_temp_data")

threads = 2

# Load data
df_ind <- read.csv("Data/casedata.csv", stringsAsFactors = FALSE)
df_frs <- read.csv("Data/longdata.csv", stringsAsFactors = FALSE)

# rename question variables to model specifications shorter

df_frs <- df_frs %>% rename(
    "Q01" = "ALSFRS_Q1",
    "Q02" = "ALSFRS_Q2",
    "Q03" = "ALSFRS_Q3",
    "Q04" = "ALSFRS_Q4",
    "Q05" = "ALSFRS_Q5",
    "Q06" = "ALSFRS_Q6",
    "Q07" = "ALSFRS_Q7",
    "Q08" = "ALSFRS_Q8",
    "Q09" = "ALSFRS_Q9",
    "Q10" = "ALSFRS_Q10",
    "Q11" = "ALSFRS_Q11",
    "Q12" = "ALSFRS_Q12"
)

# add 1 to answers as 0 not allowed in cumulative models
df_frs$Q01 <- df_frs$Q01 + 1
df_frs$Q02 <- df_frs$Q02 + 1
df_frs$Q03 <- df_frs$Q03 + 1
df_frs$Q04 <- df_frs$Q04 + 1
df_frs$Q05 <- df_frs$Q05 + 1
df_frs$Q06 <- df_frs$Q06 + 1
df_frs$Q07 <- df_frs$Q07 + 1
df_frs$Q08 <- df_frs$Q08 + 1
df_frs$Q09 <- df_frs$Q09 + 1
df_frs$Q10 <- df_frs$Q10 + 1
df_frs$Q11 <- df_frs$Q11 + 1
df_frs$Q12 <- df_frs$Q12 + 1

# Add baseline variables of potential interest to df_frs
df_frs$site <- df_ind[ match(df_frs$ID, df_ind$ID ), ]$simp_site
df_frs$sex <- df_ind[ match(df_frs$ID, df_ind$ID ), ]$sex
df_frs$age_dx <- df_ind[ match(df_frs$ID, df_ind$ID ), ]$age_dx
df_frs$dx_delay <- df_ind[ match(df_frs$ID, df_ind$ID ), ]$dx_delay

# Define a time variable relative to first ALSFRS variable
firsts <- df_frs[ df_frs$First %in% c("First", "First & Final") ,]
df_frs$t_first <- firsts[ match(df_frs$ID, firsts$ID), ]$alsfrs_dly_mnths
df_frs$time_first <- df_frs$alsfrs_dly_mnths - df_frs$t_first


# Check for and exclude rows with missing values as they can cause problems
df_frs <- df_frs[ rowMeans( is.na(df_frs) ) ==0, ] # Note this fairly brutal screen as we don't use all the columns



#############################################
# Minimise the data size  for test purposes ####
# Take 15% of individuals
samp <- sample(df_ind$ID, size = nrow(df_ind) * 0.15, replace = FALSE)
df_frs <- df_frs[df_frs$ID %in% samp, ]
# Will only use Q1 to Q6
#df_frs <- df_frs[ , !names(df_frs) %in% c("Q07", "Q08", "Q09", "Q10", "Q11", "Q12")]

#############################################


# make long version of data
frs_long <- df_frs %>%
    pivot_longer(starts_with("Q"), names_to = "question", values_to = "answer")

# save both versions of data
saveRDS(df_frs, "Data/long_format_wide12Qs.RDS")
saveRDS(frs_long, "Data/long_format_long12Qs.RDS")



## Model 6 - Multivariate probit - bgoodri's parametrization

#Code based on https://github.com/stan-dev/example-models/commit/d6f0282d64382b627dfddca6b7f9a551bda3f537 and advice on Stan Discourse by @CerulloE (Thanks everybody!)

sv <- make_stanvars_mv_probit_bgoodri(
    c("Q01", "Q02", "Q03", "Q04", "Q05", "Q06", "Q07", "Q08", "Q09", "Q10", "Q11", "Q12"))


fit_rescor_bgoodri <- brm(bf(
    mvbind(Q01, Q02, Q03, Q04, Q05, Q06,
           Q07, Q08, Q09, Q10, Q11, Q12) ~ 1 + alsfrs_dly_mnths + (alsfrs_dly_mnths | p | ID),
    family = empty_cumulative()) + set_rescor(FALSE),
    file = paste0(cache_dir, "/mvprobit_bgoodri_12Qs.rds"),
    file_refit = "on_change",
    data = df_frs, stanvars = sv, adapt_delta = 0.95,
    init = 0.005)

fit_rescor_bgoodri

mcmc_trace(fit_rescor_bgoodri, pars = c("b_Q01_alsfrs_dly_mnths", "b_Q02_Intercept[2]"))

