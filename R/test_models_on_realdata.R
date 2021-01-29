#' @title Irish ALSFRS Total Correlation Explanation
#' @author Dr. James Rooney, TCD
#' @details
#' This R Script performs a Total Correlation Explanation analysis of longtiudinal ALSFRS data from the Irish
#' ALS register.

# Load required packages
library(tidyverse)
library(tableone)
library(brms)
library(bayesplot)

cache_dir <- here::here("local_temp_data")

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

# make long version of data
frs_long <- df_frs %>%
    pivot_longer(starts_with("Q"), names_to = "question", values_to = "answer")

# save both versions of data
saveRDS(df_frs, "Data/long_format_wide.RDS")
saveRDS(frs_long, "Data/long_format_long.RDS")

# Fit first model
fit_dimensions1 <- brm(answer ~ question + alsfrs_dly_mnths + (0 + question  | ID),
                       family = cumulative("logit"), data = frs_long,
                       file = paste0(cache_dir, "/ALSFRSdimensions1_time.rds"), cores=4)
fit_dimensions1

vc <- VarCorr(fit_dimensions1)
corrs <- vc[[1]][[2]]


mcmc_trace(fit_dimensions1, pars = c("cor_ID__questionALSFRS_Q1__questionALSFRS_Q11",
                              "b_Intercept[1]",
                              "cor_ID__questionALSFRS_Q6__questionALSFRS_Q7",
                              "b_questionALSFRS_Q4") )

# Posterior predictive check
mod1_pp1 <- pp_check(fit_dimensions1, nsamples = 200)
print(mod1_pp1)




# fit second model
fit_multivariate <- brm(
    mvbind(Q01, Q02, Q03, Q04, Q05, Q06,
           Q07, Q08, Q09, Q10, Q11, Q12) ~ 1 + alsfrs_dly_mnths + (1 | p | ID),
                        data = df_frs,
                        family = cumulative("logit"),
                        file = paste0(cache_dir, "/ALSFRSmultivariate.rds"), cores = 4)
fit_multivariate




# fit third model
fit_cratio <- brm(answer ~ alsfrs_dly_mnths + cs(question) + (1 + question | ID),
                  family = cratio("logit"), data = frs_long,
                  file = paste0(cache_dir, "/ALSFRScratio.rds"), cores = 4)
fit_cratio




# fit model four
cumulative_logit_cs5 <- custom_family(
    "cumulative_logit_cs5", dpars = c("mu", "threshA", "widthB", "widthC", "widthD"),
    links = c("identity", "identity", "log", "log", "log"), lb = c(NA, NA, 0, 0, 0),
    type = "int"
)

stan_funs <- "
  // The dummy parameter is required, probably due to a bug in brms (I don't need it with fewer than 4 distributional params...)
  real cumulative_logit_cs5_lpmf(int y, real mu, real thresh1, real width2, real width3, real width4, real dummy) {
    vector[4] params = to_vector({thresh1, width2, width3, width4});
    vector[4] thresholds = cumulative_sum(params);
    return ordered_logistic_lpmf(y | mu, thresholds);
  }
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

# At least somewhat constraining priors are needed for the fit to work
prior_custom <-  c(set_prior("student_t(3, 0, 2.5)", class = "b", dpar = "threshA"),
                   set_prior("normal(0.5, 2)", class = "b", dpar = "widthB"),
                   set_prior("normal(0.5, 2)", class = "b", dpar = "widthC"),
                   set_prior("normal(0.5, 2)", class = "b", dpar = "widthD")
)


fit_custom <- brm(bf(answer ~ 0 + alsfrs_dly_mnths  + (1 + question | ID),
                     threshA ~ 0 + Intercept + question, # Avoiding centering of the intercept
                     widthB  ~ 0 + Intercept + question,
                     widthC  ~ 0 + Intercept + question,
                     widthD  ~ 0 + Intercept + question),
                  family = cumulative_logit_cs5, data = frs_long, stanvars = stanvars,
                  prior = prior_custom,
                  file = paste0(cache_dir, "/ALSFRScustom.rds"),
                  cores=4)


fit_custom





