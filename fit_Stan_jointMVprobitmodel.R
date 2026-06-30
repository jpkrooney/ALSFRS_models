library(tidyverse)
library(brms)
library(ggdist)
library(rstan)
library(patchwork)
#library(mori)
library(furrr)
library(nutpieR)
library(bayesplot)
library(splines2)



#options(mc.cores = 8, brms.backend = "cmdstanr")
cache_dir <- here::here("local_temp_data", "simulation_test")
if(!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
}
source(here::here("R", "pp_checks.R"))
source(here::here("R", "mv_probit.R"))
source(here::here("R", "simulate_data.R"))
source(here::here("R", "sampling_parallel.R"))
source(here::here("R", "brm_parallel.R"))
source(here::here("R", "nutpieR_to_rstan_obj.R"))
source(here::here("R", "effect_size_funcs.R"))
source("R/calc_basehaz_frame.R")

theme_set(cowplot::theme_cowplot())
irldat <- read_csv(here::here("Data/longdata.csv"))
iddat <- read_csv(here::here("Data/casedata.csv"))

irldat_full <- irldat %>%
    rename("subject_id" = "ID",
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
           "Q12" = "ALSFRS_Q12",
           "alsfrs_total" = "Total") %>%
    dplyr::select(subject_id, alsfrs_dly_mnths, Q01, Q02, Q03, Q04, Q05, Q06,
                  Q07, Q08, Q09, Q10, Q11, Q12, alsfrs_total)

question_col_names_irldat <- sort(names(irldat_full)[grepl("^[QR]", names(irldat_full))])
question_col_names <- sprintf("Q%02d", 1:12)
question_col_names_translation <- sprintf("Q%02d", 1:12)
names(question_col_names_translation) <- question_col_names_irldat

irldat_standardized <- irldat_full




M <- 50
n_subjects_per_group <- 50
set.seed(99)

questions_to_fit <- question_col_names#[1:4]

sv <- make_stanvars_mv_probit_bgoodri(questions_to_fit)

eff_probs <- rep(c(0.5), times = M)


irldat_formula <- as.formula(paste0(
    "mvbind(", paste0(questions_to_fit, collapse = ", "),
    ") ~ 1 + time2 * group + (1 | p | subject_id)"
))



# make data
ls_dat <- lapply(1:length(eff_probs), function(m)
    simulate_data_from_registry(irldat_standardized |>
                                    group_by(subject_id) |>
                                    filter(first(alsfrs_dly_mnths < 30)),
                                max_duration = 30,
                                min_measurements_per_subject = 4,
                                max_measurements_per_subject = 8,
                                n_subjects_per_group = n_subjects_per_group,
                                effect_prob = eff_probs[m]) |>
        group_by(subject_id) |>
        mutate(time2 = alsfrs_dly_mnths - min(alsfrs_dly_mnths)))

# Remove any longitudinal measurements occuring after survival event
ls_dat <- ls_dat |>
    map(\(dat)
        dat |>
            left_join(iddat |>
                          dplyr::select(ID, dx_surv_mnths, vital_status, age_dx),
                      by = join_by(old_id.x == ID))) #|>
#filter(dx_surv_mnths >= time2))


ls_dat <- ls_dat |>
    map(\(dat) dat |>
            mutate(across(all_of(question_col_names),
                          \(x) { x + 1 })))


stancode1 <- brms::make_stancode(bf(irldat_formula, family = empty_cumulative()) +
                                     set_rescor(FALSE), stanvars = sv, adapt_delta = 0.95, init = 0.0,
                                 data = ls_dat[[1]])
writeLines(stancode1, "tempfile.stan")
# manual edits


stanmodel1 <- stan_model(model_code = stancode1)
model1 <- nutpie_compile_model(code = stancode1)


# load joint model which has been manually edited
stancode_jm1 <- readLines("stan/stan_mspline_jm_currentvalue_assoc.stan")
stanmodel_jm1 <- stan_model(model_code = stancode_jm1)
modeljm1 <- nutpie_compile_model(code = stancode_jm1 )


# add event data needed for the joint model
ls_standat_jm1 <- ls_dat |>
    map(\(dat){
        standat <- brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                           set_rescor(FALSE), stanvars = sv, adapt_delta = 0.95, init = 0.0,
                                       data = dat)
        dat_events <- dat |>
            group_by(subject_id) |>
            filter(time2 == max(time2)) |>
            #left_join(iddat |>
            #              dplyr::select(ID, dx_surv_mnths, vital_status),
            #          by = join_by(old_id.x == ID)) |>
            mutate(event = ifelse(vital_status == "Alive", 0, 1))

        c(standat,
          calc_basehaz_frame(time = dat_events$dx_surv_mnths, n_iknots = 8),
          event_times = list(dat_events$dx_surv_mnths),
          event = list(dat_events$event),
          #Xsurv = list(ifelse(dat_events$group == "Control", 0, 1)),
          Xsurv = list(dat_events$age_dx)
        )
    })


#### Fit a test model ####
test <- nutpie_sample(modeljm1,
                      data = ls_standat_jm1[[3]],
                      num_warmup = 500,
                      num_draws = 250,
                      num_chains = 4,
                      cores = 4,
                      refresh = 50,
                      seed = 604712,
                      init_mean = 0.1,
                      pars = "u",
                      include = FALSE)
# convert nutpie fit to brms
testbrms <- nutpie_to_brms(irldat_formula, ls_dat[[1]], model = stanmodel_jm1,
                           fit = test, family = empty_cumulative(),
                           stanvars = sv)

# histogram of effect size
testbrms |>
    effect_from_sim_study( 0, 12, timevar = "time2",
                           method = "expected", summation ="sum_allQ", model = "mv-probit") |>
    hist(100)


# get posterior draws
testdraws <- as_draws_df(test)

testdraws |>
    select(starts_with('beta_surv')) |> unlist() |> hist(100)
testdraws |>
    select(starts_with('assoc')) |> summary()
testdraws |>
    select(starts_with('r_1_Q01')) |> summary()
testdraws |>
    select("eta0") |> unlist() |> hist(100)

mcmc_pairs(testdraws, pars = c("beta_surv", "wb_shape", "wb_scale", "assoc[1]"))

mcmc_pairs(testdraws, pars = c("beta_surv", "eta0", "w[1]", "w[2]", "w[3]", "w[4]", "w[5]"))
mcmc_pairs(testdraws, pars = c("beta_surv", "eta0", "assoc[1]", "theta[5]"))
mcmc_pairs(testdraws, pars = c("beta_surv", "assoc[5]", "theta[5]", "assoc[7]", "assoc[8]"))
mcmc_pairs(testdraws, pars = c("beta_surv", "eta0", "theta[3]", "w[3]", "assoc[1]", "theta[5]"))
mcmc_pairs(testdraws, pars = c("assoc[1]", "r_1_Q01_1[1]", #"sd_1[1]",
                               "Intercept_Q11[1]", "b_Q11_Intercept[2]", "w[3]"))


mcmc_pairs(testdraws, pars = c("beta_surv", "eta0", "assoc[1]", "lpred_q01"))



chk <- testdraws |>
    select(starts_with('lpred_')) |>
    slice(1) |> t() |>
    data.frame(rep(ls_standat_jm1[[1]]$event_times, 12),
               m = factor(rep(seq_len(100), each = 12))) |>
    rownames_to_column() |>
    mutate(Q = str_remove(rowname, "lpred_"),
           Q = str_remove(Q, "\\[[0-9]+\\]"))

names(chk) <- c("rowname", "lpred", "time", "m", "Q")

ggplot(chk, aes(x = time, y = lpred, group = m)) + geom_line() +
    facet_wrap(~Q, scales = "free_y")


testdraws |>
    select(starts_with('lpred_')) |>
    mutate(m = 1:n()) |>
    pivot_longer( -m) |>
    mutate(Q = stringr::str_match(name, "^lpred_q([[:digit:]]+)\\[([[:digit:]]+)\\]$")[, 2],
           id = as.numeric(stringr::str_match(name, "^lpred_q([[:digit:]]+)\\[([[:digit:]]+)\\]$")[, 3])) |>
    pivot_wider(id_cols = c(m, id), names_from = Q, values_from = value) |>
    mutate(total = `01` + `02` + `03` + `04` + `05` + `06` +
               `07` + `08` + `09` + `10` + `11` + `12`) |>
    group_by(id) |>
    summarise(median = median(total),
              q2_5 = quantile(total, probs = 0.025),
              q97_5 = quantile(total, probs = 0.975)) |>
    bind_cols(times = ls_standat_jm1[[1]]$event_times) |>
    ggplot(aes(x = times)) +
    geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.5) +
    geom_line(aes(y = median))









