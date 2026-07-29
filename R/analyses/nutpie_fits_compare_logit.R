library(tidyverse)
library(brms)
library(ggdist)
library(rstan)
library(patchwork)
library(mori)
library(furrr)
library(nutpieR)

#options(mc.cores = 8, brms.backend = "cmdstanr")
cache_dir <- here::here("local_temp_data", "simulation_test")
if(!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
}
source(here::here("R/utils", "pp_checks.R"))
source(here::here("R/utils", "mv_probit.R"))
source(here::here("R/utils", "simulate_data.R"))
source(here::here("R/utils", "sampling_parallel.R"))
source(here::here("R/utils", "brm_parallel.R"))
source(here::here("R/utils", "nutpieR_to_rstan_obj.R"))
source(here::here("R/utils", "effect_size_funcs.R"))

theme_set(cowplot::theme_cowplot())
irldat <- read_csv(here::here("Data/longdata.csv"))

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




M <- 100
set.seed(99)

questions_to_fit <- question_col_names
sv <- make_stanvars_mv_probit_bgoodri(questions_to_fit)

eff_probs <- rep(c(0.5), times = M)


irldat_formula <- as.formula(paste0(
    "mvbind(", paste0(questions_to_fit, collapse = ", "),
    ") ~ 1 + alsfrs_dly_mnths*group + (1 | p | subject_id)"
))

irldat_formula2 <- as.formula(paste0(
    "mvbind(", paste0(questions_to_fit, collapse = ", "),
    ") ~ 1 + time2*group + (1 | p | subject_id)"
))


#irldat_formula <- as.formula(paste0(
#    "mvbind(", paste0(questions_to_fit, collapse = ", "),
#    ") ~ 1 + alsfrs_dly_mnths*group + (1 | subject_id)"
#))




ls_dat <- lapply(1:length(eff_probs), function(m)
    simulate_data_from_registry(irldat_standardized |>
                                    group_by(subject_id) |>
                                    filter(first(alsfrs_dly_mnths < 30)),
                            max_duration = 30,
                            min_measurements_per_subject = 3,
                            max_measurements_per_subject = 8,
                            n_subjects_per_group = 50,
                            effect_prob = eff_probs[m]) |>
        group_by(subject_id) |>
        mutate(time2 = alsfrs_dly_mnths - min(alsfrs_dly_mnths)))


ls_dat <- ls_dat |>
    map(\(dat) dat |>
            mutate(across(all_of(question_col_names),
                  \(x) { x + 1 })))


stancode1 <- brms::make_stancode(bf(irldat_formula, family = empty_cumulative()) +
                                    set_rescor(FALSE), stanvars = sv, adapt_delta = 0.95, init = 0.0,
                                data = ls_dat[[1]])
stancode2 <- brms::make_stancode(bf(irldat_formula2, family = empty_cumulative()) +
                                     set_rescor(FALSE), stanvars = sv, adapt_delta = 0.95, init = 0.0,
                                 data = ls_dat[[1]])
# logit fits
stancode3 <- brms::make_stancode(bf(irldat_formula, family = cumulative()) +
                                    set_rescor(FALSE), adapt_delta = 0.95, init = 0.0,
                                data = ls_dat[[1]])
stancode4 <- brms::make_stancode(bf(irldat_formula2, family = cumulative()) +
                                     set_rescor(FALSE), adapt_delta = 0.95, init = 0.0,
                                 data = ls_dat[[1]])

stanmodel1 <- stan_model(model_code = stancode1)
stanmodel2 <- stan_model(model_code = stancode2)
stanmodel3 <- stan_model(model_code = stancode3)
stanmodel4 <- stan_model(model_code = stancode4)


# format data for stan
ls_standat1 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv, adapt_delta = 0.95, init = 0.0,
                            data = dat))
ls_standat2 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula2, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv, adapt_delta = 0.95, init = 0.0,
                            data = dat))
ls_standat3 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = cumulative()) +
                                set_rescor(FALSE), adapt_delta = 0.95, init = 0.0,
                            data = dat))
ls_standat4 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula2, family = cumulative()) +
                                set_rescor(FALSE), adapt_delta = 0.95, init = 0.0,
                            data = dat))

#sink(file = "modelfornutpie.stan", )
model1 <- nutpie_compile_model(code = stancode1)
model2 <- nutpie_compile_model(code = stancode2)
model3 <- nutpie_compile_model(code = stancode3)
model4 <- nutpie_compile_model(code = stancode4)

# Sample
plan(multisession, workers=3)
nutfits1 <- ls_standat1 |>
    future_map( \(dat)
         nutpie_sample(model1,
                       data = dat,
                       num_warmup = 500,
                       num_draws = 250,
                       num_chains = 4,
                       cores = 4,
                       refresh = 0,
                       seed = 60471,
                       init_mean = 0,
                       pars = "u",
                       include = FALSE), #|>
             #trim_pars(pars = "u"), # exclude nuisance pars to save memory
         .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

plan(multisession, workers=3)
nutfits2 <- ls_standat2 |>
    future_map( \(dat)
                nutpie_sample(model2,
                              data = dat,
                              num_warmup = 500,
                              num_draws = 250,
                              num_chains = 4,
                              cores = 4,
                              refresh = 0,
                              seed = 60471,
                              init_mean = 0,
                              pars = "u",
                              include = FALSE), #|>
                #trim_pars(pars = "u"), # exclude nuisance pars to save memory
                .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

plan(multisession, workers=3)
nutfits3 <- ls_standat3 |>
    future_map( \(dat)
                nutpie_sample(model3,
                              data = dat,
                              num_warmup = 500,
                              num_draws = 250,
                              num_chains = 4,
                              cores = 4,
                              refresh = 0,
                              seed = 60471,
                              init_mean = 0),
                              #pars = "u",
                              #include = FALSE), #|>
                #trim_pars(pars = "u"), # exclude nuisance pars to save memory
                .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

plan(multisession, workers=3)
nutfits4 <- ls_standat4 |>
    future_map( \(dat)
                nutpie_sample(model4,
                              data = dat,
                              num_warmup = 500,
                              num_draws = 250,
                              num_chains = 4,
                              cores = 4,
                              refresh = 0,
                              seed = 60471,
                              init_mean = 0),
                              #pars = "u",
                              #include = FALSE), #|>
                #trim_pars(pars = "u"), # exclude nuisance pars to save memory
                .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)



#### convert fits to brms and save results ####
# set filename pars
stem1 <- "_nutpiefit_"
brmsfits1 <- lapply(1: length(eff_probs), function(i) {
    nutpie_to_brms(irldat_formula, ls_dat[[i]], model = stanmodel1,
                   fit = nutfits1[[i]], family = empty_cumulative(),
                   stanvars = sv)
    })
saveRDS(brmsfits1, paste0(cache_dir, "/", Sys.Date(), stem1, "_mvprobit1_m", M,".RDS"))

brmsfits2 <- lapply(1: length(eff_probs), function(i) {
    nutpie_to_brms(irldat_formula2, ls_dat[[i]], model = stanmodel2,
                   fit = nutfits2[[i]], family = empty_cumulative(),
                   stanvars = sv)
})
saveRDS(brmsfits2, paste0(cache_dir, "/", Sys.Date(), stem1, "_mvprobit2_m", M,".RDS"))

brmsfits3 <- lapply(1: length(eff_probs), function(i) {
    nutpie_to_brms(irldat_formula, ls_dat[[i]], model = stanmodel3,
                   fit = nutfits3[[i]], family = cumulative(),
                   stanvars = NULL)
})
saveRDS(brmsfits3, paste0(cache_dir, "/", Sys.Date(), stem1, "_logit1_m", M,".RDS"))

brmsfits4 <- lapply(1: length(eff_probs), function(i) {
    nutpie_to_brms(irldat_formula2, ls_dat[[i]], model = stanmodel4,
                   fit = nutfits4[[i]], family = cumulative(),
                   stanvars = NULL)
})
saveRDS(brmsfits4, paste0(cache_dir, "/", Sys.Date(), stem1, "_logit2_m", M,".RDS"))

# reload brmsfits



#### Summarise effects ####

brmsfits1 <- share(brmsfits1)
brmsfits2 <- share(brmsfits2)
brmsfits3 <- share(brmsfits3)
brmsfits4 <- share(brmsfits4)

method1 <- "predicted"
plan(multisession, workers=10)
summaries_mv_probit1 <- brmsfits1 |>
    future_map(\(fit) effect_from_sim_study(fit, 0, 12, timevar = "alsfrs_dly_mnths",
                                     method = method1 , model = "mv-probit"),
               .progress = TRUE, .options = furrr_options(seed = 123))
summaries_mv_probit2 <- brmsfits2 |>
    future_map(\(fit) effect_from_sim_study(fit, 0, 12, timevar = "time2",
                                     method = method1 , model = "mv-probit"),
               .progress = TRUE, .options = furrr_options(seed = 123))
summaries_logit1 <- brmsfits3 |>
    future_map(\(fit) effect_from_sim_study(fit, 0, 12, timevar = "alsfrs_dly_mnths",
                                     method = method1 , model = "logit"),
               .progress = TRUE, .options = furrr_options(seed = 123))
summaries_logit2 <- brmsfits4 |>
    future_map(\(fit) effect_from_sim_study(fit, 0, 12, timevar = "time2",
                                     method = method1 , model = "logit"),
               .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)


summ_mvprob1 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "mv-probit",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_mv_probit1[[.x]],
                                    lab = "real time"))
summ_mvprob2 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "mv-probit",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_mv_probit2[[.x]],
                                    lab = "relative time"))
summ_logit1 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "logit",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_logit1[[.x]],
                                    lab = "real time"))
summ_logit2 <- map_df(1:length(eff_probs),
                      ~ data.frame(id = .x, type = "logit",
                                   eff_prob = eff_probs[[.x]], estimate = summaries_logit2[[.x]],
                                   lab = "relative time"))


all_effects_mv_probit <- bind_rows(summ_mvprob1,
                                   summ_mvprob2,
                                   summ_logit1,
                                   summ_logit2)


summs <- all_effects_mv_probit |>
    summarise(.by = c(id, type, lab),
              #eff_prob = unique(eff_prob),
              q_025 = quantile(estimate, 0.025),
              mean = mean(estimate),
              q_975 = quantile(estimate, 0.975),
              flag1 = ifelse(q_025 > 0 | q_975 < 0, 1, 0)) |>
    summarise(.by = c(type, lab),
              type1error = mean(flag1, na.rm = TRUE)) |>
    arrange(lab)
summs

plot_mv_probit <- all_effects_mv_probit %>%
    mutate(eff_prob = factor(eff_prob)) %>%
    ggplot(aes(x =  id, y = estimate)) +
    geom_hline(yintercept = 0) +
    stat_pointinterval(position = position_dodge(width = 0.3)) +
    #    expand_limits(y = shared_range) +
    ggtitle("Full predictions") + facet_wrap(~ type + lab, ncol = 2)

plot_mv_probit

