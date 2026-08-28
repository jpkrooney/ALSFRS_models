library(tidyverse)
library(brms)
library(ggdist)
library(rstan)
library(patchwork)
#library(mori)
library(furrr)
library(nutpieR)

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




M <- 2
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



##plan(multisession, workers=8)
#brmsfits <- future_map(1:length(eff_probs), \(i) {
brmsfits <- lapply(1: length(eff_probs), function(i) {
        # currently need to fit empty obejcts for each fit - this is slow
        #emptyrstan <- rstan::sampling(stanmodel, data = ls_standat[[i]],
        #                          chains = 4, iter = 1, cores = 4, refresh = 0,
        #                          init = 0.0, control = list(adapt_delta = 0.95))
        emptybrms <- brm(bf(irldat_formula, family = empty_cumulative()) +
                         set_rescor(FALSE), stanvars = sv, adapt_delta = 0.95, init = 0.1,
                     data = ls_dat[[i]], empty = TRUE)
        #
        rstanobj <- nutpieR_to_rstan_obj(nutfits[[i]], stanmodel1)
        brmsfit <- emptybrms
        brmsfit$fit <- rstanobj
        brmsfit <- rename_pars(brmsfit)
        brmsfit}) #, .progress = TRUE, .options = furrr_options(seed = 123))

saveRDS(brmsfits, paste0(cache_dir, "/", "20260402_nutpiemvprobit_min3scores.RDS"))


###### Reload fits for plotting ######
brmsfits1 <- readRDS(paste0(cache_dir, "/", "20260402_nutpiemvprobit_fits_18mnths_n50.RDS"))
brmsfits2 <- readRDS(paste0(cache_dir, "/", "20260402_nutpiemvprobit_fits_24mnths_n50.RDS"))
brmsfits3 <- readRDS(paste0(cache_dir, "/", "20260402_nutpiemvprobit_fits_30mnths_n50.RDS"))
brmsfits4 <- readRDS(paste0(cache_dir, "/", "20260402_nutpiemvprobit_fits_18mnths_n50_no_p.RDS"))
brmsfits5 <- readRDS(paste0(cache_dir, "/", "20260402_nutpiemvprobit_fits_24mnths_n50_no_p.RDS"))
brmsfits6 <- readRDS(paste0(cache_dir, "/", "20260402_nutpiemvprobit_fits_30mnths_n50_no_p.RDS"))


#### define an effect function ####
effect_from_sim_study <- function(fit, time1, time2) {
    predict_data <- data.frame(group = c("Control", "Control", "Treatment", "Treatment"),
                               alsfrs_dly_mnths = c(time1, time2, time1, time2))

    pred <- posterior_epred(fit, newdata = predict_data, re_formula = NA)

    expected_sum <- function(x) {
        col_values <- as.integer(colnames(x))
        multiplied <- sweep(x, MARGIN = 2, STATS = col_values, FUN = "*")
        rowSums(multiplied)
    }

    start_control <- expected_sum(pred[,1,])
    end_control <- expected_sum(pred[,2,])
    start_treatment <- expected_sum(pred[,3,])
    end_treatment <- expected_sum(pred[,4,])

    (end_treatment - start_treatment) - (end_control - start_control)
}


#### Summarise effects ####
summaries_mv_probit1 <- brmsfits1 |>
    map(\(fit) effect_from_sim_study(fit))
summaries_mv_probit2 <- brmsfits2 |>
    map(\(fit) effect_from_sim_study(fit))
summaries_mv_probit3 <- brmsfits3 |>
    map(\(fit) effect_from_sim_study(fit))
summaries_mv_probit4 <- brmsfits4 |>
    map(\(fit) effect_from_sim_study(fit))
summaries_mv_probit5 <- brmsfits5 |>
    map(\(fit) effect_from_sim_study(fit))
summaries_mv_probit6 <- brmsfits6 |>
    map(\(fit) effect_from_sim_study(fit))

all_effects_mv_probit1 <- purrr::map_df(1:length(eff_probs),
                                       ~ data.frame(id = .x, type = "mv_probit",
                                                    eff_prob = 0.5, estimate = summaries_mv_probit1[[.x]],
                                                    maxdelay = "18", lab = "p"))
all_effects_mv_probit2 <- purrr::map_df(1:length(eff_probs),
                                        ~ data.frame(id = .x, type = "mv_probit",
                                                     eff_prob = 0.5, estimate = summaries_mv_probit2[[.x]],
                                                     maxdelay = "24", lab = "p"))
all_effects_mv_probit3 <- purrr::map_df(1:length(eff_probs),
                                        ~ data.frame(id = .x, type = "mv_probit",
                                                     eff_prob = 0.5, estimate = summaries_mv_probit3[[.x]],
                                                     maxdelay = "30", lab = "p"))
all_effects_mv_probit4 <- purrr::map_df(1:length(eff_probs),
                                        ~ data.frame(id = .x, type = "mv_probit",
                                                     eff_prob = 0.5, estimate = summaries_mv_probit4[[.x]],
                                                     maxdelay = "18", lab = "no p"))
all_effects_mv_probit5 <- purrr::map_df(1:length(eff_probs),
                                        ~ data.frame(id = .x, type = "mv_probit",
                                                     eff_prob = 0.5, estimate = summaries_mv_probit5[[.x]],
                                                     maxdelay = "24", lab = "no p"))
all_effects_mv_probit6 <- purrr::map_df(1:length(eff_probs),
                                        ~ data.frame(id = .x, type = "mv_probit",
                                                     eff_prob = 0.5, estimate = summaries_mv_probit6[[.x]],
                                                     maxdelay = "30", lab = "no p"))

all_effects_mv_probit <- bind_rows(all_effects_mv_probit1,
                                   all_effects_mv_probit2,
                                   all_effects_mv_probit3,
                                   all_effects_mv_probit4,
                                   all_effects_mv_probit5,
                                   all_effects_mv_probit6)


summs <- all_effects_mv_probit |>
    group_by(id, maxdelay, lab) |>
    summarise(#eff_prob = unique(eff_prob),
              q_025 = quantile(estimate, 0.025),
              mean = mean(estimate),
              q_975 = quantile(estimate, 0.975),
              flag1 = ifelse(q_025 < 0 & q_975 >0, 1, 0)) |>
    group_by(maxdelay, lab) |>
    summarise(coverage = mean(flag1, na.rm = TRUE)) |>
    arrange(lab, maxdelay)
summs


plot_mv_probit <- all_effects_mv_probit %>%
    mutate(eff_prob = factor(eff_prob),
           maxdelay = factor(maxdelay)) %>%
    ggplot(aes(x =  id, y = estimate)) +
    geom_hline(yintercept = 0) +
    stat_pointinterval(position = position_dodge(width = 0.3)) +
    #    expand_limits(y = shared_range) +
    ggtitle("MV probit") + facet_wrap(~ maxdelay + lab, ncol = 2)

plot_mv_probit

