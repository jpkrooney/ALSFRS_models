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
source(here::here("R", "pp_checks.R"))
source(here::here("R", "mv_probit_het.R"))
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




M <- 200
n_subjects_per_group <- 50
set.seed(99)

questions_to_fit <- question_col_names#[1:4]

sv <- make_stanvars_mv_probit_bgoodri(questions_to_fit)

eff_probs <- rep(c(0.5, 0.85), times = M)


irldat_formula <- as.formula(paste0(
    #"mvbind(", paste0(questions_to_fit, collapse = ", "),
    "mvbind(", paste0(questions_to_fit, collapse = ", "),
    ") ~ 1 + time2 * group + (time2 | p | subject_id)"
))



ls_dat <- lapply(1:length(eff_probs), function(m)
    simulate_data_from_registry(irldat_standardized |>
                                    group_by(subject_id) |>
                                    filter(first(alsfrs_dly_mnths < 30)),
                                max_duration = 30,
                                min_measurements_per_subject = 2,
                                max_measurements_per_subject = 8,
                                n_subjects_per_group = n_subjects_per_group,
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

stanmodel1 <- stan_model(model_code = stancode1)


# format data for stan
ls_standat1 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv, adapt_delta = 0.95, init = 0.0,
                            data = dat))

#sink(file = "modelfornutpie.stan", )
model1 <- nutpie_compile_model(code = stancode1)


# Sample
plan(multisession, workers=10)

brmsnutfits1 <- list(ls_dat,
                 ls_standat1) |>
    future_pmap( \(dat, standat){
        fit <- nutpie_sample(model1,
                      data = standat,
                      num_warmup = 500,
                      num_draws = 250,
                      num_chains = 4,
                      cores = 1,
                      refresh = 0,
                      seed = 60471,
                      init = 0.1,
                      pars = c("u", "r_1", "z_1", "L_1", "sd_1"),
                      include = FALSE)
        nutpie_to_brms(irldat_formula, dat, model = stanmodel1,
                       fit = fit, family = empty_cumulative(),
                       stanvars = sv)
        }, .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

#### convert fits to brms and save results ####
# set filename pars
stem1 <- paste0("_nutpiefit_mvprob_het", n_subjects_per_group)
#brmsfits1 <- lapply(1: length(eff_probs), function(i) {
#    nutpie_to_brms(irldat_formula, ls_dat[[i]], model = stanmodel1,
#                   fit = nutfits1[[i]], family = empty_cumulative(),
#                   stanvars = sv)
#})
saveRDS(brmsnutfits1, paste0(cache_dir, "/", Sys.Date(), stem1, "mvp_defaultpriors_", M,".RDS"))
rm(nutfits1) ; gc()
#rm(brmsfits1) ; gc()


brmsfits0 <- readRDS(paste0(cache_dir, "/", "2026-07-12", "_nutpiefit_mvprob_het50", "mvp_defaultpriors_100", ".RDS"))


plan(multisession, workers=10)
summaries_mvprobhet50 <- brmsfits0 |>
    future_map(\(fit) effect_from_sim_study(fit, 0, 12, timevar = "time2",
                                            method = "expected" , model = "mv-probit"),
               .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

rm(brmsfits0) ; gc()

brmsfits1 <- readRDS(paste0(cache_dir, "/", "2026-07-13", "_nutpiefit_mvprob_het70", "mvp_defaultpriors_", M,".RDS"))

plan(multisession, workers=10)
summaries_mvprobhet70 <- brmsfits1 |>
    future_map(\(fit) effect_from_sim_study(fit, 0, 12, timevar = "time2",
                                            method = "expected" , model = "mv-probit"),
               .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

rm(brmsfits1) ; gc()




summ_mvphet <- map_df(1:length(eff_probs),
                   ~ data.frame(id = .x, type = "mv-probit",
                                eff_prob = eff_probs[[.x]], estimate = summaries_mvprobhet[[.x]],
                                lab = "1 | p | id"))


summs <- summ_mvphet |>
    summarise(.by = c(id, type, lab),
              eff_prob = unique(eff_prob),
              q_025 = quantile(estimate, 0.025),
              mean = mean(estimate),
              q_975 = quantile(estimate, 0.975),
              flag1 = ifelse(eff_prob == 0.5 & q_025 > 0 | q_975 < 0, 1, 0),
              flag2 = ifelse(q_025 > 0, 1, 0)) |>
    summarise(.by = c(type, lab, eff_prob),
              type1error = mean(flag1, na.rm = TRUE),
              power = mean(flag2)) |>
    arrange(eff_prob, type)
summs

plot_mv_probit <- summ_mvphet %>%
    mutate(eff_prob = factor(eff_prob)) %>%
    ggplot(aes(x =  id, y = estimate)) +
    geom_hline(yintercept = 0) +
    stat_pointinterval(position = position_dodge(width = 0.3)) +
    #    expand_limits(y = shared_range) +
    ggtitle(label = "Difference in expected value by group",
            subtitle = paste0(M, " simulations,", n_subjects_per_group, " subjects per group")) +
    facet_wrap(~ eff_prob + type + lab, ncol = 2)#, scales = "free_y")

plot_mv_probit




