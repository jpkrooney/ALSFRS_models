library(tidyverse)
library(glue)
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
source(here::here("R", "define_empty_cumulative_family.R"))
source(here::here("R", "mv_stanchunks.R"))
source(here::here("R", "mv_probit_studt_heterosked.R"))
source(here::here("R", "simulate_data.R"))
source(here::here("R", "sampling_parallel.R"))
source(here::here("R", "brm_parallel.R"))
source(here::here("R", "nutpieR_to_rstan_obj.R"))
source(here::here("R", "effect_size_funcs.R"))
source(here::here("R", "post_epred_hetsked.R"))

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




M <- 400
n_subjects_per_group <- 50
set.seed(99)

questions_to_fit <- question_col_names[1:4]

# vary the prior on d.f. of student-t distribution
sv1 <- make_stanvars_mv_probit_bgoodri(questions_to_fit, hetsked = FALSE)
sv2 <- make_stanvars_mv_probit_bgoodri(questions_to_fit, hetsked = TRUE)
sv3 <- make_stanvars_mv_studt(questions_to_fit, nu_as_par = FALSE, nu = 3, hetsked = FALSE)
sv4 <- make_stanvars_mv_studt(questions_to_fit, nu_as_par = TRUE, nu_beta = 0.1, hetsked = FALSE)
sv5 <- make_stanvars_mv_studt(questions_to_fit, nu_as_par = TRUE, nu_beta = 0.9, hetsked = FALSE)
sv6 <- make_stanvars_mv_studt(questions_to_fit, nu_as_par = FALSE, nu = 3, hetsked = TRUE)



eff_probs <- rep(c(0.5, 0.75), times = M)


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
                                effect_prob = eff_probs[m],
                                include_cols = questions_to_fit) |>
        group_by(subject_id) |>
        mutate(time2 = alsfrs_dly_mnths - min(alsfrs_dly_mnths)))

# save data
saveRDS(ls_dat, paste0(cache_dir, "/sim_irldata_4Q.RDS"))

ls_dat |>
    map(\(dat)
        dat |>
            group_by(subject_id) |>
            filter(visit_id == max(visit_id)) |>
            ungroup() |>
            select(visit_id)) |>
    unlist() |>
    hist()





ls_dat <- ls_dat |>
    map(\(dat) dat |>
            mutate(across(all_of(question_col_names),
                          \(x) { x + 1 })))


stancode1 <- brms::make_stancode(bf(irldat_formula, family = empty_cumulative()) +
                                     set_rescor(FALSE), stanvars = sv1,
                                 data = ls_dat[[1]])
stancode2 <- brms::make_stancode(bf(irldat_formula, family = empty_cumulative()) +
                                     set_rescor(FALSE), stanvars = sv2,
                                 data = ls_dat[[1]])
stancode3 <- brms::make_stancode(bf(irldat_formula, family = empty_cumulative()) +
                                     set_rescor(FALSE), stanvars = sv3,
                                 data = ls_dat[[1]])
stancode4 <- brms::make_stancode(bf(irldat_formula, family = empty_cumulative()) +
                                     set_rescor(FALSE), stanvars = sv4,
                                 data = ls_dat[[1]])
stancode5 <- brms::make_stancode(bf(irldat_formula, family = empty_cumulative()) +
                                     set_rescor(FALSE), stanvars = sv5,
                                 data = ls_dat[[1]])
stancode6 <- brms::make_stancode(bf(irldat_formula, family = empty_cumulative()) +
                                     set_rescor(FALSE), stanvars = sv6,
                                 data = ls_dat[[1]])

stanmodel1 <- stan_model(model_code = stancode1)
stanmodel2 <- stan_model(model_code = stancode2)
stanmodel3 <- stan_model(model_code = stancode3)
stanmodel4 <- stan_model(model_code = stancode4)
stanmodel5 <- stan_model(model_code = stancode5)
stanmodel6 <- stan_model(model_code = stancode6)


# format data for stan
ls_standat1 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv1,
                            data = dat))
ls_standat2 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv2,
                            data = dat))
ls_standat3 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv3,
                            data = dat))
ls_standat4 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv4,
                            data = dat))
ls_standat5 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv5,
                            data = dat))
ls_standat6 <- ls_dat |>
    map(\(dat)
        brms::make_standata(bf(irldat_formula, family = empty_cumulative()) +
                                set_rescor(FALSE), stanvars = sv6,
                            data = dat))


#sink(file = "modelfornutpie.stan", )
model1 <- nutpie_compile_model(code = stancode1)
model2 <- nutpie_compile_model(code = stancode2)
model3 <- nutpie_compile_model(code = stancode3)
model4 <- nutpie_compile_model(code = stancode4)
model5 <- nutpie_compile_model(code = stancode5)
model6 <- nutpie_compile_model(code = stancode6)


# Sample with first prior
plan(multisession, workers=10)

brmsnutfits1 <- list(ls_dat,
                     ls_standat1) |>
    future_pmap( \(dat, standat){
        fit <- nutpie_sample(model1,
                             data = standat,
                             num_warmup = 500,
                             num_draws = 250,
                             num_chains = 3,
                             cores = 1,
                             refresh = 0,
                             seed = 60471,
                             init = 0.01,
                             pars = c("u", "r_1", "z_1", "L_1", "sd_1"),
                             include = FALSE)
        nutpie_to_brms(irldat_formula, dat, model = stanmodel1,
                       fit = fit, family = empty_cumulative(),
                       stanvars = sv1)
    }, .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

#### convert fits to brms and save results ####
# set filename pars
stem1 <- paste0("_compare_copula_Q4_n", n_subjects_per_group)
saveRDS(brmsnutfits1, paste0(cache_dir, "/", Sys.Date(), stem1, "_mvprob_", M,".RDS"))
rm(brmsnutfits1) ; gc()




# Sample with second prior
plan(multisession, workers=10)

brmsnutfits2 <- list(ls_dat,
                     ls_standat2) |>
    future_pmap( \(dat, standat){
        fit <- nutpie_sample(model2,
                             data = standat,
                             num_warmup = 500,
                             num_draws = 250,
                             num_chains = 3,
                             cores = 1,
                             refresh = 0,
                             seed = 60471,
                             init = 0.01,
                             pars = c("u", "r_1", "z_1", "L_1", "sd_1"),
                             include = FALSE)
        nutpie_to_brms(irldat_formula, dat, model = stanmodel2,
                       fit = fit, family = empty_cumulative(),
                       stanvars = sv2)
    }, .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

#### convert fits to brms and save results ####
# set filename pars
#stem1 <- paste0("_nutpiefit_mvprob_het", n_subjects_per_group)
saveRDS(brmsnutfits2, paste0(cache_dir, "/", Sys.Date(), stem1, "_mvprobhet_", M,".RDS"))
rm(brmsnutfits2) ; gc()





# Sample with second prior
plan(multisession, workers=10)

brmsnutfits3 <- list(ls_dat,
                     ls_standat3) |>
    future_pmap( \(dat, standat){
        fit <- nutpie_sample(model3,
                             data = standat,
                             num_warmup = 500,
                             num_draws = 250,
                             num_chains = 3,
                             cores = 1,
                             refresh = 0,
                             seed = 60471,
                             init = 0.01,
                             pars = c("u", "r_1", "z_1", "L_1", "sd_1", "w_raw"),
                             include = FALSE)
        nutpie_to_brms(irldat_formula, dat, model = stanmodel3,
                       fit = fit, family = empty_cumulative(),
                       stanvars = sv3)
    }, .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

#### convert fits to brms and save results ####
# set filename pars
#stem1 <- paste0("_nutpiefit_mvprob_het", n_subjects_per_group)
saveRDS(brmsnutfits3, paste0(cache_dir, "/", Sys.Date(), stem1, "_mvstudt_n3_", M,".RDS"))
rm(brmsnutfits3) ; gc()



# Sample with second prior
plan(multisession, workers=10)
brmsnutfits4 <- list(ls_dat,
                     ls_standat4) |>
    future_pmap( \(dat, standat){
        fit <- nutpie_sample(model4,
                             data = standat,
                             num_warmup = 500,
                             num_draws = 250,
                             num_chains = 3,
                             cores = 1,
                             refresh = 0,
                             seed = 60471,
                             init = 0.01,
                             pars = c("u", "r_1", "z_1", "L_1", "sd_1", "w_raw"),
                             include = FALSE)
        nutpie_to_brms(irldat_formula, dat, model = stanmodel4,
                       fit = fit, family = empty_cumulative(),
                       stanvars = sv4)
    }, .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

#### convert fits to brms and save results ####
# set filename pars
saveRDS(brmsnutfits4, paste0(cache_dir, "/", Sys.Date(), stem1, "_mvstudt_nubeta0.1_", M,".RDS"))
rm(brmsnutfits4) ; gc()



# Sample with second prior
plan(multisession, workers=10)
brmsnutfits5 <- list(ls_dat,
                     ls_standat5) |>
    future_pmap( \(dat, standat){
        fit <- nutpie_sample(model5,
                             data = standat,
                             num_warmup = 500,
                             num_draws = 250,
                             num_chains = 3,
                             cores = 1,
                             refresh = 0,
                             seed = 60471,
                             init = 0.01,
                             pars = c("u", "r_1", "z_1", "L_1", "sd_1", "w_raw"),
                             include = FALSE)
        nutpie_to_brms(irldat_formula, dat, model = stanmodel5,
                       fit = fit, family = empty_cumulative(),
                       stanvars = sv5)
    }, .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

#### convert fits to brms and save results ####
# set filename pars
saveRDS(brmsnutfits5, paste0(cache_dir, "/", Sys.Date(), stem1, "_mvstudt_nubeta0.9_", M,".RDS"))
rm(brmsnutfits5) ; gc()



# Sample with second prior
plan(multisession, workers=10)
brmsnutfits6 <- list(ls_dat,
                     ls_standat6) |>
    future_pmap( \(dat, standat){
        fit <- nutpie_sample(model6,
                             data = standat,
                             num_warmup = 500,
                             num_draws = 250,
                             num_chains = 3,
                             cores = 1,
                             refresh = 0,
                             seed = 60471,
                             init = 0.01,
                             pars = c("u", "r_1", "z_1", "L_1", "sd_1", "w_raw"),
                             include = FALSE)
        nutpie_to_brms(irldat_formula, dat, model = stanmodel6,
                       fit = fit, family = empty_cumulative(),
                       stanvars = sv6)
    }, .progress = TRUE, .options = furrr_options(seed = 123))
plan(sequential)

#### convert fits to brms and save results ####
# set filename pars
saveRDS(brmsnutfits6, paste0(cache_dir, "/", Sys.Date(), stem1, "_mvstudthetsked_nu3_", M,".RDS"))
rm(brmsnutfits6) ; gc()







##### Posterior predictions #####

# data for predictions
new_data <- data.frame(time2 = c(0, 12, 0, 12),
                       group = c(0, 0, 1, 1)) |>
    mutate(interaction = group * time2)


# Load fits and predict effect sizes
brmsnutfits1 <- readRDS(paste0(cache_dir, "/", "2026-08-31_compare_copula_Q6_n50_mvprob_400.RDS"))
summaries_mvprobhet1 <- brmsnutfits1 |>
    map(\(fit) {
        epreds <- post_epred_mvprobit_hetsked(fit, new_data, mode = "homoskedastic")
        sum_epreds <- apply(epreds, c(1, 2), sum) -12
        (sum_epreds[,4] - sum_epreds[,2]) - (sum_epreds[,3] - sum_epreds[,1])
    })
rm(brmsnutfits1) ; gc()

brmsnutfits2 <- readRDS(paste0(cache_dir, "/", "2026-08-31_compare_copula_Q6_n50_mvprobhet_400.RDS"))
summaries_mvprobhet2 <- brmsnutfits2 |>
    map(\(fit) {
        epreds <- post_epred_mvprobit_hetsked(fit, new_data, mode = c("heteroskedastic"))
        sum_epreds <- apply(epreds, c(1, 2), sum) -12
        (sum_epreds[,4] - sum_epreds[,2]) - (sum_epreds[,3] - sum_epreds[,1])
    })
rm(brmsnutfits2) ; gc()

brmsnutfits3 <- readRDS(paste0(cache_dir, "/", "2026-08-31_compare_copula_Q6_n50_mvstudt_n3_400.RDS"))
summaries_mvprobhet3 <- brmsnutfits3 |>
    map(\(fit) {
        epreds <- post_epred_studthetsked(fit, new_data, nu =3, mode = c("homoskedastic"))
        sum_epreds <- apply(epreds, c(1, 2), sum) -12
        (sum_epreds[,4] - sum_epreds[,2]) - (sum_epreds[,3] - sum_epreds[,1])
    })
rm(brmsnutfits3) ; gc()

brmsnutfits4 <- readRDS(paste0(cache_dir, "/", "2026-08-31_compare_copula_Q6_n50_mvstudt_nubeta0.1_400.RDS"))
summaries_mvprobhet4 <- brmsnutfits4 |>
    map(\(fit) {
        epreds <- post_epred_studthetsked(fit, new_data, mode = c("homoskedastic"))
        sum_epreds <- apply(epreds, c(1, 2), sum) -12
        (sum_epreds[,4] - sum_epreds[,2]) - (sum_epreds[,3] - sum_epreds[,1])
    })
rm(brmsnutfits4) ; gc()

brmsnutfits5 <- readRDS(paste0(cache_dir, "/", "2026-08-31_compare_copula_Q6_n50_mvstudt_nubeta0.9_400.RDS"))
summaries_mvprobhet5 <- brmsnutfits5 |>
    map(\(fit) {
        epreds <- post_epred_studthetsked(fit, new_data,mode = c("homoskedastic"))
        sum_epreds <- apply(epreds, c(1, 2), sum) -12
        (sum_epreds[,4] - sum_epreds[,2]) - (sum_epreds[,3] - sum_epreds[,1])
    })
rm(brmsnutfits5) ; gc()

brmsnutfits6 <- readRDS(paste0(cache_dir, "/", "2026-08-31_compare_copula_Q6_n50_mvstudthetsked_nu3_400.RDS"))
summaries_mvprobhet6 <- brmsnutfits6 |>
    map(\(fit) {
        epreds <- post_epred_studthetsked(fit, new_data, nu = 3.0, mode = c("het2par"))
        sum_epreds <- apply(epreds, c(1, 2), sum) -12
        (sum_epreds[,4] - sum_epreds[,2]) - (sum_epreds[,3] - sum_epreds[,1])
    })
rm(brmsnutfits6) ; gc()



##### summarise predictions
summ_mvphet1 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "mv_probit",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_mvprobhet1[[.x]],
                                    lab = "homoskedastic"))
summ_mvphet2 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "mv_probit-het",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_mvprobhet2[[.x]],
                                    lab = "heteroskedastic"))
summ_mvphet3 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "mv_studt",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_mvprobhet3[[.x]],
                                    lab = "homoskedastic, nu = 3"))
summ_mvphet4 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "mv_studt",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_mvprobhet4[[.x]],
                                    lab = "homoskedastic, nu-beta = 0.1"))
summ_mvphet5 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "mv_studt",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_mvprobhet5[[.x]],
                                    lab = "homoskedastic, nu-beta = 0.9"))
summ_mvphet6 <- map_df(1:length(eff_probs),
                       ~ data.frame(id = .x, type = "mv_studt_nu 3 het",
                                    eff_prob = eff_probs[[.x]], estimate = summaries_mvprobhet6[[.x]],
                                    lab = "heteroskedastic, nu = 3"))

# combine summaries
all_summs <- rbind(summ_mvphet1, summ_mvphet2, summ_mvphet3,
                   summ_mvphet4, summ_mvphet5, summ_mvphet6)



summs <- all_summs |>
    summarise(.by = c(id, type, lab),
              eff_prob = unique(eff_prob),
              q_025 = quantile(estimate, 0.025),
              mean = mean(estimate),
              q_975 = quantile(estimate, 0.975),
              flag1 = ifelse(eff_prob == 0.5 & q_025 > 0 | q_975 < 0, 1, 0),
              flag2 = ifelse(q_025 > 0, 1, 0)) |>
    summarise(.by = c(type, lab, eff_prob),
              type1error = mean(flag1, na.rm = TRUE),
              power = mean(flag2))
summs |>
    filter(eff_prob == 0.5) |>
    dplyr::select(type, lab, eff_prob, type1error) |>
    arrange(type1error)

summs |>
    filter(eff_prob == 0.75) |>
    dplyr::select(type, lab, eff_prob, power) |>
    arrange(power)

plot_mv_probit <- all_summs %>%
    mutate(eff_prob = factor(eff_prob)) %>%
    ggplot(aes(x =  id, y = estimate)) +
    geom_hline(yintercept = 0) +
    stat_pointinterval(position = position_dodge(width = 0.3)) +
    #    expand_limits(y = shared_range) +
    ggtitle(label = "Difference in expected value by group",
            subtitle = paste0(M, " simulations,", n_subjects_per_group, " subjects per group")) +
    facet_wrap(~ eff_prob + lab, ncol = 6)#, scales = "free_y")

plot_mv_probit





# Posterior nu
nu3 <- brmsnutfits3 |>
    map(\(fit)
        as_draws_df(fit)$nu )
nu4 <- brmsnutfits4 |>
    map(\(fit)
        as_draws_df(fit)$nu )
nu5 <- brmsnutfits5 |>
    map(\(fit)
        as_draws_df(fit)$nu )


# posterior het slops
slopes2 <- brmsnutfits2 |>
    imap(\(fit, idx)
         as_draws_df(fit) |> as_tibble() |>
             select(starts_with( c("scale_intercept", "scale_slope", "nu") )) |>
             mutate(m = idx,
                    model = "mvp-het")) |>
    bind_rows()
slopes6 <- brmsnutfits6 |>
    imap(\(fit, idx)
         as_draws_df(fit) |> as_tibble() |>
             select(starts_with( c("scale_intercept", "scale_slope", "nu") )) |>
             mutate(m = idx,
                    model = "mvstudt-het")) |>
    bind_rows()

all_slopes <- bind_rows(slopes2,  slopes6) |>
    janitor::clean_names()

all_slopes |> summary()

all_slopes |>
    pivot_longer(cols = c(-model, -m),
                 names_prefix = "scale_slope_",
                 names_to = "slope") |>
    ggplot(aes(x = value, col = as.factor(model))) + geom_density() +
    facet_wrap(~slope, scales = "free")


all_slopes |>
    filter(m == 10) |>
    ggplot(aes(x = log(nu), y = scale_slope_9)) + geom_point() + facet_wrap(~model)


