#FUNCTION FOR THE LMER SIMPLE EFFECTS
effects_lmer <- function(fit) {
    # Estimating conf ints. for a linear combination
    # Based on https://bookdown.org/ccolonescu/RPoE4/intervalest.html
    vc <- vcov(fit)

    lincomb <- rep(0, nrow(vc))
    names(lincomb) <- rownames(vc)
    lincomb["groupTreatment"] <- 0 #changed, it was 1 previously
    lincomb["alsfrs_dly_mnths:groupTreatment"] <- 12

    alpha <- c(0.5, 0.05)
    estimate <- as.numeric(lincomb %*% fit@beta)
    df <- df.residual(fit)
    tcr <- qt(1-alpha/2, df)
    estimate_se <- as.numeric(sqrt( lincomb %*% (vc %*% lincomb)))
    lower <- estimate - tcr*estimate_se
    upper <- estimate + tcr*estimate_se
    data.frame(estimate = estimate, lower95 = lower[2], lower50 = lower[1], upper50 = upper[1],
               upper95 = upper[2])
}




#FUNCTION FOR THE LMER SPLINES EFFECTS
effects_lmer_splines <- function(fit) {
    # Estimating conf ints. for a linear combination
    # Based on https://bookdown.org/ccolonescu/RPoE4/intervalest.html
    vc <- vcov(fit)

    spline_coeffs_0 <- get_spline_coeffs(0)
    spline_coeffs_12 <- get_spline_coeffs(12)
    spline_coeffs_diff <- spline_coeffs_12 - spline_coeffs_0

    lincomb <- rep(0, nrow(vc))
    names(lincomb) <- rownames(vc)
    lincomb["groupTreatment"] <- 0 #changed, it was 1 previously
    lincomb["dly_mnths_spline_1:groupTreatment"] <- spline_coeffs_diff[1]
    lincomb["groupTreatment:dly_mnths_spline_2"] <- spline_coeffs_diff[2]

    alpha <- c(0.5, 0.05)
    estimate <- as.numeric(lincomb %*% fit@beta)
    df <- df.residual(fit)
    tcr <- qt(1-alpha/2, df)
    estimate_se <- as.numeric(sqrt( lincomb %*% (vc %*% lincomb)))
    lower <- estimate - tcr*estimate_se
    upper <- estimate + tcr*estimate_se
    data.frame(estimate = estimate, lower95 = lower[2], lower50 = lower[1], upper50 = upper[1],
               upper95 = upper[2])
}


#FUNCTION FOR THE LOGIT EFFECTS
#(the only change in the original function is the -1 in the col_values)
effects_logit <- function(fit) {
    predict_data <- data.frame(group = c("Control", "Control", "Treatment", "Treatment"),
                               alsfrs_dly_mnths = c(0, 12, 0, 12))

    pred <- posterior_epred(fit, newdata = predict_data, re_formula = NA)

    expected_sum <- function(x) {
        col_values <- as.integer(colnames(x)) - 1 #back to the real values of the categories
        multiplied <- sweep(x, MARGIN = 2, STATS = col_values, FUN = "*")
        rowSums(multiplied)
    }

    start_control <- expected_sum(pred[,1,])
    end_control <- expected_sum(pred[,2,])
    start_treatment <- expected_sum(pred[,3,])
    end_treatment <- expected_sum(pred[,4,])

    (end_treatment - start_treatment) - (end_control - start_control)
}

effect_from_sim_study <- function(fit) {
    predict_data <- data.frame(group = c("Control", "Control", "Treatment", "Treatment"),
                               alsfrs_dly_mnths = c(0, 12, 0, 12))

    pred <- brms::posterior_epred(fit, newdata = predict_data, re_formula = NA)

    expected_sum <- function(x) {
        col_values <- as.integer(colnames(x))
        multiplied <- sweep(x, MARGIN = 2, STATS = col_values, FUN = "*")
        rowSums(multiplied)
    }

    start_control <- expected_sum(pred[,1,])
    end_control <- expected_sum(pred[,2,])
    start_treatment <- expected_sum(pred[,3,])
    end_treatment <- expected_sum(pred[,4,])

    n_divergent <- sum(rstan::get_divergent_iterations(fit$fit))
    ss <- posterior::summarise_draws(posterior::as_draws(fit))
    max_rhat <- max(ss$rhat, na.rm = TRUE)
    min_ess_bulk <- min(ss$ess_bulk, na.rm = TRUE)
    min_ess_tail <- min(ss$ess_tail, na.rm = TRUE)

    list(
        effect_samples = (end_treatment - start_treatment) - (end_control - start_control),
        diags = data.frame(n_divergent, max_rhat, min_ess_bulk, min_ess_tail)
    )

}


#COLLECT THE PER-FIT DIAGNOSTICS FROM ONE summary_* LIST INTO A TIDY DATA FRAME
# summary_list: output of brm_parallel(summarise_fun = effect_from_sim_study)
# type:         label for the model, as used in the effects_* data frames
# eff_prob:     vector of effect probabilities, one per fit (defaults to the
#               `eff_probs` object the simulation script keeps in scope)
# Fits that failed (NULL, or no $diags) yield a row of NAs so that the rows stay
# aligned with eff_prob.
collect_diags <- function(summary_list, type, eff_prob = eff_probs) {
    stopifnot(length(summary_list) == length(eff_prob))

    empty <- data.frame(n_divergent = NA_real_, max_rhat = NA_real_,
                        min_ess_bulk = NA_real_, min_ess_tail = NA_real_)

    d <- purrr::map_df(summary_list, function(s) {
        if (is.null(s) || is.null(s$diags)) empty else s$diags
    })

    cbind(data.frame(id = seq_along(summary_list), type = type, eff_prob = eff_prob), d)
}


#COLLECT DIAGNOSTICS FROM SEVERAL summary_* LISTS AT ONCE
# summaries: a *named* list of summary_* objects; names become the `type` column
collect_diags_all <- function(summaries, eff_prob = eff_probs) {
    purrr::imap_dfr(summaries, function(s, nm) collect_diags(s, nm, eff_prob))
}


#SUMMARISE THE PER-FIT DIAGNOSTICS BY MODEL AND EFFECT PROBABILITY
# A mean is the wrong summary for most of these: what matters is whether *any*
# fit misbehaved, so the tail statistics (worst rhat, smallest ESS, share of
# fits with divergences) are reported alongside the means.
summarise_diags <- function(diags, ess_threshold = 400, rhat_threshold = 1.01) {
    diags %>%
        dplyr::group_by(type, eff_prob) %>%
        dplyr::summarise(
            n_fits        = dplyr::n(),
            n_failed      = sum(is.na(max_rhat)),
            pct_divergent = mean(n_divergent > 0, na.rm = TRUE),
            mean_divergent= mean(n_divergent, na.rm = TRUE),
            max_divergent = max(n_divergent, na.rm = TRUE),
            worst_rhat    = max(max_rhat, na.rm = TRUE),
            pct_bad_rhat  = mean(max_rhat > rhat_threshold, na.rm = TRUE),
            worst_ess_bulk= min(min_ess_bulk, na.rm = TRUE),
            median_ess_bulk = median(min_ess_bulk, na.rm = TRUE),
            pct_low_ess_bulk = mean(min_ess_bulk < ess_threshold, na.rm = TRUE),
            worst_ess_tail= min(min_ess_tail, na.rm = TRUE),
            pct_low_ess_tail = mean(min_ess_tail < ess_threshold, na.rm = TRUE),
            .groups = "drop"
        )
}
