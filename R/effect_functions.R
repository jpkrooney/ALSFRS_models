#FUNCTION FOR THE LMER SIMPLE EFFECTS
effects_lmer <- function(fit) {
    # Estimating conf ints. for a linear combination
    # Based on https://bookdown.org/ccolonescu/RPoE4/intervalest.html
    vc <- vcov(fit)

    lincomb <- rep(0, nrow(vc))
    names(lincomb) <- rownames(vc)
    lincomb["groupTreatment"] <- 1
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
    lincomb["groupTreatment"] <- 1
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
