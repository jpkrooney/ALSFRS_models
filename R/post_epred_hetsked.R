post_epred_hetsked <- function(fit, new_data, beta_names = c("_time2", "_groupTreatment", "_time2:groupTreatment")){

    # Extract posterior samples
    post <- as.data.frame(fit)

    # Extract required metadata
    n_samps <- nrow(post)
    n_obs <- nrow(new_data)
    d <- length(fit$formula$responses)  # num response vars
    respnames <- fit$formula$responses  # names of responses
    thresnames <- fit$fit@par_dims[ paste0("Intercept_", respnames)]
    nthres_idx <- unlist(thresnames)
    b_threscolnames <- lapply(1:length(respnames), function(i)
        paste0(paste0("b_", respnames[i], "_Intercept"), "[", 1:nthres_idx[i], "]"))

    # Extract estimated thresholds/cutpoints for ordinal outcomes
    thresholds_list <- lapply(b_threscolnames, function(i)
        data.frame(-Inf, post[, i], Inf)
        ) # is +Inf needed?

    # Extract beta coefficients for each outcome
    betas_list <- lapply(1:d, function(d)
        post[, paste0("b_Q", sprintf(fmt = "%02d", d), beta_names)])

    # Extract sigma coefficients for each outcome
    sig_intercepts_list  <- lapply(1:d, function(d)
        post[, paste0("scale_intercept[", d, "]")])
    sig_slopes_list <- lapply(1:d, function(d)
        post[, paste0("scale_slope[", d, "]")])


    # Result array: [samples, observations, dimensions]
    epred_results <- array(0, dim = c(n_samps, n_obs, d))

    for (d in 1:12) {
        # Extract parameters for this dimension
        sig_int  <- post[, paste0("scale_intercept[", d, "]")]
        sig_slope <- post[, paste0("scale_slope[", d, "]")]

        # Current thresholds for this dimension
        thr_mat <- thresholds_list[[d]]

        mu_val_all <- as.matrix(betas_list[[d]]) %*% t(new_data[, ])
        sigma_val_all <- exp(sig_int + outer(sig_slope, new_data[, 1], "*"))

        # copy matrices into 3rd dimension for calculations
        thr_mat3d <- array( as.matrix(thr_mat),
                           dim = c(n_samps, ncol(thr_mat), n_obs) )
        mu_val3d <- aperm( array(as.matrix(mu_val_all),
                                 dim = c(n_samps, n_obs, ncol(thr_mat))),
                           c(1, 3, 2) )
        sig_val3d <- aperm( array(as.matrix(sigma_val_all),
                                 dim = c(n_samps, n_obs, ncol(thr_mat))),
                           c(1, 3, 2) )


        # calculate probabilities for ordinal thresholds
        eta3d <- (thr_mat3d - mu_val3d) / sig_val3d
        prob_cum3d <- pnorm(eta3d)
        #prob_k3d <- aperm(apply(prob_cum3d, c(1, 3), diff), c(2, 1, 3))
        prob_k3d <- apply(prob_cum3d, c(1, 3), diff)


        # calculate expected values
        pdims <- dim(prob_k3d ) # store 3D dims
        dim(prob_k3d ) <- c(pdims[1], pdims[2]*pdims[3]) #flatten to 2D
        expected_y <- t(prob_k3d) %*% (1:(nthres_idx[d] + 1))
        dim(expected_y) <- c(n_samps, pdims[3])


        # pack into result array
        epred_results[, , d] <- expected_y
    }
    return(epred_results)
}




post_epred_studthetsked <- function(fit, new_data, beta_names = c("_time2", "_groupTreatment", "_time2:groupTreatment"),
                                    nu = NULL, mode = c("homoskedastic", "het1par", "het2par")){

    # Extract posterior samples
    post <- as.data.frame(fit)

    # if nu = NULL, does a nu paramater exist in post?
    if(is.null(nu) & !any(names(post) == "nu")) {
        stop ("If a value for nu is not provided if should be a parameter of the model")
    }

    # Extract required metadata
    n_samps <- nrow(post)
    n_obs <- nrow(new_data)
    D <- length(fit$formula$responses)  # num response vars
    respnames <- fit$formula$responses  # names of responses
    thresnames <- fit$fit@par_dims[ paste0("Intercept_", respnames)]
    nthres_idx <- unlist(thresnames)
    b_threscolnames <- lapply(1:length(respnames), function(i)
        paste0(paste0("b_", respnames[i], "_Intercept"), "[", 1:nthres_idx[i], "]"))

    # Extract estimated thresholds/cutpoints for ordinal outcomes
    thresholds_list <- lapply(b_threscolnames, function(i)
        data.frame(-Inf, post[, i], Inf)
    ) # is +Inf needed?

    # Extract beta coefficients for each outcome
    betas_list <- lapply(1:D, function(d)
        post[, paste0("b_Q", sprintf(fmt = "%02d", d), beta_names)])

    # Extract sigma coefficients for each outcome
    if(mode != "homoskedastic"){
        if(mode == "het2par"){
            sig_intercepts_list  <- lapply(1:D, function(d)
                post[, paste0("scale_intercept[", d, "]")])
        }
        sig_slopes_list <- lapply(1:D, function(d)
            post[, paste0("scale_slope[", d, "]")])
    }

    # Extract L_rescor
    L_rescor <- post[, grep("L_rescor", names(post), value = TRUE)]

    # Extract student-t df related params
    if(is.null(nu)){
        nu <- post$nu
    } else {
        nu <- rep(nu, n_samps)
    }

    # Result array: [samples, observations, dimensions]
    epred_results <- array(0, dim = c(n_samps, n_obs, D))

    for (d in 1:D) {
        # Extract parameters for this dimension
        #sig_int  <- post[, paste0("scale_intercept[", d, "]")]
        #sig_slope <- post[, paste0("scale_slope[", d, "]")]

        # Current thresholds for this dimension
        thr_mat <- thresholds_list[[d]]

        mu_val_all <- as.matrix(betas_list[[d]]) %*% t(new_data[, ])

        # prepare variances
        if(mode == "homoskedastic"){
            sigma_val_all <- rep(1, n_samps) # marginal
        } else if(mode == "het1par"){
            log_het_scale <- outer(sig_slopes_list[[d]], new_data[, 1], "*")
            sigma_val_all <- exp(log_het_scale) # rep(1, n_samps)) - should add but since log(1) = 0 can skip
        } else if(mode == "het2par"){
            log_het_scale <- outer(sig_slopes_list[[d]], new_data[, 1], "*")
            sigma_val_all <- exp(sig_intercepts_list[[d]] + log_het_scale) # rep(1, n_samps)) - should add but since log(1) = 0 can skip
        }

        # copy matrices into 3rd dimension for calculations
        thr_mat3d <- array( as.matrix(thr_mat),
                            dim = c(n_samps, ncol(thr_mat), n_obs) )
        mu_val3d <- aperm( array(as.matrix(mu_val_all),
                                 dim = c(n_samps, n_obs, ncol(thr_mat))),
                           c(1, 3, 2) )
        sig_val3d <- aperm( array(as.matrix(sigma_val_all),
                                  dim = c(n_samps, n_obs, ncol(thr_mat))),
                            c(1, 3, 2) )


        # calculate probabilities for ordinal thresholds
        eta3d <- (thr_mat3d - mu_val3d) / sig_val3d
        #prob_cum3d <- pnorm(eta3d)
        prob_cum3d <- pt(eta3d, nu)
        #prob_k3d <- aperm(apply(prob_cum3d, c(1, 3), diff), c(2, 1, 3))
        prob_k3d <- apply(prob_cum3d, c(1, 3), diff)


        # calculate expected values
        pdims <- dim(prob_k3d ) # store 3D dims
        dim(prob_k3d ) <- c(pdims[1], pdims[2]*pdims[3]) #flatten to 2D
        expected_y <- t(prob_k3d) %*% (1:(nthres_idx[d] + 1))
        dim(expected_y) <- c(n_samps, pdims[3])


        # pack into result array
        epred_results[, , d] <- expected_y
    }
    return(epred_results)
}






post_epred_mvprobit_hetsked <- function(fit, new_data, beta_names = c("_time2", "_groupTreatment", "_time2:groupTreatment"),
                                    mode = c("homoskedastic", "heteroskedastic")){

    # Extract posterior samples
    post <- as.data.frame(fit)

    # Extract required metadata
    n_samps <- nrow(post)
    n_obs <- nrow(new_data)
    D <- length(fit$formula$responses)  # num response vars
    respnames <- fit$formula$responses  # names of responses
    thresnames <- fit$fit@par_dims[ paste0("Intercept_", respnames)]
    nthres_idx <- unlist(thresnames)
    b_threscolnames <- lapply(1:length(respnames), function(i)
        paste0(paste0("b_", respnames[i], "_Intercept"), "[", 1:nthres_idx[i], "]"))

    # Extract estimated thresholds/cutpoints for ordinal outcomes
    thresholds_list <- lapply(b_threscolnames, function(i)
        data.frame(-Inf, post[, i], Inf)
    ) # is +Inf needed?

    # Extract beta coefficients for each outcome
    betas_list <- lapply(1:D, function(d)
        post[, paste0("b_Q", sprintf(fmt = "%02d", d), beta_names)])

    # Extract sigma coefficients for each outcome
    if(mode != "homoskedastic"){
        sig_intercepts_list  <- lapply(1:D, function(d)
            post[, paste0("scale_intercept[", d, "]")])
        sig_slopes_list <- lapply(1:D, function(d)
            post[, paste0("scale_slope[", d, "]")])
    }

    # Extract L_rescor
    L_rescor <- post[, grep("L_rescor", names(post), value = TRUE)]


    # Result array: [samples, observations, dimensions]
    epred_results <- array(0, dim = c(n_samps, n_obs, D))

    for (d in 1:D) {
        # Extract parameters for this dimension
        #sig_int  <- post[, paste0("scale_intercept[", d, "]")]
        #sig_slope <- post[, paste0("scale_slope[", d, "]")]

        # Current thresholds for this dimension
        thr_mat <- thresholds_list[[d]]

        mu_val_all <- as.matrix(betas_list[[d]]) %*% t(new_data[, ])

        # prepare variances
        if(mode == "homoskedastic"){
            sigma_val_all <- rep(1, n_samps) # marginal
        } else if(mode == "heteroskedastic"){
            log_het_scale <- outer(sig_slopes_list[[d]], new_data[, 1], "*")
            sigma_val_all <- exp(sig_intercepts_list[[d]] + log_het_scale) # rep(1, n_samps)) - should add but since log(1) = 0 can skip
        }

        # copy matrices into 3rd dimension for calculations
        thr_mat3d <- array( as.matrix(thr_mat),
                            dim = c(n_samps, ncol(thr_mat), n_obs) )
        mu_val3d <- aperm( array(as.matrix(mu_val_all),
                                 dim = c(n_samps, n_obs, ncol(thr_mat))),
                           c(1, 3, 2) )
        sig_val3d <- aperm( array(as.matrix(sigma_val_all),
                                  dim = c(n_samps, n_obs, ncol(thr_mat))),
                            c(1, 3, 2) )


        # calculate probabilities for ordinal thresholds
        eta3d <- (thr_mat3d - mu_val3d) / sig_val3d
        prob_cum3d <- pnorm(eta3d)
        #prob_k3d <- aperm(apply(prob_cum3d, c(1, 3), diff), c(2, 1, 3))
        #prob_k3d <- apply(prob_cum3d, c(1, 3), diff)
        prob_k3d <- apply(prob_cum3d, c(1, 3), diff)

        # calculate expected values
        pdims <- dim(prob_k3d ) # store 3D dims
        dim(prob_k3d ) <- c(pdims[1], pdims[2]*pdims[3]) #flatten to 2D
        expected_y <- t(prob_k3d) %*% (1:(nthres_idx[d] + 1))
        dim(expected_y) <- c(n_samps, pdims[3])


        # pack into result array
        epred_results[, , d] <- expected_y
    }
    return(epred_results)
}









