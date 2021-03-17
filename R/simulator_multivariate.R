simulate_multivariate_probit <- function(N_obs, N_cat, N_dim) {
   for(step in 1:50) {
       # Rejection sampling loop to achieve coverage for all responses
       corr <- rlkj(N_dim, 1, cholesky = FALSE)
       noise <- mvtnorm::rmvnorm(n = N_obs, sigma = corr)

       beta <- rnorm(N_dim, sd = 1)

       X <- rnorm(N_obs, sd = 5)

       thresholds <- matrix(nrow = N_cat - 1, ncol = N_dim)
       for(dim in 1:N_dim) {
           thresholds[,dim] <- sort(rnorm(N_cat - 1, sd = 3))
       }

       Y <- matrix(nrow = N_obs, ncol = N_dim)
       for(n in 1:N_obs) {
        Y_latent <- X[n] * beta + noise[n,]
        Y[n, ] <- colSums(sweep(thresholds, MARGIN = 2, STATS = Y_latent, FUN = "<")) + 1
       }

       col_mins <- apply(Y, MARGIN = 2, FUN = min)
       col_maxs <- apply(Y, MARGIN = 2, FUN = max)
       if(all(col_mins == 1) && all(col_maxs == N_cat)) {
           break;
       }
   }

   if(step > 1) {
       cat("Required", step, "rejection steps.\n")
   }

   colnames(Y) <- paste0("q", 1:N_dim)

   observed_df <- cbind(as.data.frame(Y),
                        data.frame(X = X))

   corr_upper <- numeric(choose(N_dim, 2))
   for (k in 1:N_dim) {
       for (j in 1:(k - 1)) {
           corr_upper[choose(k - 1, 2) + j] = corr[j, k];
       }
   }

   list(
       true = loo::nlist(
           thresholds,
           beta,
           corr,
           corr_upper
       ), observed_df =
           observed_df
   )
}


generator_multivariate_probit_approx <- function(N_obs, N_cat, N_dim, disc = 100) {
    raw <- simulate_multivariate_probit(N_obs, N_cat, N_dim)
    raw$observed_df$obs_id <- 1:nrow(raw$observed_df)
    observed_long <- pivot_longer(raw$observed_df,
                                  starts_with("q"),
                                  names_to = "question",
                                  values_to = "answer")

    f <- brmsformula(
        answer | thres(gr = question) ~ X:question + (0 + question | obs_id),
        disc ~ 1,
        family = cumulative(link = "logit",  link_disc = "identity"))

    # By default, the SD would be 1, but it is instead 1 + the effect of the
    # varying intercept
    # threshold_sd_multiplier <- sqrt(sd_obs_id ^ 2 + 1)

    priors <- c(
        # prior_string(paste0("constant(",sd_obs_id,")"), class = "sd", group = "obs_id"),
        # prior_string(paste0("normal(0, ", threshold_sd_multiplier * 3,")"), class = "Intercept"),
        prior(constant(1), class = "sd", group = "obs_id"),
        prior(normal(0, 3), class = "Intercept"),
        prior(normal(0, 1), class = "b"),
        prior(lkj(1), class = "cor"),
        prior(constant(my_disc), class = "Intercept", dpar = "disc")
    )

    stanvars <- stanvar(name = "my_disc", x= disc)

    stancode <- make_stancode(f,
        prior = priors,
        stanvars = stanvars,
        data = observed_long)

    standata <- make_standata(f,
                              prior = priors,
                              stanvars = stanvars,
                              data = observed_long)
    class(standata) <- NULL

    true = list(
        b = raw$true$beta,
        cor_1 = raw$true$corr_upper
    )

    for(i in 1:N_dim) {
        true[[paste0("Intercept_", i)]] <- raw$true$thresholds[,i]# * threshold_sd_multiplier
    }

    list(
        true = true,
        stancode = stancode,
        observed = standata
    )
}


# Generate a random sample from the LKJ distribution
#
# Function taken from  https://github.com/biobakery/banocc/blob/master/R/rlkj.R
# MIT License.
# This function is based on code from
# \url{https://groups.google.com/forum/#!msg/stan-users/3gDvAs_qwN8/Xpgi2rPlx68J}.
#
# @param d The dimension of the correlation matrix
# @param eta The scaling parameter of the LKJ distribution; must be > 1
#   (eta=1 means the distribution is uniform over d by d correlation matrices)
# @param cholesky Boolean: return the cholesky decomposition?

rlkj <- function(d, eta = 1, cholesky = FALSE) {
    if (d < 2){
        stop("Dimension of correlation matrix must be >= 2")
    }
    if (eta < 1){
        stop("The value of eta must be >= 1")
    }
    alpha <- eta + (d - 2) / 2
    L <- matrix(0, d, d)
    L[1,1] <- 1
    L[-1,1] <- partials <- rgbeta(d - 1, alpha)
    if(d == 2) {
        L[2,2] <- sqrt(1 - L[2,1]^2)
        if(cholesky) return(L)
        Sigma <- tcrossprod(L)
        return(Sigma)
    }
    W <- log(1 - partials^2)
    for(i in 2:(d - 1)) {
        gap <- (i+1):d
        gap1 <- i:(d-1)
        alpha <- alpha - 0.5
        partials <- rgbeta(d - i, alpha)
        L[i,i] <- exp(0.5 * W[i-1])
        L[gap,i] <- partials * exp(0.5 * W[gap1])
        W[gap1] <- W[gap1] + log(1 - partials^2)
    }
    L[d,d] <- exp(0.5 * W[d-1])
    if(cholesky) return(L)
    Sigma <- tcrossprod(L)
    return(Sigma)
}

rgbeta <- function(d, shape) {
    if(shape == Inf)     rep(0, d)
    else if(shape > 0)  -1 + 2 * rbeta(d, shape, shape)
    else if(shape == 0) -1 + 2 * rbinom(d, 1, 0.5)
    else stop("shape must be non-negative")
}
