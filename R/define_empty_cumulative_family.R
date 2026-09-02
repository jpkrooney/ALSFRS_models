empty_cumulative <- function() {
  custom_family(
    "empty_cumulative",
    dpars = c("mu"),
    links = c("identity"),
    lb = 1,
    type = "int",
    threshold = "flexible",
    specials = c("ordinal", "ordered_thres", "thres_minus_eta"),
    posterior_epred = posterior_epred_empty_cumulative
  )
}

posterior_epred_empty_cumulative <- function(prep) {
  prep$family$family <- "cumulative"
  prep$family$link <- "probit"
  brms:::posterior_epred_ordinal(prep)
}

posterior_predict_mv_probit <- function(
  fit,
  nsamples = NULL,
  subset = NULL,
  ...
) {
  subset <- brms:::validate_draw_ids(fit, subset, nsamples)
  linpred <- posterior_linpred(fit, transform = FALSE, subset = subset, ...)

  prep <- prepare_predictions(fit, subset = subset, ...)
  N_dims <- length(prep$resps)

  rescor <- get_rescor(fit, size = N_dims, subset = subset)

  out <- array(NA_integer_, c(prep$ndraws, prep$nobs, N_dims))
  for (s in seq_len(prep$ndraws)) {
    for (o in seq_len(prep$nobs)) {
      mu_noise <- brms::rmulti_normal(
        1,
        mu = linpred[s, o, ],
        Sigma = rescor[s, , ]
      )
      for (resp_id in 1:N_dims) {
        thres <- prep$resps[[resp_id]]$thres
        out[s, o, resp_id] <- sum(thres$thres[s, ] < mu_noise[resp_id])
      }
    }
  }
  out + 1
}

get_rescor <- function(fit, size, subset = NULL) {
  rescor <- as.matrix(
    fit,
    variable = "^rescor\\[",
    regex = TRUE,
    subset = subset
  )

  nsamples <- dim(rescor)[1]
  out <- array(NA_real_, dim = c(nsamples, size, size))
  for (i in seq_len(size)) {
    out[, i, i] <- 1
  }
  stopifnot(min(rescor) >= -1, max(rescor) <= 1)
  stopifnot(ncol(rescor) == size * (size - 1) / 2)
  k <- 0
  for (i in seq_len(size)[-1]) {
    for (j in seq_len(i - 1)) {
      k <- k + 1
      out[, j, i] <- out[, i, j] <- rescor[, paste0("rescor[", k, "]")]
    }
  }
  stopifnot(all(!is.na(out)))
  out
}


rmulti_normal_custom <- function(mu, Sigma, check = FALSE) {
  n <- dim(mu)[1]
  p <- dim(mu)[2]
  if (check) {
    if (!(is_wholenumber(n) && n > 0)) {
      stop2("n must be a positive integer.")
    }
    if (!all(dim(Sigma) == c(p, p))) {
      stop2("Dimension of Sigma is incorrect.")
    }
    if (!is_symmetric(Sigma)) {
      stop2("Sigma must be a symmetric matrix.")
    }
  }
  samples <- matrix(rnorm(n * p), nrow = n, ncol = p, byrow = TRUE)
  mu + samples %*% chol(Sigma)
}


posterior_predict_mv_probit2 <- function(
  fit,
  nsamples = NULL,
  subset = NULL,
  ...
) {
  subset <- brms:::validate_draw_ids(fit, subset, nsamples)
  linpred <- posterior_linpred(fit, transform = FALSE, subset = subset, ...)
  prep <- prepare_predictions(fit, subset = subset, ...)
  N_dims <- length(prep$resps)
  rescor <- get_rescor(fit, size = N_dims, subset = subset)

  out <- array(NA_integer_, c(prep$ndraws, prep$nobs, N_dims))
  for (s in seq_len(prep$ndraws)) {
    mu_noise <- rmulti_normal_custom(mu = linpred[s, , ], Sigma = rescor[s, , ])
    for (resp_id in 1:N_dims) {
      thres <- prep$resps[[resp_id]]$thres
      temp <- rep(thres$thres[s, ], each = prep$nobs) < mu_noise[, resp_id]
      dim(temp) <- c(prep$nobs, length(thres$thres[s, ]))
      out[s, , resp_id] <- rowSums(temp)
    }
  }
  out + 1
}
