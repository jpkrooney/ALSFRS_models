empty_cumulative <- function() {
  custom_family(
    "empty_cumulative", dpars = c("mu"),
    links = c("identity"), lb = 1,
    type = "int", threshold = "flexible", specials = c("ordinal", "ordered_thres", "thres_minus_eta"))
}

make_stanvars_mv_probit_base <- function(column_names, rescor_prior_eta = 1) {
  N_dims <- length(column_names)


  stan_funs <- "
      real empty_cumulative_lpmf(int y, real mu, vector intercept) {
        return 0;
      }

    "

  stan_tdata <- paste0("
           int N_thres = nthres_", column_names[1], ";",

    paste0("
         if(N != N_", column_names, ") { reject(\"Requiring equal sample size in all dimensions.\"); }" ,
           collapse = "\n")
  )

  stan_params <- paste0("
      cholesky_factor_corr[", N_dims, "] L_rescor;
    ")


  stan_priors <- paste0("
      target += lkj_corr_cholesky_lpdf(L_rescor | ", rescor_prior_eta, ");
    ")

  stan_thresholds <- paste0("
         vector[N_thres] thresholds[", N_dims, "];",
                            paste0("
         thresholds[", 1:N_dims, "] = Intercept_", column_names, ";", collapse = "")
    )

  stan_genquants <- paste0("
     corr_matrix[", N_dims, "] Rescor = multiply_lower_tri_self_transpose(L_rescor);
     vector<lower=-1,upper=1>[", choose(N_dims, 2) ,"] rescor;
     // extract upper diagonal of rescor matrix
     for (k in 1:", N_dims, ") {
        for (j in 1:(k - 1)) {
          rescor[choose(k - 1, 2) + j] = Rescor[j, k];
        }
      }
    ")

  stanvars_mult_probit <- stanvar(scode = stan_funs, block = "functions") +
    stanvar(scode = stan_tdata, block = "tdata") +
    stanvar(scode = stan_params, block = "parameters") +
    stanvar(scode = stan_thresholds, block = "model", position = "start") +
    stanvar(scode = stan_priors, block = "model", position = "start") +
    stanvar(scode = stan_genquants, block = "genquant")


  stanvars_mult_probit
}

make_stanvars_mv_probit_bgoodri <- function(column_names, rescor_prior_eta = 1) {

  N_dims <- length(column_names)


  stanvars_base <- make_stanvars_mv_probit_base(column_names, rescor_prior_eta = rescor_prior_eta)

  stan_funs <- "
      real approx_Phi(real x) {
        return inv_logit(x * 1.702);
      }

      real approx_inv_Phi(real x) {
        return logit(x) / 1.702;
      }
    "

  stan_tdata <- paste0("
         if(nthres_", column_names[1]," != nthres_", column_names[2:length(column_names)], ") { reject(\"Requiring equal number of thresholds in all dimensions.\"); }" ,
                              collapse = "\n")


  stan_params <- paste0("
      real<lower=0, upper=1> u[N, ", N_dims, "]; // raw residuals
    ")

  stan_likelihood <- paste0("
       for(n in 1:N) {
            real mus[", N_dims, "] = {", paste0("mu_", column_names, "[n]", collapse = ", "), "};
            int Ys[", N_dims, "] = {", paste0("Y_", column_names, "[n]", collapse = ", "), "};

            vector[", N_dims, "] z;
            real prev;
            prev = 0;
            for (d in 1:", N_dims, ") {
              real t; // threshold at which utility = 0
              if (Ys[d] == 1){
                real ub = approx_Phi((thresholds[d, 1] -(mus[d] + prev)) / L_rescor[d,d]);
                t = ub * u[n,d];
                target += log(ub);  // Jacobian adjustment
              } else if (Ys[d] == N_thres + 1) {
                real lb = approx_Phi((thresholds[d, N_thres] -(mus[d] + prev)) / L_rescor[d,d]);
                t = lb + (1 - lb) * u[n,d];
                target += log1m(lb);  // Jacobian adjustment
              } else {
                real lb = approx_Phi((thresholds[d, Ys[d] - 1] -(mus[d] + prev)) / L_rescor[d,d]);
                real ub = approx_Phi((thresholds[d, Ys[d]    ] -(mus[d] + prev)) / L_rescor[d,d]);
                t = lb + (ub - lb) * u[n,d];
                target += log(ub - lb);
              }
              z[d] = approx_inv_Phi(t);
              if (d < ", N_dims, ") prev = L_rescor[d+1,1:d] * head(z, d);
              // Jacobian adjustments imply z is truncated standard normal
              // thus utility --- mu + L_rescor * z --- is truncated multivariate normal
            }
        }
    ")


  stanvars_mult_probit <- stanvars_base +
    stanvar(scode = stan_funs, block = "functions") +
    stanvar(scode = stan_tdata, block = "tdata") +
    stanvar(scode = stan_params, block = "parameters") +
    stanvar(scode = stan_likelihood, block = "likelihood", position = "end")


  stanvars_mult_probit
}


make_stanvars_mv_probit_augmented <- function(column_names, rescor_prior_eta = 1) {

  stanvars_base <- make_stanvars_mv_probit_base(column_names, rescor_prior_eta = rescor_prior_eta)

  N_dims <- length(column_names)


  stan_funs <- "
      real constrain_residual_lp(real z, real mu, int observed, vector thresholds) {
        int N_cat = num_elements(thresholds) + 1;
        if(observed == 1) {
          real ub = thresholds[1]  - mu;
          target += z;
          return(ub - exp(z));
        } else if(observed == N_cat) {
          real lb = thresholds[N_cat - 1] - mu;
          target += z;
          return(lb + exp(z));
        } else {
          real lb = thresholds[observed - 1] - mu;
          real diff = thresholds[observed]  - thresholds[observed - 1];
          real inv_logit_z = inv_logit(z);
          //target += log(diff) + log(inv_logit_z) + log1m(inv_logit_z);
          target += log(diff) - fabs(z) + 2.0 * log1p_exp(-fabs(z));
          return(lb + diff * inv_logit_z);
        }
      }
    "

  stan_params <- paste0("
      real z_rescor[N, ", N_dims, "]; // raw residuals
    ")


  constrain_residuals_code <- paste0(
    "residuals[n, ", 1:N_dims, "] = constrain_residual_lp(z_rescor[n, ", 1:N_dims, "], mu_", column_names, "[n], Y_", column_names, "[n], Intercept_", column_names,");",    collapse = "\n          ")

  stan_likelihood <- paste0("
        {
          vector[", N_dims, "] residuals[N];
          for(n in 1:N) {
              ", constrain_residuals_code , "
          }
          target += multi_normal_cholesky_lpdf(residuals | rep_vector(0, ", N_dims, "), L_rescor);
        }
    ")

  stanvars_mult_probit <-
    stanvars_base +
    stanvar(scode = stan_funs, block = "functions") +
    stanvar(scode = stan_params, block = "parameters") +
    stanvar(scode = stan_likelihood, block = "likelihood", position = "end")


  stanvars_mult_probit
}

posterior_predict_mv_probit <- function(fit, nsamples = NULL, subset = NULL, ...) {
  subset <- brms:::subset_samples(fit, subset, nsamples)
  linpred <- posterior_linpred(fit, transform = FALSE, subset = subset, ...)

  prep <- prepare_predictions(fit, subset = subset, ...)
  N_dims <- length(prep$resps)



  rescor <- get_rescor(fit, size = N_dims, subset = subset)

  out <- array(NA_integer_, c(prep$nsamples, prep$nobs, N_dims))
  for(s in seq_len(prep$nsamples)) {
    for(o in seq_len(prep$nobs)) {
        mu_noise <- brms::rmulti_normal(1, mu = linpred[s, o, ], Sigma = rescor[s, , ])
        for(resp_id in 1:N_dims) {
          thres <- prep$resps[[resp_id]]$thres
          out[s, o, resp_id] <- sum(thres$thres[s, ] < mu_noise[resp_id])
        }
    }
  }
  out + 1
}

get_rescor <- function(fit, size, subset = NULL) {
  rescor  <- as.matrix(fit, pars = "^rescor\\[", subset = subset)

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
      k = k + 1
      out[, j, i] <- out[, i, j] <- rescor[, paste0("rescor[",k,"]")]
    }
  }
  stopifnot(all(!is.na(out)))
  out
}
