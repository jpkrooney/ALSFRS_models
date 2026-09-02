##### MV probit #####
make_stanvars_mv_probit_base <- function(column_names, rescor_prior_eta = 1) {
  N_dims <- length(column_names)

  # compose the base model
  stanvars_mult_probit_base <- stanvar(
    scode = stan_funcs_base,
    block = "functions"
  ) +
    stanvar(scode = stan_tdata_base(column_names), block = "tdata") +
    stanvar(
      scode = stan_params_base(model = "mvprobit", N_dims, nu_as_par = FALSE),
      block = "parameters"
    ) +
    stanvar(
      scode = stan_thresholds_base(N_dims, column_names),
      block = "model",
      position = "start"
    ) +
    stanvar(
      scode = stan_priors_base(
        model = "mvprobit",
        N_dims,
        nu_as_par = FALSE,
        rescor_prior_eta
      ),
      block = "model",
      position = "start"
    ) +
    stanvar(scode = stan_genquants_base(N_dims), block = "genquant")

  stanvars_mult_probit_base
}


make_stanvars_mv_probit_bgoodri <- function(
  column_names,
  rescor_prior_eta = 1,
  hetsked = FALSE
) {
  # Based on code upoaded by Ben Goodrich which uses the
  # GHK algorithm for generating TruncMVN.
  N_dims <- length(column_names)

  stanvars_base <- make_stanvars_mv_probit_base(
    column_names,
    rescor_prior_eta = rescor_prior_eta
  )

  sv <- stanvars_base +
    stanvar(scode = stan_tdata_step2(column_names), block = "tdata") +
    stanvar(scode = stan_params_step2(N_dims, hetsked), block = "parameters")
  if (hetsked == TRUE) {
    sv <- sv +
      stanvar(scode = stan_priors_scales, block = "model", position = "start")
  }
  sv +
    stanvar(
      scode = stan_likelihood(model = "mvprobit", N_dims, column_names),
      block = "likelihood",
      position = "end"
    )
}


##### MV student-t #####
make_stanvars_mv_studt_base <- function(
  column_names,
  rescor_prior_eta = 1,
  nu_as_par
) {
  N_dims <- length(column_names)

  # compose the base model
  stanvars_mult_studt_base <- stanvar(
    scode = stan_funcs_base,
    block = "functions"
  ) +
    stanvar(scode = stan_tdata_base(column_names), block = "tdata") +
    stanvar(
      scode = stan_params_base(model = "mvstudt", N_dims, nu_as_par),
      block = "parameters"
    ) +
    stanvar(
      scode = stan_thresholds_base(N_dims, column_names),
      block = "model",
      position = "start"
    ) +
    stanvar(
      scode = stan_priors_base(
        model = "mvstudt",
        N_dims,
        nu_as_par,
        rescor_prior_eta
      ),
      block = "model",
      position = "start"
    ) +
    stanvar(scode = stan_genquants_base(N_dims), block = "genquant")

  stanvars_mult_studt_base
}


make_stanvars_mv_studt <- function(
  column_names,
  rescor_prior_eta = 1,
  nu_as_par = FALSE,
  nu_beta = NULL,
  nu = NULL,
  hetsked = FALSE
) {
  #if (nu_as_par == TRUE & is.null(nu_beta)) {
  #  stop("If nu is a parameter, you must provide a value for nu_beta")
  #}
  if (nu_as_par == TRUE) {
    if (is.null(nu_beta)) {
      stop(
        "Error: nu_as_par is TRUE. You MUST provide 'nu_beta' (the hyperprior) in the data block."
      )
    }
  } else {
    if (is.null(nu)) {
      stop(
        "Error: nu_as_par is FALSE. You MUST provide 'nu' (the constant) to be passed in the data block."
      )
    }
  }

  # Based on code upoaded by Ben Goodrich which uses the
  # GHK algorithm for generating TruncMVN.
  N_dims <- length(column_names)

  stanvars_base <- make_stanvars_mv_studt_base(
    column_names,
    rescor_prior_eta = rescor_prior_eta,
    nu_as_par = nu_as_par
  )

  sv <- if (nu_as_par == TRUE) {
    # nu is a parameter, so we only pass the hyperprior to data
    sv <- stanvars_base + stanvar(x = nu_beta, name = "nu_beta", block = "data")
  } else {
    # nu is a constant, so we pass the value of nu to data
    sv <- stanvars_base + stanvar(x = nu, name = "nu", block = "data")
  }
  sv <- sv +
    stanvar(scode = stan_tdata_step2(column_names), block = "tdata") +
    stanvar(scode = stan_params_step2(N_dims, hetsked), block = "parameters")
  if (hetsked == TRUE) {
    sv <- sv +
      stanvar(scode = stan_priors_scales, block = "model", position = "start")
  }
  sv +
    stanvar(
      scode = stan_likelihood(model = "mvstudt", N_dims, column_names),
      block = "likelihood",
      position = "end"
    )

  #stanvars_mult_studt <- if (nu_as_par == FALSE) {
  #  stanvars_mult_studt <- stanvars_mult_studt_base +
  #    stanvar(x = nu, name = "nu", block = "data") +
  #    stanvar(scode = stan_tdata_step2(column_names), block = "tdata") +
  #    stanvar(scode = stan_params_step2(N_dims, hetsked), block = "parameters")
  #  if (hetsked == TRUE) {
  #    stanvars_mult_studt <- stanvars_mult_studt +
  #      stanvar(scode = stan_priors_scales, block = "model", position = "start")
  #  }
  #  stanvars_mult_studt <- stanvars_mult_studt +
  #    stanvar(
  #      scode = stan_likelihood(model = "mvstudt", N_dims, column_names),
  #      block = "likelihood",
  #      position = "end"
  #    )
  #} else {
  #  stanvars_mult_studt <- stanvars_mult_studt_base +
  #    stanvar(x = nu_beta, name = "nu_beta", block = "data") +
  #    stanvar(scode = stan_tdata_step2(column_names), block = "tdata") +
  #    stanvar(scode = stan_params_step2(N_dims, hetsked), block = "parameters")
  #  if (hetsked == TRUE) {
  #    stanvars_mult_studt <- stanvars_mult_studt +
  #      stanvar(scode = stan_priors_scales, block = "model", position = "start")
  #  }
  #  stanvars_mult_studt <- stanvars_mult_studt +
  #    stanvar(
  #      scode = stan_likelihood(model = "mvstudt", N_dims, column_names),
  #      block = "likelihood",
  #      position = "end"
  #    )
  #}

  #stanvars_mult_studt
}

#make_stanvars_mv_probit_bgoodri2 <- function(
#  column_names,
#  rescor_prior_eta = 1
#) {
#  # Based on code upoaded by Ben Goodrich which uses the
#  # GHK algorithm for generating TruncMVN.
#  N_dims <- length(column_names)
#
#  stanvars_base <- make_stanvars_mv_probit_base(
#    column_names,
#    rescor_prior_eta = rescor_prior_eta
#  )
#
#  stan_likelihood <- paste0(
#    "
#       for(n in 1:N) {
#            array[",
#    N_dims,
#    "] real mus = {",
#    paste0("mu_", column_names, "[n]", collapse = ", "),
#    "};
#            array[",
#    N_dims,
#    "] int Ys = {",
#    paste0("Y_", column_names, "[n]", collapse = ", "),
#    "};
#
#            vector[",
#    N_dims,
#    "] z;
#            real prev;
#            prev = 0;
#            for (d in 1:",
#    N_dims,
#    ") {
#              real t; // threshold at which utility = 0
#              if (Ys[d] == 1){
#                real ub = approx_Phi((thresholds[d, 1] -(mus[d] + prev)) / L_rescor[d,d]);
#                t = ub * u[n,d];
#                target += log(ub);  // Jacobian adjustment
#              } else if (Ys[d] == N_thres + 1) {
#                real lb = approx_Phi((thresholds[d, N_thres] -(mus[d] + prev)) / L_rescor[d,d]);
#                t = lb + (1 - lb) * u[n,d];
#                target += log1m(lb);  // Jacobian adjustment
#              } else {
#                real lb = approx_Phi((thresholds[d, Ys[d] - 1] -(mus[d] + prev)) / L_rescor[d,d]);
#                real ub = approx_Phi((thresholds[d, Ys[d]    ] -(mus[d] + prev)) / L_rescor[d,d]);
#                t = lb + (ub - lb) * u[n,d];
#                target += log(ub - lb);
#              }
#              z[d] = approx_inv_Phi(t);
#              if (d < ",
#    N_dims,
#    ") prev = L_rescor[d+1,1:d] * head(z, d);
#              // Jacobian adjustments imply z is truncated standard normal
#              // thus utility --- mu + L_rescor * z --- is truncated multivariate normal
#            }
#        }
#    "
#  )
#
#  stanvars_mult_probit <- stanvars_base +
#    stanvar(scode = stan_tdata_step2(column_names), block = "tdata") +
#    stanvar(scode = stan_params_step2(N_dims), block = "parameters") +
#    stanvar(scode = stan_likelihood, block = "likelihood", position = "end")
#
#  stanvars_mult_probit
#}
#
#make_stanvars_mv_studt2 <- function(
#  column_names,
#  rescor_prior_eta = 1,
#  nu_as_par = FALSE,
#  nu_beta = NULL,
#  nu = NULL
#) {
#  if (nu_as_par == TRUE & is.null(nu_beta)) {
#    stop("If nu is a parameter, you must provide a value for nu_beta")
#  }
#
#  # Based on code upoaded by Ben Goodrich which uses the
#  # GHK algorithm for generating TruncMVN.
#  N_dims <- length(column_names)
#
#  stanvars_mult_studt_base <- make_stanvars_mv_studt_base(
#    column_names,
#    rescor_prior_eta = rescor_prior_eta,
#    nu_as_par = nu_as_par
#  )
#
#  stan_likelihood <- paste0(
#    "
#  vector[N] w = (nu / 2.0) / w_raw; //student t d.f.
#
#       for(n in 1:N) {
#            array[",
#    N_dims,
#    "] real mus = {",
#    paste0("mu_", column_names, "[n]", collapse = ", "),
#    "};
#            array[",
#    N_dims,
#    "] int Ys = {",
#    paste0("Y_", column_names, "[n]", collapse = ", "),
#    "};
#
#            vector[",
#    N_dims,
#    "] z;
#            real prev;
#            prev = 0;
#            for (d in 1:",
#    N_dims,
#    ") {
#              real t; // threshold at which utility = 0
#
#              real scale_d = exp(log(L_rescor[d,d]) + log(sqrt(w[n]))); // variance for dim from residual correlaiton + student t component
#
#              if (Ys[d] == 1){
#                real ub = approx_Phi((thresholds[d, 1] - (mus[d] + prev)) / scale_d);
#                t = ub * u[n,d];
#                target += log(ub);  // Jacobian adjustment
#              } else if (Ys[d] == N_thres + 1) {
#                real lb = approx_Phi((thresholds[d, N_thres] -(mus[d] + prev)) / scale_d);
#                t = lb + (1 - lb) * u[n,d];
#                target += log1m(lb);  // Jacobian adjustment
#              } else {
#                real lb = approx_Phi((thresholds[d, Ys[d] - 1] -(mus[d] + prev)) / scale_d);
#                real ub = approx_Phi((thresholds[d, Ys[d] ] - (mus[d] + prev)) / scale_d);
#                t = lb + (ub - lb) * u[n,d];
#                target += log(ub - lb);
#              }
#              z[d] = approx_inv_Phi(t);
#              if (d < ",
#    N_dims,
#    ") prev = L_rescor[d+1,1:d] * head(z, d)* sqrt(w[n]);
#              // Jacobian adjustments imply z is truncated standard normal
#              // thus utility --- mu + L_rescor * z --- is truncated multivariate normal
#            }
#        }
#    "
#  )
#
#  stanvars_mult_studt <- if (nu_as_par == FALSE) {
#    stanvars_mult_studt_base +
#      stanvar(x = nu, name = "nu", block = "data") +
#      stanvar(scode = stan_tdata_step2(column_names), block = "tdata") +
#      stanvar(scode = stan_params_step2(N_dims), block = "parameters") +
#      stanvar(scode = stan_likelihood, block = "likelihood", position = "end")
#  } else {
#    stanvars_mult_studt_base +
#      stanvar(x = nu_beta, name = "nu_beta", block = "data") +
#      stanvar(scode = stan_tdata_step2(column_names), block = "tdata") +
#      stanvar(scode = stan_params_step2(N_dims), block = "parameters") +
#      stanvar(scode = stan_likelihood, block = "likelihood", position = "end")
#  }
#
#  stanvars_mult_studt
#}
#
