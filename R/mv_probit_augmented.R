empty_ordinal <- function(N_cat) {
  custom_family(
    "empty_ordinal", dpars = c("mu"),
    links = c("identity"), lb = 1, ub = N_cat,
    type = "int")
}

make_stanvars_ordinal_augmented <- function(column_names, N_cat, threshold_prior_family, threshold_prior_args) {

    N_dims <- length(column_names)


    stan_funs <- "
      real empty_ordinal_lpmf(int y, real mu) {
        return 0;
      }

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

    stan_tdata <- paste0("
       if(N != N_", column_names, ") { reject(\"Requiring equal sample size in all dimensions.\"); }" ,
                         collapse = "\n")

    stan_params <- paste0("
      cholesky_factor_corr[", N_dims, "] L_rescor;
      real z_rescor[N, ", N_dims, "]; // raw residuals
      ordered[", N_cat - 1,"] thresholds[", N_dims, "];
    ")


    stan_priors <- paste0("
      target += lkj_corr_cholesky_lpdf(L_rescor | 1);
      for(q in 1:", N_dims, ") {
        target += ", threshold_prior_family, "_lpdf(thresholds[q] | ", threshold_prior_args, ");
      }
    ")

    constrain_residuals_code <- paste0(
        "residuals[n, ", 1:N_dims, "] = constrain_residual_lp(z_rescor[n, ", 1:N_dims, "], mu_", column_names, "[n], Y_", column_names, "[n], thresholds[", 1:N_dims,"]);",    collapse = "\n          ")

    stan_likelihood <- paste0("
        {
          vector[", N_dims, "] residuals[N];
          for(n in 1:N) {
              ", constrain_residuals_code , "
          }
          target += multi_normal_cholesky_lpdf(residuals | rep_vector(0, ", N_dims, "), L_rescor);
        }
    ")

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
        stanvar(scode = stan_priors, block = "model", position = "start") +
        stanvar(scode = stan_likelihood, block = "likelihood", position = "end") +
        stanvar(scode = stan_genquants, block = "genquant")


    stanvars_mult_probit
}
