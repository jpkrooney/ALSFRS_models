make_stanvars_mv_probit_bgoodri <- function(column_names, N_cat, threshold_prior_family, threshold_prior_args) {

    N_dims <- length(column_names)


    stan_funs <- "
      real empty_ordinal_lpmf(int y, real mu) {
        return 0;
      }

      real approx_Phi(real x) {
        return inv_logit(x * 1.702);
      }

      real approx_inv_Phi(real x) {
        return logit(x) / 1.702;
      }
    "

    stan_tdata <- paste0("
       if(N != N_", column_names, ") { reject(\"Requiring equal sample size in all dimensions.\"); }" ,
                         collapse = "\n")

    stan_params <- paste0("
      cholesky_factor_corr[", N_dims, "] L_rescor;
      real<lower=0, upper=1> u[N, ", N_dims, "]; // raw residuals
      ordered[", N_cat - 1,"] thresholds[", N_dims, "];
    ")


    stan_priors <- paste0("
      target += lkj_corr_cholesky_lpdf(L_rescor | 1);
      for(q in 1:", N_dims, ") {
        target += ", threshold_prior_family, "_lpdf(thresholds[q] | ", threshold_prior_args, ");
      }
    ")

    stan_likelihood <- paste0("
        {
          for(n in 1:N) {
              real mus[", N_dims, "] = {", paste0("mu_", column_names, "[n]", collapse = ", "), "};
              int Ys[", N_dims, "] = {", paste0("Y_", column_names, "[n]", collapse = ", "), "};

              vector[", N_dims, "] z;
              real prev;
              prev = 0;
              for (d in 1:", N_dims, ") { // Phi and inv_Phi may overflow and / or be numerically inaccurate
                real t; // threshold at which utility = 0
                if (Ys[d] == 1){
                  real ub = approx_Phi((thresholds[d, 1] -(mus[d] + prev)) / L_rescor[d,d]);
                  t = ub * u[n,d];
                  target += log(ub);  // Jacobian adjustment
                } else if (Ys[d] == ", N_cat, ") {
                  real lb = approx_Phi((thresholds[d, ", N_cat - 1, "] -(mus[d] + prev)) / L_rescor[d,d]);
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
