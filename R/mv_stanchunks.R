#### Stan code is assembled in two steps ####
# The first step enables multivariate probit or multivariate student t copula
# The second step allows for homoskedastic or heteroskedastic correlations between multivariate outcome variables

#### Define base model code chunks and functions ####
stan_funcs_base <- "
      real empty_cumulative_lpmf(int y, real mu, vector intercept) {
        return 0;
      }

      real approx_Phi(real x) {
        return inv_logit(x * 1.702);
      }

      real approx_inv_Phi(real x) {
        return logit(x) / 1.702;
      }
    "

stan_tdata_base <- function(column_names) {
  paste0(
    "
           int N_thres = nthres_",
    column_names[1],
    ";",

    paste0(
      "
         if(N != N_",
      column_names,
      ") { reject(\"Requiring equal sample size in all dimensions.\"); }",
      collapse = "\n"
    )
  )
}


stan_params_base <- function(
  model = c("mvprobit", "mvstudt"),
  N_dims,
  nu_as_par = FALSE
) {
  # define optional code lines
  optioncode <- ifelse(
    model == "mvprobit",
    "",
    ifelse(
      nu_as_par == FALSE,
      "  vector<lower=0> [N] w_raw;\n",
      "  real<lower=2> nu; // student-t df lower bound of 2 to ensure finite variance
  vector<lower=0> [N] w_raw;\n"
    )
  )

  # compose stan code
  stan_string <- glue(
    "\n
      cholesky_factor_corr[ {N_dims} ] L_rescor;
    {optioncode}
    "
  )
  stan_string
}


stan_priors_base <- function(
  model = c("mvprobit", "mvstudt"),
  N_dims,
  nu_as_par = FALSE,
  rescor_prior_eta
) {
  # define optional code lines
  optioncode <- ifelse(
    model == "mvprobit",
    "",
    ifelse(
      nu_as_par == FALSE,
      "// student t degrees of freedom parameter prior
  target += gamma_lpdf(w_raw | nu/2, nu/2);",
      "// student t degrees of freedom parameter prior
  target += gamma_lpdf(w_raw | nu/2, nu/2);
  target += gamma_lpdf(nu | 2, nu_beta);\n"
    )
  )

  glue(
    "\n
      target += lkj_corr_cholesky_lpdf(L_rescor | {rescor_prior_eta} );\n
      {optioncode}
    "
  )
}


stan_thresholds_base <- function(N_dims, column_names) {
  paste0(
    "
         array [",
    N_dims,
    "] vector[N_thres] thresholds;",
    paste0(
      "
         thresholds[",
      1:N_dims,
      "] = Intercept_",
      column_names,
      ";",
      collapse = ""
    )
  )
}

stan_genquants_base <- function(N_dims) {
  paste0(
    "
     corr_matrix[",
    N_dims,
    "] Rescor = multiply_lower_tri_self_transpose(L_rescor);
     vector<lower=-1,upper=1>[",
    choose(N_dims, 2),
    "] rescor;
     // extract upper diagonal of rescor matrix
     for (k in 1:",
    N_dims,
    ") {
        for (j in 1:(k - 1)) {
          rescor[choose(k - 1, 2) + j] = Rescor[j, k];
        }
      }
    "
  )
}


#### Second step stan chunks ####
stan_tdata_step2 <- function(column_names) {
  paste0(
    "
         if(nthres_",
    column_names[1],
    " != nthres_",
    column_names[2:length(column_names)],
    ") { reject(\"Requiring equal number of thresholds in all dimensions.\"); }",
    collapse = "\n"
  )
}

stan_params_step2 <- function(N_dims, hetsked = FALSE) {
  if (hetsked == FALSE) {
    paste0(
      "
      array[N, ",
      N_dims,
      "] real<lower=0, upper=1> u; // raw residuals
    "
    )
  } else {
    paste0(
      "
      array[N, ",
      N_dims,
      "] real<lower=0, upper=1> u; // raw residuals
    vector[",
      N_dims,
      "] scale_intercept; // The baseline scale for each dimension
    vector[",
      N_dims,
      "] scale_slope;     // The rate of change of scale over dimensions

    "
    )
  }
}

stan_priors_scales <- paste0(
  "
  target += normal_lpdf(scale_intercept | 0, 0.5);
  target += normal_lpdf(scale_slope | 0, 0.5);
      "
)


#### Main stan likelihood code ####
stan_likelihood <- function(
  model = c("mvprobit", "mvstudt"),
  N_dims,
  column_names,
  hetsked = FALSE
) {
  # define optional code lines
  wparam_chk <- ifelse(
    model == "mvstudt",
    "  vector[N] w = (nu / 2.0) / w_raw; //student t d.f.
",
    ""
  )
  scale_chk <- if (model == "mvprobit" & hetsked == FALSE) {
    "real scale_d = L_rescor[d,d];"
  } else if (model == "mvprobit" & hetsked == TRUE) {
    "real scale_d = exp(log(L_rescor[d,d]) + log(scale_intercept[d] - scale_slope[d]* exp(- X_Q01[n, 1]));"
  } else if (model == "mvstudt" & hetsked == FALSE) {
    "real scale_d = exp(log(L_rescor[d,d]) + log(sqrt(w[n]))); // variance for dim from residual correlaiton + student t component"
  } else if (model == "mvstudt" & hetsked == TRUE) {
    "real scale_d = exp(log(L_rescor[d,d]) + log(scale_intercept[d] - scale_slope[d]* exp(- X_Q01[n, 1]) +
      log(sqrt(w[n]))); // variance for dim from residual correlaiton + student t component"
  }

  # define required code chunks
  initial_chk <- paste0(
    "
       for(n in 1:N) {
            array[",
    N_dims,
    "] real mus = {",
    paste0("mu_", column_names, "[n]", collapse = ", "),
    "};
            array[",
    N_dims,
    "] int Ys = {",
    paste0("Y_", column_names, "[n]", collapse = ", "),
    "};

            vector[",
    N_dims,
    "] z;
            real prev;
            prev = 0;
            for (d in 1:",
    N_dims,
    ") {
              real t; // threshold at which utility = 0
              "
  )

  mvloop_chk <- paste0(
    "
                if (Ys[d] == 1){
                real ub = approx_Phi((thresholds[d, 1] -(mus[d] + prev)) / scale_d);
                t = ub * u[n,d];
                target += log(ub);  // Jacobian adjustment
              } else if (Ys[d] == N_thres + 1) {
                real lb = approx_Phi((thresholds[d, N_thres] -(mus[d] + prev)) / scale_d);
                t = lb + (1 - lb) * u[n,d];
                target += log1m(lb);  // Jacobian adjustment
              } else {
                real lb = approx_Phi((thresholds[d, Ys[d] - 1] -(mus[d] + prev)) / scale_d);
                real ub = approx_Phi((thresholds[d, Ys[d]    ] -(mus[d] + prev)) / scale_d);
                t = lb + (ub - lb) * u[n,d];
                target += log(ub - lb);
              }
              z[d] = approx_inv_Phi(t);
              if (d < ",
    N_dims,
    ") prev = L_rescor[d+1,1:d] * head(z, d);
              // Jacobian adjustments imply z is truncated standard normal
              // thus utility --- mu + L_rescor * z --- is truncated multivariate normal
      }
  }
  "
  )

  # compose chunks
  paste0(wparam_chk, initial_chk, scale_chk, mvloop_chk)
}
