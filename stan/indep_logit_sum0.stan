// generated with brms 2.16.3
functions {

  vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    Q_r = Q_r * inv_sqrt(1 - inv(N));
    return Q_r;
  }

  row_vector sum_to_zero_QR(row_vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    row_vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }

 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  /* cumulative-logit log-PDF for a single response
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   disc: discrimination parameter
   *   thres: ordinal thresholds
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real cumulative_logit_lpmf(int y, real mu, real disc, vector thres) {
     int nthres = num_elements(thres);
     if (y == 1) {
       return log_inv_logit(disc * (thres[1] - mu));
     } else if (y == nthres + 1) {
       return log1m_inv_logit(disc * (thres[nthres] - mu));
     } else {
       return log_diff_exp(
         log_inv_logit(disc * (thres[y] - mu)),
         log_inv_logit(disc * (thres[y - 1] - mu))
       );
     }
   }
  /* cumulative-logit log-PDF for a single response and merged thresholds
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   disc: discrimination parameter
   *   thres: vector of merged ordinal thresholds
   *   j: start and end index for the applid threshold within 'thres'
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real cumulative_logit_merged_lpmf(int y, real mu, real disc, vector thres, int[] j) {
     return cumulative_logit_lpmf(y | mu, disc, thres[j[1]:j[2]]);
   }
  /* ordered-logistic log-PDF for a single response and merged thresholds
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   thres: vector of merged ordinal thresholds
   *   j: start and end index for the applid threshold within 'thres'
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real ordered_logistic_merged_lpmf(int y, real mu, vector thres, int[] j) {
     return ordered_logistic_lpmf(y | mu, thres[j[1]:j[2]]);
   }
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_Q01;  // number of observations
  int Y_Q01[N_Q01];  // response variable
  int<lower=2> nthres_Q01;  // number of thresholds
  int<lower=1> K_Q01;  // number of population-level effects
  matrix[N_Q01, K_Q01] X_Q01;  // population-level design matrix
  int<lower=1> N_Q02;  // number of observations
  int Y_Q02[N_Q02];  // response variable
  int<lower=2> nthres_Q02;  // number of thresholds
  int<lower=1> K_Q02;  // number of population-level effects
  matrix[N_Q02, K_Q02] X_Q02;  // population-level design matrix
  int<lower=1> N_Q03;  // number of observations
  int Y_Q03[N_Q03];  // response variable
  int<lower=2> nthres_Q03;  // number of thresholds
  int<lower=1> K_Q03;  // number of population-level effects
  matrix[N_Q03, K_Q03] X_Q03;  // population-level design matrix
  int<lower=1> N_Q04;  // number of observations
  int Y_Q04[N_Q04];  // response variable
  int<lower=2> nthres_Q04;  // number of thresholds
  int<lower=1> K_Q04;  // number of population-level effects
  matrix[N_Q04, K_Q04] X_Q04;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1_Q01[N_Q01];  // grouping indicator per observation
  int<lower=1> J_1_Q02[N_Q02];  // grouping indicator per observation
  int<lower=1> J_1_Q03[N_Q03];  // grouping indicator per observation
  int<lower=1> J_1_Q04[N_Q04];  // grouping indicator per observation
  // group-level predictor values
  vector[N_Q01] Z_1_Q01_1;
  vector[N_Q02] Z_1_Q02_2;
  vector[N_Q03] Z_1_Q03_3;
  vector[N_Q04] Z_1_Q04_4;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_Q01 = K_Q01;
  matrix[N_Q01, Kc_Q01] Xc_Q01;  // centered version of X_Q01
  vector[Kc_Q01] means_X_Q01;  // column means of X_Q01 before centering
  int Kc_Q02 = K_Q02;
  matrix[N_Q02, Kc_Q02] Xc_Q02;  // centered version of X_Q02
  vector[Kc_Q02] means_X_Q02;  // column means of X_Q02 before centering
  int Kc_Q03 = K_Q03;
  matrix[N_Q03, Kc_Q03] Xc_Q03;  // centered version of X_Q03
  vector[Kc_Q03] means_X_Q03;  // column means of X_Q03 before centering
  int Kc_Q04 = K_Q04;
  matrix[N_Q04, Kc_Q04] Xc_Q04;  // centered version of X_Q04
  vector[Kc_Q04] means_X_Q04;  // column means of X_Q04 before centering

  vector[2 * N_1] groups_Q_r = Q_sum_to_zero_QR(N_1);

  for (i in 1:K_Q01) {
    means_X_Q01[i] = mean(X_Q01[, i]);
    Xc_Q01[, i] = X_Q01[, i] - means_X_Q01[i];
  }
  for (i in 1:K_Q02) {
    means_X_Q02[i] = mean(X_Q02[, i]);
    Xc_Q02[, i] = X_Q02[, i] - means_X_Q02[i];
  }
  for (i in 1:K_Q03) {
    means_X_Q03[i] = mean(X_Q03[, i]);
    Xc_Q03[, i] = X_Q03[, i] - means_X_Q03[i];
  }
  for (i in 1:K_Q04) {
    means_X_Q04[i] = mean(X_Q04[, i]);
    Xc_Q04[, i] = X_Q04[, i] - means_X_Q04[i];
  }
}
parameters {
  vector[Kc_Q01] b_Q01;  // population-level effects
  ordered[nthres_Q01] Intercept_Q01;  // temporary thresholds for centered predictors
  vector[Kc_Q02] b_Q02;  // population-level effects
  ordered[nthres_Q02] Intercept_Q02;  // temporary thresholds for centered predictors
  vector[Kc_Q03] b_Q03;  // population-level effects
  ordered[nthres_Q03] Intercept_Q03;  // temporary thresholds for centered predictors
  vector[Kc_Q04] b_Q04;  // population-level effects
  ordered[nthres_Q04] Intercept_Q04;  // temporary thresholds for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  matrix[M_1, N_1 - 1] z_1_s0;  // standardized group-level effects

}
transformed parameters {
  real<lower=0> disc_Q01 = 1;  // discrimination parameters
  real<lower=0> disc_Q02 = 1;  // discrimination parameters
  real<lower=0> disc_Q03 = 1;  // discrimination parameters
  real<lower=0> disc_Q04 = 1;  // discrimination parameters
  matrix[M_1, N_1] z_1;  // standardized group-level effects

  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_Q01_1;
  vector[N_1] r_1_Q02_2;
  vector[N_1] r_1_Q03_3;
  vector[N_1] r_1_Q04_4;

  for(m in 1:M_1) {
    z_1[m, ] = sum_to_zero_QR(z_1_s0[m, ], groups_Q_r);
  }
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_Q01_1 = r_1[, 1];
  r_1_Q02_2 = r_1[, 2];
  r_1_Q03_3 = r_1[, 3];
  r_1_Q04_4 = r_1[, 4];
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_Q01] mu_Q01 = Xc_Q01 * b_Q01;
    // initialize linear predictor term
    vector[N_Q02] mu_Q02 = Xc_Q02 * b_Q02;
    // initialize linear predictor term
    vector[N_Q03] mu_Q03 = Xc_Q03 * b_Q03;
    // initialize linear predictor term
    vector[N_Q04] mu_Q04 = Xc_Q04 * b_Q04;
    for (n in 1:N_Q01) {
      // add more terms to the linear predictor
      mu_Q01[n] += r_1_Q01_1[J_1_Q01[n]] * Z_1_Q01_1[n];
    }
    for (n in 1:N_Q02) {
      // add more terms to the linear predictor
      mu_Q02[n] += r_1_Q02_2[J_1_Q02[n]] * Z_1_Q02_2[n];
    }
    for (n in 1:N_Q03) {
      // add more terms to the linear predictor
      mu_Q03[n] += r_1_Q03_3[J_1_Q03[n]] * Z_1_Q03_3[n];
    }
    for (n in 1:N_Q04) {
      // add more terms to the linear predictor
      mu_Q04[n] += r_1_Q04_4[J_1_Q04[n]] * Z_1_Q04_4[n];
    }
    for (n in 1:N_Q01) {
      target += ordered_logistic_lpmf(Y_Q01[n] | mu_Q01[n], Intercept_Q01);
    }
    for (n in 1:N_Q02) {
      target += ordered_logistic_lpmf(Y_Q02[n] | mu_Q02[n], Intercept_Q02);
    }
    for (n in 1:N_Q03) {
      target += ordered_logistic_lpmf(Y_Q03[n] | mu_Q03[n], Intercept_Q03);
    }
    for (n in 1:N_Q04) {
      target += ordered_logistic_lpmf(Y_Q04[n] | mu_Q04[n], Intercept_Q04);
    }
  }
  // priors including constants
  target += student_t_lpdf(Intercept_Q01 | 3, 0, 2.5);
  target += student_t_lpdf(Intercept_Q02 | 3, 0, 2.5);
  target += student_t_lpdf(Intercept_Q03 | 3, 0, 2.5);
  target += student_t_lpdf(Intercept_Q04 | 3, 0, 2.5);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 4 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
}
generated quantities {
  // compute actual thresholds
  vector[nthres_Q01] b_Q01_Intercept = Intercept_Q01 + dot_product(means_X_Q01, b_Q01);
  // compute actual thresholds
  vector[nthres_Q02] b_Q02_Intercept = Intercept_Q02 + dot_product(means_X_Q02, b_Q02);
  // compute actual thresholds
  vector[nthres_Q03] b_Q03_Intercept = Intercept_Q03 + dot_product(means_X_Q03, b_Q03);
  // compute actual thresholds
  vector[nthres_Q04] b_Q04_Intercept = Intercept_Q04 + dot_product(means_X_Q04, b_Q04);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
