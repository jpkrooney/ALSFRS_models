// generated with brms 2.23.1
functions {
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

      real empty_cumulative_lpmf(int y, real mu, vector intercept) {
        return 0;
      }

      real approx_Phi(real x) {
        return inv_logit(x * 1.702);
      }

      real approx_inv_Phi(real x) {
        return logit(x) / 1.702;
      }



      real jm_mvprobit_cumhazard_gk(
                    row_vector Bobs,
                    matrix Bquad,
                    row_vector quad_times,
                    vector GKw,
                    vector theta,
                    real time,
                    //real Z_1_1,
                    real beta_surv,
                    real x,
                    vector assoc, //surv params
                    //real trt,
                    vector b_q01, real rInt_q01,// long params
                    vector b_q02, real rInt_q02,
                    vector b_q03, real rInt_q03,
                    vector b_q04, real rInt_q04,
                    vector b_q05, real rInt_q05,
                    vector b_q06, real rInt_q06,
                    vector b_q07, real rInt_q07,
                    vector b_q08, real rInt_q08,
                    vector b_q09, real rInt_q09,
                    vector b_q10, real rInt_q10,
                    vector b_q11, real rInt_q11,
                    vector b_q12, real rInt_q12) {
        // assumes effects constant over time
        int Q = rows(Bquad);
        real H = 0;

        for(q in 1:Q){
            real h0_q = dot_product( theta, to_vector(Bquad[q]')); //baseline hazard
            real m_q01 = eval_linpred(quad_times[q], x, b_q01, rInt_q01);
            //real m_q02 = eval_linpred(quad_times[q], x, b_q02, rInt_q02);
            //real m_q03 = eval_linpred(quad_times[q], x, b_q03, rInt_q03);
            //real m_q04 = eval_linpred(quad_times[q], x, b_q04, rInt_q04);
            //real m_q05 = eval_linpred(quad_times[q], x, b_q05, rInt_q05);
            //real m_q06 = eval_linpred(quad_times[q], x, b_q06, rInt_q06);
            //real m_q07 = eval_linpred(quad_times[q], x, b_q07, rInt_q07);
            //real m_q08 = eval_linpred(quad_times[q], x, b_q08, rInt_q08);
            //real m_q09 = eval_linpred(quad_times[q], x, b_q09, rInt_q09);
            //real m_q10 = eval_linpred(quad_times[q], x, b_q10, rInt_q10);
            //real m_q11 = eval_linpred(quad_times[q], x, b_q11, rInt_q11);
            //real m_q12 = eval_linpred(quad_times[q], x, b_q12, rInt_q12);

            // current value association
            H += GKw[q] * h0_q * exp(beta_surv * x + assoc[1]*m_q01);// + assoc[2]*m_q02 +
                //assoc[3] *m_q03 + assoc[4] *m_q04 + assoc[5] *m_q05 + assoc[6] *m_q06 +
                //assoc[7] *m_q07 + assoc[8] *m_q08 + assoc[9] *m_q09 + assoc[10] *m_q10 +
                //assoc[11] *m_q11 + assoc[12] *m_q12); // need multiple m_q's
        }
        H *= time/2;
        return H;
    }

    real eval_linpred(real time, real trt, vector b, real rInt){
        // population level effects
        real linpred = time * b[1] + trt* b[2] + time * trt * b[3] + rInt * 1;
        return linpred;
    }


}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_Q01;  // number of observations
  array[N_Q01] int Y_Q01;  // response variable
  int<lower=2> nthres_Q01;  // number of thresholds
  int<lower=1> K_Q01;  // number of population-level effects
  matrix[N_Q01, K_Q01] X_Q01;  // population-level design matrix
  int<lower=1> Kc_Q01;  // number of population-level effects after centering
  int<lower=1> N_Q02;  // number of observations
  array[N_Q02] int Y_Q02;  // response variable
  int<lower=2> nthres_Q02;  // number of thresholds
  int<lower=1> K_Q02;  // number of population-level effects
  matrix[N_Q02, K_Q02] X_Q02;  // population-level design matrix
  int<lower=1> Kc_Q02;  // number of population-level effects after centering
  int<lower=1> N_Q03;  // number of observations
  array[N_Q03] int Y_Q03;  // response variable
  int<lower=2> nthres_Q03;  // number of thresholds
  int<lower=1> K_Q03;  // number of population-level effects
  matrix[N_Q03, K_Q03] X_Q03;  // population-level design matrix
  int<lower=1> Kc_Q03;  // number of population-level effects after centering
  int<lower=1> N_Q04;  // number of observations
  array[N_Q04] int Y_Q04;  // response variable
  int<lower=2> nthres_Q04;  // number of thresholds
  int<lower=1> K_Q04;  // number of population-level effects
  matrix[N_Q04, K_Q04] X_Q04;  // population-level design matrix
  int<lower=1> Kc_Q04;  // number of population-level effects after centering
  int<lower=1> N_Q05;  // number of observations
  array[N_Q05] int Y_Q05;  // response variable
  int<lower=2> nthres_Q05;  // number of thresholds
  int<lower=1> K_Q05;  // number of population-level effects
  matrix[N_Q05, K_Q05] X_Q05;  // population-level design matrix
  int<lower=1> Kc_Q05;  // number of population-level effects after centering
  int<lower=1> N_Q06;  // number of observations
  array[N_Q06] int Y_Q06;  // response variable
  int<lower=2> nthres_Q06;  // number of thresholds
  int<lower=1> K_Q06;  // number of population-level effects
  matrix[N_Q06, K_Q06] X_Q06;  // population-level design matrix
  int<lower=1> Kc_Q06;  // number of population-level effects after centering
  int<lower=1> N_Q07;  // number of observations
  array[N_Q07] int Y_Q07;  // response variable
  int<lower=2> nthres_Q07;  // number of thresholds
  int<lower=1> K_Q07;  // number of population-level effects
  matrix[N_Q07, K_Q07] X_Q07;  // population-level design matrix
  int<lower=1> Kc_Q07;  // number of population-level effects after centering
  int<lower=1> N_Q08;  // number of observations
  array[N_Q08] int Y_Q08;  // response variable
  int<lower=2> nthres_Q08;  // number of thresholds
  int<lower=1> K_Q08;  // number of population-level effects
  matrix[N_Q08, K_Q08] X_Q08;  // population-level design matrix
  int<lower=1> Kc_Q08;  // number of population-level effects after centering
  int<lower=1> N_Q09;  // number of observations
  array[N_Q09] int Y_Q09;  // response variable
  int<lower=2> nthres_Q09;  // number of thresholds
  int<lower=1> K_Q09;  // number of population-level effects
  matrix[N_Q09, K_Q09] X_Q09;  // population-level design matrix
  int<lower=1> Kc_Q09;  // number of population-level effects after centering
  int<lower=1> N_Q10;  // number of observations
  array[N_Q10] int Y_Q10;  // response variable
  int<lower=2> nthres_Q10;  // number of thresholds
  int<lower=1> K_Q10;  // number of population-level effects
  matrix[N_Q10, K_Q10] X_Q10;  // population-level design matrix
  int<lower=1> Kc_Q10;  // number of population-level effects after centering
  int<lower=1> N_Q11;  // number of observations
  array[N_Q11] int Y_Q11;  // response variable
  int<lower=2> nthres_Q11;  // number of thresholds
  int<lower=1> K_Q11;  // number of population-level effects
  matrix[N_Q11, K_Q11] X_Q11;  // population-level design matrix
  int<lower=1> Kc_Q11;  // number of population-level effects after centering
  int<lower=1> N_Q12;  // number of observations
  array[N_Q12] int Y_Q12;  // response variable
  int<lower=2> nthres_Q12;  // number of thresholds
  int<lower=1> K_Q12;  // number of population-level effects
  matrix[N_Q12, K_Q12] X_Q12;  // population-level design matrix
  int<lower=1> Kc_Q12;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N_Q01] int<lower=1> J_1_Q01;  // grouping indicator per observation
  array[N_Q02] int<lower=1> J_1_Q02;  // grouping indicator per observation
  array[N_Q03] int<lower=1> J_1_Q03;  // grouping indicator per observation
  array[N_Q04] int<lower=1> J_1_Q04;  // grouping indicator per observation
  array[N_Q05] int<lower=1> J_1_Q05;  // grouping indicator per observation
  array[N_Q06] int<lower=1> J_1_Q06;  // grouping indicator per observation
  array[N_Q07] int<lower=1> J_1_Q07;  // grouping indicator per observation
  array[N_Q08] int<lower=1> J_1_Q08;  // grouping indicator per observation
  array[N_Q09] int<lower=1> J_1_Q09;  // grouping indicator per observation
  array[N_Q10] int<lower=1> J_1_Q10;  // grouping indicator per observation
  array[N_Q11] int<lower=1> J_1_Q11;  // grouping indicator per observation
  array[N_Q12] int<lower=1> J_1_Q12;  // grouping indicator per observation
  // group-level predictor values
  vector[N_Q01] Z_1_Q01_1;
  vector[N_Q02] Z_1_Q02_2;
  vector[N_Q03] Z_1_Q03_3;
  vector[N_Q04] Z_1_Q04_4;
  vector[N_Q05] Z_1_Q05_5;
  vector[N_Q06] Z_1_Q06_6;
  vector[N_Q07] Z_1_Q07_7;
  vector[N_Q08] Z_1_Q08_8;
  vector[N_Q09] Z_1_Q09_9;
  vector[N_Q10] Z_1_Q10_10;
  vector[N_Q11] Z_1_Q11_11;
  vector[N_Q12] Z_1_Q12_12;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?


  /// survival data
  vector[N_1] event_times;
  vector[N_1] event;
  vector[N_1] Xsurv;
  int<lower=1> K;
  matrix[N_1,K] B_obs;
  int<lower=1> Q;
  array[N_1] matrix[Q,K] B_quad;
  matrix[N_1, Q] quad_times;
  vector[Q] GK_weights;
}
transformed data {
  matrix[N_Q01, Kc_Q01] Xc_Q01;  // centered version of X_Q01
  vector[Kc_Q01] means_X_Q01;  // column means of X_Q01 before centering
  matrix[N_Q02, Kc_Q02] Xc_Q02;  // centered version of X_Q02
  vector[Kc_Q02] means_X_Q02;  // column means of X_Q02 before centering
  matrix[N_Q03, Kc_Q03] Xc_Q03;  // centered version of X_Q03
  vector[Kc_Q03] means_X_Q03;  // column means of X_Q03 before centering
  matrix[N_Q04, Kc_Q04] Xc_Q04;  // centered version of X_Q04
  vector[Kc_Q04] means_X_Q04;  // column means of X_Q04 before centering
  matrix[N_Q05, Kc_Q05] Xc_Q05;  // centered version of X_Q05
  vector[Kc_Q05] means_X_Q05;  // column means of X_Q05 before centering
  matrix[N_Q06, Kc_Q06] Xc_Q06;  // centered version of X_Q06
  vector[Kc_Q06] means_X_Q06;  // column means of X_Q06 before centering
  matrix[N_Q07, Kc_Q07] Xc_Q07;  // centered version of X_Q07
  vector[Kc_Q07] means_X_Q07;  // column means of X_Q07 before centering
  matrix[N_Q08, Kc_Q08] Xc_Q08;  // centered version of X_Q08
  vector[Kc_Q08] means_X_Q08;  // column means of X_Q08 before centering
  matrix[N_Q09, Kc_Q09] Xc_Q09;  // centered version of X_Q09
  vector[Kc_Q09] means_X_Q09;  // column means of X_Q09 before centering
  matrix[N_Q10, Kc_Q10] Xc_Q10;  // centered version of X_Q10
  vector[Kc_Q10] means_X_Q10;  // column means of X_Q10 before centering
  matrix[N_Q11, Kc_Q11] Xc_Q11;  // centered version of X_Q11
  vector[Kc_Q11] means_X_Q11;  // column means of X_Q11 before centering
  matrix[N_Q12, Kc_Q12] Xc_Q12;  // centered version of X_Q12
  vector[Kc_Q12] means_X_Q12;  // column means of X_Q12 before centering

           int N_thres = nthres_Q01;
         if(N != N_Q01) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q02) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q03) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q04) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q05) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q06) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q07) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q08) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q09) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q10) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q11) { reject("Requiring equal sample size in all dimensions."); }

         if(N != N_Q12) { reject("Requiring equal sample size in all dimensions."); }

         if(nthres_Q01 != nthres_Q02) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q03) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q04) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q05) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q06) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q07) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q08) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q09) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q10) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q11) { reject("Requiring equal number of thresholds in all dimensions."); }

         if(nthres_Q01 != nthres_Q12) { reject("Requiring equal number of thresholds in all dimensions."); }
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
  for (i in 1:K_Q05) {
    means_X_Q05[i] = mean(X_Q05[, i]);
    Xc_Q05[, i] = X_Q05[, i] - means_X_Q05[i];
  }
  for (i in 1:K_Q06) {
    means_X_Q06[i] = mean(X_Q06[, i]);
    Xc_Q06[, i] = X_Q06[, i] - means_X_Q06[i];
  }
  for (i in 1:K_Q07) {
    means_X_Q07[i] = mean(X_Q07[, i]);
    Xc_Q07[, i] = X_Q07[, i] - means_X_Q07[i];
  }
  for (i in 1:K_Q08) {
    means_X_Q08[i] = mean(X_Q08[, i]);
    Xc_Q08[, i] = X_Q08[, i] - means_X_Q08[i];
  }
  for (i in 1:K_Q09) {
    means_X_Q09[i] = mean(X_Q09[, i]);
    Xc_Q09[, i] = X_Q09[, i] - means_X_Q09[i];
  }
  for (i in 1:K_Q10) {
    means_X_Q10[i] = mean(X_Q10[, i]);
    Xc_Q10[, i] = X_Q10[, i] - means_X_Q10[i];
  }
  for (i in 1:K_Q11) {
    means_X_Q11[i] = mean(X_Q11[, i]);
    Xc_Q11[, i] = X_Q11[, i] - means_X_Q11[i];
  }
  for (i in 1:K_Q12) {
    means_X_Q12[i] = mean(X_Q12[, i]);
    Xc_Q12[, i] = X_Q12[, i] - means_X_Q12[i];
  }

}
parameters {
  vector[Kc_Q01] b_Q01;  // regression coefficients
  ordered[nthres_Q01] Intercept_Q01;  // temporary thresholds for centered predictors
  vector[Kc_Q02] b_Q02;  // regression coefficients
  ordered[nthres_Q02] Intercept_Q02;  // temporary thresholds for centered predictors
  vector[Kc_Q03] b_Q03;  // regression coefficients
  ordered[nthres_Q03] Intercept_Q03;  // temporary thresholds for centered predictors
  vector[Kc_Q04] b_Q04;  // regression coefficients
  ordered[nthres_Q04] Intercept_Q04;  // temporary thresholds for centered predictors
  vector[Kc_Q05] b_Q05;  // regression coefficients
  ordered[nthres_Q05] Intercept_Q05;  // temporary thresholds for centered predictors
  vector[Kc_Q06] b_Q06;  // regression coefficients
  ordered[nthres_Q06] Intercept_Q06;  // temporary thresholds for centered predictors
  vector[Kc_Q07] b_Q07;  // regression coefficients
  ordered[nthres_Q07] Intercept_Q07;  // temporary thresholds for centered predictors
  vector[Kc_Q08] b_Q08;  // regression coefficients
  ordered[nthres_Q08] Intercept_Q08;  // temporary thresholds for centered predictors
  vector[Kc_Q09] b_Q09;  // regression coefficients
  ordered[nthres_Q09] Intercept_Q09;  // temporary thresholds for centered predictors
  vector[Kc_Q10] b_Q10;  // regression coefficients
  ordered[nthres_Q10] Intercept_Q10;  // temporary thresholds for centered predictors
  vector[Kc_Q11] b_Q11;  // regression coefficients
  ordered[nthres_Q11] Intercept_Q11;  // temporary thresholds for centered predictors
  vector[Kc_Q12] b_Q12;  // regression coefficients
  ordered[nthres_Q12] Intercept_Q12;  // temporary thresholds for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix

      cholesky_factor_corr[12] L_rescor;


      array[N, 12] real<lower=0, upper=1> u; // raw residuals


          // survival parameters
    real beta_surv;
    real eta0;
    simplex[K] w;
    //real assoc;
    vector[1] assoc; // increase size if wish to use more longituudinal submodels

}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_Q01_1;
  vector[N_1] r_1_Q02_2;
  vector[N_1] r_1_Q03_3;
  vector[N_1] r_1_Q04_4;
  vector[N_1] r_1_Q05_5;
  vector[N_1] r_1_Q06_6;
  vector[N_1] r_1_Q07_7;
  vector[N_1] r_1_Q08_8;
  vector[N_1] r_1_Q09_9;
  vector[N_1] r_1_Q10_10;
  vector[N_1] r_1_Q11_11;
  vector[N_1] r_1_Q12_12;
  // prior contributions to the log posterior
  real lprior = 0;
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_Q01_1 = r_1[, 1];
  r_1_Q02_2 = r_1[, 2];
  r_1_Q03_3 = r_1[, 3];
  r_1_Q04_4 = r_1[, 4];
  r_1_Q05_5 = r_1[, 5];
  r_1_Q06_6 = r_1[, 6];
  r_1_Q07_7 = r_1[, 7];
  r_1_Q08_8 = r_1[, 8];
  r_1_Q09_9 = r_1[, 9];
  r_1_Q10_10 = r_1[, 10];
  r_1_Q11_11 = r_1[, 11];
  r_1_Q12_12 = r_1[, 12];
  lprior += student_t_lpdf(Intercept_Q01 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q02 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q03 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q04 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q05 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q06 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q07 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q08 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q09 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q10 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q11 | 3, 0, 2.5);
  lprior += student_t_lpdf(Intercept_Q12 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 12 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);


    // survival transformed params
  vector[K] theta;
  theta = exp(eta0) * w;

}
model {

         array [12] vector[N_thres] thresholds;
         thresholds[1] = Intercept_Q01;
         thresholds[2] = Intercept_Q02;
         thresholds[3] = Intercept_Q03;
         thresholds[4] = Intercept_Q04;
         thresholds[5] = Intercept_Q05;
         thresholds[6] = Intercept_Q06;
         thresholds[7] = Intercept_Q07;
         thresholds[8] = Intercept_Q08;
         thresholds[9] = Intercept_Q09;
         thresholds[10] = Intercept_Q10;
         thresholds[11] = Intercept_Q11;
         thresholds[12] = Intercept_Q12;

      target += lkj_corr_cholesky_lpdf(L_rescor | 1);

  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_Q01] mu_Q01 = rep_vector(0.0, N_Q01);
    // initialize linear predictor term
    vector[N_Q02] mu_Q02 = rep_vector(0.0, N_Q02);
    // initialize linear predictor term
    vector[N_Q03] mu_Q03 = rep_vector(0.0, N_Q03);
    // initialize linear predictor term
    vector[N_Q04] mu_Q04 = rep_vector(0.0, N_Q04);
    // initialize linear predictor term
    vector[N_Q05] mu_Q05 = rep_vector(0.0, N_Q05);
    // initialize linear predictor term
    vector[N_Q06] mu_Q06 = rep_vector(0.0, N_Q06);
    // initialize linear predictor term
    vector[N_Q07] mu_Q07 = rep_vector(0.0, N_Q07);
    // initialize linear predictor term
    vector[N_Q08] mu_Q08 = rep_vector(0.0, N_Q08);
    // initialize linear predictor term
    vector[N_Q09] mu_Q09 = rep_vector(0.0, N_Q09);
    // initialize linear predictor term
    vector[N_Q10] mu_Q10 = rep_vector(0.0, N_Q10);
    // initialize linear predictor term
    vector[N_Q11] mu_Q11 = rep_vector(0.0, N_Q11);
    // initialize linear predictor term
    vector[N_Q12] mu_Q12 = rep_vector(0.0, N_Q12);
    mu_Q01 += Xc_Q01 * b_Q01;
    mu_Q02 += Xc_Q02 * b_Q02;
    mu_Q03 += Xc_Q03 * b_Q03;
    mu_Q04 += Xc_Q04 * b_Q04;
    mu_Q05 += Xc_Q05 * b_Q05;
    mu_Q06 += Xc_Q06 * b_Q06;
    mu_Q07 += Xc_Q07 * b_Q07;
    mu_Q08 += Xc_Q08 * b_Q08;
    mu_Q09 += Xc_Q09 * b_Q09;
    mu_Q10 += Xc_Q10 * b_Q10;
    mu_Q11 += Xc_Q11 * b_Q11;
    mu_Q12 += Xc_Q12 * b_Q12;
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
    for (n in 1:N_Q05) {
      // add more terms to the linear predictor
      mu_Q05[n] += r_1_Q05_5[J_1_Q05[n]] * Z_1_Q05_5[n];
    }
    for (n in 1:N_Q06) {
      // add more terms to the linear predictor
      mu_Q06[n] += r_1_Q06_6[J_1_Q06[n]] * Z_1_Q06_6[n];
    }
    for (n in 1:N_Q07) {
      // add more terms to the linear predictor
      mu_Q07[n] += r_1_Q07_7[J_1_Q07[n]] * Z_1_Q07_7[n];
    }
    for (n in 1:N_Q08) {
      // add more terms to the linear predictor
      mu_Q08[n] += r_1_Q08_8[J_1_Q08[n]] * Z_1_Q08_8[n];
    }
    for (n in 1:N_Q09) {
      // add more terms to the linear predictor
      mu_Q09[n] += r_1_Q09_9[J_1_Q09[n]] * Z_1_Q09_9[n];
    }
    for (n in 1:N_Q10) {
      // add more terms to the linear predictor
      mu_Q10[n] += r_1_Q10_10[J_1_Q10[n]] * Z_1_Q10_10[n];
    }
    for (n in 1:N_Q11) {
      // add more terms to the linear predictor
      mu_Q11[n] += r_1_Q11_11[J_1_Q11[n]] * Z_1_Q11_11[n];
    }
    for (n in 1:N_Q12) {
      // add more terms to the linear predictor
      mu_Q12[n] += r_1_Q12_12[J_1_Q12[n]] * Z_1_Q12_12[n];
    }
    for (n in 1:N_Q01) {
      target += empty_cumulative_lpmf(Y_Q01[n] | mu_Q01[n], Intercept_Q01);
    }
    for (n in 1:N_Q02) {
      target += empty_cumulative_lpmf(Y_Q02[n] | mu_Q02[n], Intercept_Q02);
    }
    for (n in 1:N_Q03) {
      target += empty_cumulative_lpmf(Y_Q03[n] | mu_Q03[n], Intercept_Q03);
    }
    for (n in 1:N_Q04) {
      target += empty_cumulative_lpmf(Y_Q04[n] | mu_Q04[n], Intercept_Q04);
    }
    for (n in 1:N_Q05) {
      target += empty_cumulative_lpmf(Y_Q05[n] | mu_Q05[n], Intercept_Q05);
    }
    for (n in 1:N_Q06) {
      target += empty_cumulative_lpmf(Y_Q06[n] | mu_Q06[n], Intercept_Q06);
    }
    for (n in 1:N_Q07) {
      target += empty_cumulative_lpmf(Y_Q07[n] | mu_Q07[n], Intercept_Q07);
    }
    for (n in 1:N_Q08) {
      target += empty_cumulative_lpmf(Y_Q08[n] | mu_Q08[n], Intercept_Q08);
    }
    for (n in 1:N_Q09) {
      target += empty_cumulative_lpmf(Y_Q09[n] | mu_Q09[n], Intercept_Q09);
    }
    for (n in 1:N_Q10) {
      target += empty_cumulative_lpmf(Y_Q10[n] | mu_Q10[n], Intercept_Q10);
    }
    for (n in 1:N_Q11) {
      target += empty_cumulative_lpmf(Y_Q11[n] | mu_Q11[n], Intercept_Q11);
    }
    for (n in 1:N_Q12) {
      target += empty_cumulative_lpmf(Y_Q12[n] | mu_Q12[n], Intercept_Q12);
    }

         for(n in 1:N) {
              array[12] real mus = {mu_Q01[n], mu_Q02[n], mu_Q03[n], mu_Q04[n], mu_Q05[n], mu_Q06[n], mu_Q07[n], mu_Q08[n], mu_Q09[n], mu_Q10[n], mu_Q11[n], mu_Q12[n]};
              array[12] int Ys = {Y_Q01[n], Y_Q02[n], Y_Q03[n], Y_Q04[n], Y_Q05[n], Y_Q06[n], Y_Q07[n], Y_Q08[n], Y_Q09[n], Y_Q10[n], Y_Q11[n], Y_Q12[n]};

              vector[12] z;
              real prev;
              prev = 0;
              for (d in 1:12) {
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
                if (d < 12) prev = L_rescor[d+1,1:d] * head(z, d);
                // Jacobian adjustments imply z is truncated standard normal
                // thus utility --- mu + L_rescor * z --- is truncated multivariate normal
              }
          }

  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));

    // Surv model
  beta_surv ~ normal(0, 1);
  assoc ~ normal(0, 0.5);
  w ~ dirichlet(rep_vector(3.0, K));
  eta0 ~ normal(0, 0.25);


  vector[N_1] eta = beta_surv * Xsurv;

  for(i in 1:N_1){
    real lpred_q01 = eval_linpred(event_times[i], Xc_Q01[i, 2], b_Q01, r_1_Q01_1[i]);
    //real lpred_q02 = eval_linpred(event_times[i], Xc_Q02[i, 2], b_Q02, r_1_Q02_2[i]);
    //real lpred_q03 = eval_linpred(event_times[i], Xc_Q03[i, 2], b_Q03, r_1_Q03_3[i]);
    //real lpred_q04 = eval_linpred(event_times[i], Xc_Q04[i, 2], b_Q04, r_1_Q04_4[i]);
    //real lpred_q05 = eval_linpred(event_times[i], Xc_Q05[i, 2], b_Q05, r_1_Q05_5[i]);
    //real lpred_q06 = eval_linpred(event_times[i], Xc_Q06[i, 2], b_Q06, r_1_Q06_6[i]);
    //real lpred_q07 = eval_linpred(event_times[i], Xc_Q07[i, 2], b_Q07, r_1_Q07_7[i]);
    //real lpred_q08 = eval_linpred(event_times[i], Xc_Q08[i, 2], b_Q08, r_1_Q08_8[i]);
    //real lpred_q09 = eval_linpred(event_times[i], Xc_Q09[i, 2], b_Q09, r_1_Q09_9[i]);
    //real lpred_q10 = eval_linpred(event_times[i], Xc_Q10[i, 2], b_Q10, r_1_Q10_10[i]);
    //real lpred_q11 = eval_linpred(event_times[i], Xc_Q11[i, 2], b_Q11, r_1_Q11_11[i]);
    //real lpred_q12 = eval_linpred(event_times[i], Xc_Q12[i, 2], b_Q12, r_1_Q12_12[i]);

    real h0 = dot_product(theta, to_vector(B_obs[i]'));
    real H0 = jm_mvprobit_cumhazard_gk(
                    B_obs[i],
                    B_quad[i],
                    quad_times[i],
                    GK_weights,
                    theta,
                    event_times[i],
                    //1.0,
                    beta_surv,
                    Xsurv[i],
                    assoc,
                    //event[i], //surv params
                    b_Q01, r_1_Q01_1[i] ,// long params
                    b_Q02, r_1_Q02_2[i] ,
                    b_Q03, r_1_Q03_3[i] ,
                    b_Q04, r_1_Q04_4[i] ,
                    b_Q05, r_1_Q05_5[i] ,
                    b_Q06, r_1_Q06_6[i] ,
                    b_Q07, r_1_Q07_7[i] ,
                    b_Q08, r_1_Q08_8[i] ,
                    b_Q09, r_1_Q09_9[i] ,
                    b_Q10, r_1_Q10_10[i],
                    b_Q11, r_1_Q11_11[i],
                    b_Q12, r_1_Q12_12[i]
                    );

    // current value association
    real eta_tot = eta[i] + assoc[1]*lpred_q01;// + assoc[2]*lpred_q02 +
                //assoc[3]*lpred_q03 + assoc[4]*lpred_q04 + assoc[5]*lpred_q05 +
                //assoc[6]*lpred_q06 + assoc[7]*lpred_q07 + assoc[8]*lpred_q08 +
                //assoc[9]*lpred_q09 + assoc[10]*lpred_q10 + assoc[11]*lpred_q11 +
                //assoc[12]*lpred_q12;

    // Survival likelihood
    target += event[i] * (log(h0) + eta_tot) - H0;



  }

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
  // compute actual thresholds
  vector[nthres_Q05] b_Q05_Intercept = Intercept_Q05 + dot_product(means_X_Q05, b_Q05);
  // compute actual thresholds
  vector[nthres_Q06] b_Q06_Intercept = Intercept_Q06 + dot_product(means_X_Q06, b_Q06);
  // compute actual thresholds
  vector[nthres_Q07] b_Q07_Intercept = Intercept_Q07 + dot_product(means_X_Q07, b_Q07);
  // compute actual thresholds
  vector[nthres_Q08] b_Q08_Intercept = Intercept_Q08 + dot_product(means_X_Q08, b_Q08);
  // compute actual thresholds
  vector[nthres_Q09] b_Q09_Intercept = Intercept_Q09 + dot_product(means_X_Q09, b_Q09);
  // compute actual thresholds
  vector[nthres_Q10] b_Q10_Intercept = Intercept_Q10 + dot_product(means_X_Q10, b_Q10);
  // compute actual thresholds
  vector[nthres_Q11] b_Q11_Intercept = Intercept_Q11 + dot_product(means_X_Q11, b_Q11);
  // compute actual thresholds
  vector[nthres_Q12] b_Q12_Intercept = Intercept_Q12 + dot_product(means_X_Q12, b_Q12);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;

     corr_matrix[12] Rescor = multiply_lower_tri_self_transpose(L_rescor);
     vector<lower=-1,upper=1>[66] rescor;
     // extract upper diagonal of rescor matrix
     for (k in 1:12) {
        for (j in 1:(k - 1)) {
          rescor[choose(k - 1, 2) + j] = Rescor[j, k];
        }
      }

  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }

    vector[N_1] lpred_q01;
    vector[N_1] lpred_q02;
    vector[N_1] lpred_q03;
    vector[N_1] lpred_q04;
    vector[N_1] lpred_q05;
    vector[N_1] lpred_q06;
    vector[N_1] lpred_q07;
    vector[N_1] lpred_q08;
    vector[N_1] lpred_q09;
    vector[N_1] lpred_q10;
    vector[N_1] lpred_q11;
    vector[N_1] lpred_q12;

    for(i in 1:N_1){
        lpred_q01[i] = eval_linpred(event_times[i], Xc_Q01[i, 2], b_Q01, r_1_Q01_1[i]);
        lpred_q02[i] = eval_linpred(event_times[i], Xc_Q02[i, 2], b_Q02, r_1_Q02_2[i]);
        lpred_q03[i] = eval_linpred(event_times[i], Xc_Q03[i, 2], b_Q03, r_1_Q03_3[i]);
        lpred_q04[i] = eval_linpred(event_times[i], Xc_Q04[i, 2], b_Q04, r_1_Q04_4[i]);
        lpred_q05[i] = eval_linpred(event_times[i], Xc_Q05[i, 2], b_Q05, r_1_Q05_5[i]);
        lpred_q06[i] = eval_linpred(event_times[i], Xc_Q06[i, 2], b_Q06, r_1_Q06_6[i]);
        lpred_q07[i] = eval_linpred(event_times[i], Xc_Q07[i, 2], b_Q07, r_1_Q07_7[i]);
        lpred_q08[i] = eval_linpred(event_times[i], Xc_Q08[i, 2], b_Q08, r_1_Q08_8[i]);
        lpred_q09[i] = eval_linpred(event_times[i], Xc_Q09[i, 2], b_Q09, r_1_Q09_9[i]);
        lpred_q10[i] = eval_linpred(event_times[i], Xc_Q10[i, 2], b_Q10, r_1_Q10_10[i]);
        lpred_q11[i] = eval_linpred(event_times[i], Xc_Q11[i, 2], b_Q11, r_1_Q11_11[i]);
        lpred_q12[i] = eval_linpred(event_times[i], Xc_Q12[i, 2], b_Q12, r_1_Q12_12[i]);
    }



}

