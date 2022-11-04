functions {
#include strapp_funs.stan
  // vector strapp_cure(vector beta0, real int0, matrix X0, matrix inv_sqrt_X0tX0, real pcured) {
  //   int n0 = rows(X0);
  //   vector[n0] W = inv_logit(int0 + X0 * beta0);
  //   W = W .* (1 - W);
  //   return (1 - pcured) * inv_sqrt_X0tX0 * matrix_sqrt(X0' * (diag_pre_multiply(W, X0)) * beta0);
  // }
  // vector strapp_ph(vector beta0, real int0, matrix X0, matrix inv_sqrt_X0tX0) {
  //   int n0 = rows(X0);
  //   vector[n0] W = inv_logit(int0 + X0 * beta0);
  //   W = W .* (1 - W);
  //   return inv_sqrt_X0tX0 * matrix_sqrt(X0' * (diag_pre_multiply(W, X0)) * beta0);
  // }
}
data {
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    p;             // number of covariates (EXCLUDING intercept term)
  int<lower=0,upper=1>            y0[n0];        // binary historical data
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
}
transformed data {
  matrix[p,p] inv_sqrt_X0tX0 = inverse_spd( matrix_sqrt( crossprod(X0) ) );
}
// The parameter accepted by the model
parameters {
  real                  intercept_hist;   // intercept for historical data
  vector[p]             beta_hist;        // reg. coefs for logistic reg model (includes intercept)
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // power prior: logistic regression for historical data
  target += bernoulli_logit_glm_lpmf(y0 | X0, intercept_hist, beta_hist);
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta_hist | 0, 10);
}
// 
// generated quantities {
//   real pcure = uniform_rng(0, 1);
//   vector[p] beta_cur_cure;
//   beta_cur_cure = strapp_cure(beta_hist, intercept_hist, X0, inv_sqrt_X0tX0, pcure);
// }

