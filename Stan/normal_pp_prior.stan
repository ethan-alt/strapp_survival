functions {
#include strapp_funs.stan
  // vector strapp_ph(vector beta0, real int0, matrix X0, matrix inv_sqrt_X0tX0, real sigma) {
  //   int n0 = rows(X0);
  //   return inv_sqrt_X0tX0 * matrix_sqrt(X0' * (diag_pre_multiply(W, X0)) * beta0);
  // }
  // vector strapp_cure(vector beta0, real int0, matrix X0, matrix inv_sqrt_X0tX0, real sigma, real pcured) {
  //   int n0 = rows(X0);
  //   vector[n0] W = inv_logit(int0 + X0 * beta0);
  //   W = W .* (1 - W);
  //   return (1 - pcured) * inv_sqrt_X0tX0 * matrix_sqrt(X0' * (diag_pre_multiply(W, X0)) * beta0);
  // }
}
data {
  int<lower=0>                    n0;              // historical data sample size
  int<lower=0>                    p;               // number of covariates (EXCLUDING intercept term)
  vector[n0]                      y0;
  matrix[n0,p]                    X0;              // historical data design matrix (EXCLUDING intercept term)
  real<lower=0>                   precision_shape; // shape hyperparameter for gamma prior on precision
  real<lower=0>                   precision_rate;  // rate hyperparameter for gamma prior on precision
}

// The parameter accepted by the model
parameters {
  real                  intercept_hist;   // intercept for historical data
  vector[p]             beta_hist;        // reg. coefs for logistic reg model (includes intercept)
  real<lower=0>         precision;
}

transformed parameters {
  real sigma = inv_sqrt(precision);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // power prior: normal lm regression for historical data
  target += normal_id_glm_lpdf(y0 | X0, intercept_hist, beta_hist, sigma);
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta_hist | 0, 10);
  target += gamma_lpdf(precision | precision_shape, precision_rate);
}
generated quantities {
  real sigmasq = square(sigma);
  // real pcure = uniform_rng(0, 1);
  // vector[p] beta_cur_ph;
  // vector[p] beta_cur_cure;
  // beta_cur_ph   = sqrt(precision) * beta_hist;
  // beta_cur_cure = sqrt(precision) * inv_sqrt(1-pcure) * beta_hist;
}

