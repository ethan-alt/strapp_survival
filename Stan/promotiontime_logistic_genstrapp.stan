
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of covariates
  vector[n1]                      y1;            // time to event current data
  int<lower=0,upper=1>            y0[n0];        // binary historical data
  matrix[n1,p]                    X1;            // current data design matrix (no intercept for current data likelihood)
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
  int<lower=1,upper=J>            intindx[n1];   // index giving interval into which obs i failed / was censored
  int<lower=0,upper=1>            death_ind[n1]; // event indicator (1 = event; 0 = censored)
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  real<lower=0,upper=1>           a0;            // power prior parameter
  vector<lower=0>[J]              hazard_shape;  // shape parameter for gamma prior on hazards
  vector<lower=0>[J]              hazard_rate;   // rate (inverse scale) parameter for gamma prior on hazards
  real<lower=0>                   rel_tol;
  real<lower=0>                   f_tol;
  int<lower=0>                    max_steps;
}

transformed data {
  // For nonlinear equation solver
  int  x_i[2]    = { n0, p };         // integer array giving rows and columns of matrix
  real x_r[n0*p] = to_array_1d(X0);   // convert matrix to 1-dim array
  
  
  real logcens[n1];
  matrix[p,p] sqrt_crossprod_X0_inv;
  
  // compute censoring indicator
  for ( i in 1:n1 )
    logcens[i] = log(1 - death_ind[i]);
  
  // compute inverse of the square root of X0'X0, i.e., (X0'X0)^(-1/2)
  sqrt_crossprod_X0_inv = inverse_spd( matrix_sqrt( crossprod(X0) ) );
}

// The parameter accepted by the model
parameters {
  vector<lower=0>[J]    lambda;         // the J hazards for each interval
  real                  intercept_hist; // intercept for historical data set
  vector[p]             beta_hist;      // reg. coefs for current data model
  vector[p]             c0;             // gen strapp location parameter
  real<lower=0>         c0_sd;          // gen strapp SD of location parameter
  real                  intercept;      // intercept for current data
}

// current beta is a function of historical beta
transformed parameters {
  vector[p] beta = algebra_solver(
    system_genstrapp, beta_hist, append_row(beta_hist, append_row(intercept_hist, append_row(c0, intercept))), x_r, x_i
    , rel_tol, f_tol, max_steps
  );
}

// The model to be estimated
model {
  vector[J]   log_lambda = log(lambda);
  // vector[J-1] lambda_d_breaks;
  vector[J-1] cumblhaz;
  vector[n1]  eta;
  
  // // compute lambda[j] * (s[j] * s[j-1])
  // for ( j in 2:J )
  //   lambda_d_breaks[j-1] = lambda[j-1] * (breaks[j] - breaks[j-1]);
  cumblhaz = cumulative_sum(lambda[1:(J-1)] .* (breaks[2:J] - breaks[1:(J-1) ] ) );
    
  // gamma prior on hazards
  target += gamma_lpdf(lambda | hazard_shape, hazard_rate);
  
  // normal prior on c0
  target += normal_lpdf(c0 | 0, c0_sd);
  
  // half std normal prior on c0_sd
  target += std_normal_lpdf(c0_sd) + 0.69314718056;  // add on log(2)
  
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta_hist | 0, 10);
  target += normal_lpdf(intercept | 0, 10);
  
  // power prior: logistic regression for historical data
  target += a0 * bernoulli_logit_glm_lpmf(y0 | X0, intercept_hist, beta_hist);
  
  // Current data likelihood:
  eta = intercept + X1 * beta;
  for ( i in 1:n1 ) {
    target += promotion_time_cure_rate_lpdf(
      y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz
    );
  }
}


