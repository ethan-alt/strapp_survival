
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of covariates (including intercept term)
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
}

transformed data {
  real logcens[n1];
  matrix[p,p] invsqrt_fisher_pwe = matrix_sqrt( inverse_spd( crossprod(X0) ) );
  
  // compute censoring indicator
  for ( i in 1:n1 )
    logcens[i] = log(1 - death_ind[i]);
}

// The parameter accepted by the model
parameters {
  vector<lower=0>[J]    lambda;       // the J hazards for each interval
  real                  intercept_hist;   // intercept for historical data
  vector[p]             beta_hist;    // reg. coefs for logistic reg model (includes intercept)
}

// current beta is a function of historical beta
transformed parameters {
  vector[p] beta;
  beta =   invsqrt_fisher_pwe 
           * ( sqrt_fisher_logistic(X0, beta_hist, intercept_hist) * beta_hist );
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz;
  vector[n1]  eta;
  
  // compute lambda[j] * (s[j] * s[j-1])
  cumblhaz = cumulative_sum(lambda[1:(J-1)] .* (breaks[2:J] - breaks[1:(J-1) ] ) );
    
  // gamma prior on hazards
  target += gamma_lpdf(lambda | hazard_shape, hazard_shape);
  
  // power prior: logistic regression for historical data
  target += a0 * bernoulli_logit_glm_lpmf(y0 | X0, intercept_hist, beta_hist);
  
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta_hist | 0, 10);
  
  // Current data likelihood:
  eta = X1 * beta;
  for ( i in 1:n1 )
    target += pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz);
}







