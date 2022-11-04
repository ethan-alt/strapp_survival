
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n1;              // current data sample size
  int<lower=0>                    n0;              // historical data sample size
  int<lower=0>                    J;               // number of time intervals
  int<lower=0>                    p;               // number of covariates
  vector[n1]                      y1;              // time to event current data
  vector[n0]                      y0;              // normal historical data
  matrix[n1,p]                    X1;              // current data design matrix (no intercept for current data likelihood)
  matrix[n0,p]                    X0;              // historical data design matrix (EXCLUDING intercept term)
  int<lower=1,upper=J>            intindx[n1];     // index giving interval into which obs i failed / was censored
  int<lower=0,upper=1>            death_ind[n1];   // event indicator (1 = event; 0 = censored)
  vector[J+1]                     breaks;          // J+1-dim interval of breaks
  real<lower=0,upper=1>           a0;              // power prior parameter
  vector<lower=0>[J]              hazard_shape;    // shape parameter for gamma prior on hazards
  vector<lower=0>[J]              hazard_rate;     // rate (inverse scale) parameter for gamma prior on hazards
  real<lower=0>                   precision_shape; // shape hyperparameter for gamma prior on precision for historical data
  real<lower=0>                   precision_rate;  // rate hyperparameter for gamma prior on precision for historical data
  
}

transformed data {
  matrix[p,p] inv_sqrt_X0tX0 = inverse_spd( matrix_sqrt( crossprod(X0) ) );
}

// The parameter accepted by the model
parameters {
  vector<lower=0>[J]    lambda;          // the J hazards for each interval
  real                  intercept_hist;  // historical intercept
  vector[p]             beta_hist;       // reg. coefs for historical data
  vector[p]             c0;              // gen strapp parameter
  real<lower=0>         c0_sd;           // gen strapp SD
  real<lower=0>         precision;       // precision parameter for historical data
}

transformed parameters {
  // SD for historical data
  real sigma = inv_sqrt(precision);
  
  // regression coefficients
  vector[p] beta = inv(sigma) * beta_hist - inv_sqrt_X0tX0 * c0;
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
  target += gamma_lpdf(lambda | hazard_shape, hazard_rate);
  
  // gamma prior on precision for historical data
  target += gamma_lpdf(precision | precision_shape, precision_rate);
  
  // Normal prior on c0
  target += normal_lpdf(c0 | 0, c0_sd);
  
  // half std normal prior on c0_sd
  target += std_normal_lpdf(c0_sd) + 0.69314718056;  // add log(2)
  
  // power prior: linear regression for historical data
  target += a0 * normal_id_glm_lpdf(y0 | X0, intercept_hist, beta_hist, sigma);
  
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta_hist | 0, 10);
  
  eta = X1 * beta;
  for ( i in 1:n1 ) {
    target += pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz);
  }
}


