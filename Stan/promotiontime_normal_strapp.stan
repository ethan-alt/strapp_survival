
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of covariates (including intercept)
  vector[n1]                      y1;            // time to event current data
  vector[n0]                      y0;            // normal historical data
  matrix[n1,p]                    X1;            // current data design matrix (EXCLUDING intercept)
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
  int<lower=1,upper=J>            intindx[n1];   // index giving interval into which obs i failed / was censored
  int<lower=0,upper=1>            death_ind[n1]; // event indicator (1 = event; 0 = censored)
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  real<lower=0,upper=1>           a0;            // power prior parameter
  vector<lower=0>[J]              hazard_shape;  // shape parameter for gamma prior on hazards
  vector<lower=0>[J]              hazard_rate;   // rate (inverse scale) parameter for gamma prior on hazards
  real<lower=0>                   precision_shape; // shape hyperparameter for gamma prior on precision for historical data
  real<lower=0>                   precision_rate;  // rate hyperparameter for gamma prior on precision for historical data
}
transformed data {
  matrix[p,p] invsqrt_X0tX0 = inverse_spd( matrix_sqrt( crossprod( X0 ) ) );
  real log_det_invsqrt_XtX = log_determinant(invsqrt_X0tX0);
}
// The parameter accepted by the model
parameters {
  vector<lower=0>[J]    lambda;       // the J hazards for each interval
  real                  intercept_hist;
  real                  intercept;
  vector[p]             beta;         // reg. coefs
  real<lower=0>         precision;    // precision parameter for historical data
}

// current beta is a function of historical beta
transformed parameters {
  // SD for historical data
  real sigma = inv_sqrt(precision);
  vector[p] beta_hist = sigma * invsqrt_X0tX0 * ( sqrt_fisher_promotiontime(X0, beta, intercept ) * beta ); 
}

// The model to be estimated
model {
  vector[J]   log_lambda = log(lambda);
  vector[J-1] lambda_d_breaks;
  vector[n1]  eta;
  
  // compute lambda[j] * (s[j] * s[j-1])
  for ( j in 2:J )
    lambda_d_breaks[j-1] = lambda[j-1] * (breaks[j] - breaks[j-1]);
    
  // gamma prior on hazards
  target += gamma_lpdf(lambda | hazard_shape, hazard_rate);
  
  // gamma prior on precision for historical data
  target += gamma_lpdf(precision | precision_shape, precision_rate);
  
  // power prior: linear regression for historical data
  target += a0 * normal_id_glm_lpdf(y0 | X0, intercept_hist, beta_hist, sigma);
  
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta_hist | 0, 10);
  target += normal_lpdf(intercept | 0, 10);

  // likelihood for current data
  eta = intercept + X1 * beta;
  for ( i in 1:n1 ) {
    target += promotion_time_cure_rate_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], lambda_d_breaks);
  }
  
  // Jacobian
  target += log_abs_det_jacob_promotion(beta, intercept, sigma, X0, log_det_invsqrt_XtX );
}


