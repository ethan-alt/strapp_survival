
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of covariates (including intercept term)
  vector[n1]                      y1;            // time to event current data
  vector[n0]                      y0;            // normal historical data
  matrix[n1,p]                    X1;            // current data design matrix (no intercept for current data likelihood)
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
  int<lower=1,upper=J>            intindx[n1];   // index giving interval into which obs i failed / was censored
  int<lower=0,upper=1>            death_ind[n1]; // event indicator (1 = event; 0 = censored)
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  real<lower=0,upper=1>           a0;            // power prior parameter
  vector<lower=0>[J]              hazard_shape;  // shape parameter for gamma prior on hazards
  vector<lower=0>[J]              hazard_rate;   // rate (inverse scale) parameter for gamma prior on hazards
  real<lower=0>                   precision_shape;
  real<lower=0>                   precision_rate;
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
  real                  intercept_hist;   // intercept for historical data
  real<lower=0,upper=1> p_cured;      // proportion of cured individuals
  vector<lower=0>[J]    lambda;       // the J hazards for each interval
  vector[p]             beta_hist;    // reg. coefs for logistic reg model (includes intercept)
  real<lower=0>         precision;    // precision parameter for normal 
  vector[p]             c0;           // genstrapp parameter
  real<lower=0>         c0_sd;        // SD for genstrapp
}

// current beta is a function of historical beta
transformed parameters {
  real sigma = inv_sqrt(precision);
  // current beta = sqrt(cov(beta)) * inv(sqrt(cov(beta0[2:p]))) * beta0[2:p]
  vector[p] beta;
  beta = inv_sqrt(1 - p_cured) * ( inv(sigma) * beta_hist - invsqrt_fisher_pwe * c0) ;
}

model {
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz;
  vector[n1]  eta;
  
  // cumulative basline hazard
  cumblhaz = cumulative_sum(lambda[1:(J-1)] .* (breaks[2:J] - breaks[1:(J-1) ] ) );
    
  // gamma prior on hazards
  target += gamma_lpdf(lambda | hazard_shape, hazard_rate);
  
  // gamma prior on precision
  target += gamma_lpdf(precision | precision_shape, precision_rate);
  
  // prior for c0 is normal
  target += normal_lpdf(c0 | 0, c0_sd);
  
  // half std normal prior on c0_sd
  target += std_normal_lpdf(c0_sd) + 0.69314718056;  // add on log(2)
  
  // power prior: logistic regression for historical data
  target += a0 * normal_id_glm_lpdf(y0 | X0, intercept_hist, beta_hist, sigma);
  
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta_hist | 0, 10);
  
  // Current data likelihood:
  eta = X1 * beta;
  for ( i in 1:n1 ) {
    target += log_mix(
      p_cured,
      logcens[i],
      pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz)
    );
  }
}







