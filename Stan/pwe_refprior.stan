
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of covariates
  vector[n1]                      y1;            // time to event current data
  matrix[n1,p]                    X1;            // current data design matrix
  int<lower=1,upper=J>            intindx[n1];   // index giving interval into which obs i failed / was censored
  int<lower=0,upper=1>            death_ind[n1]; // event indicator (1 = event; 0 = censored)
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  vector<lower=0>[J]              hazard_shape;  // shape parameter for gamma prior on hazards
  vector<lower=0>[J]              hazard_rate;   // rate (inverse scale) parameter for gamma prior on hazards
}

transformed data {
  real logcens[n1];
  matrix[p,p] sqrt_crossprod_X0_inv;
  
  // compute censoring indicator
  for ( i in 1:n1 )
    logcens[i] = log(1 - death_ind[i]);
}

// The parameter accepted by the model
parameters {
  vector<lower=0>[J]    lambda;       // the J hazards for each interval
  vector[p]             beta;
}

model {
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz;
  vector[n1]  eta;
  
  // cumulative basline hazard
  cumblhaz = cumulative_sum(lambda[1:(J-1)] .* (breaks[2:J] - breaks[1:(J-1) ] ) );
    
  // gamma prior on hazards
  target += gamma_lpdf(lambda | hazard_shape, hazard_rate);
  
  // noninformative normal prior on regression coefficients
  target += normal_lpdf(beta | 0, 10);
  
  // Current data likelihood:
  eta = X1 * beta;
  for ( i in 1:n1 )
    target += pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz);
}


