
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
}

transformed data {
  real logcens[n1];
  
  // compute censoring indicator
  for ( i in 1:n1 )
    logcens[i] = log(1 - death_ind[i]);
}

// The parameter accepted by the model
parameters {
  vector<lower=0>[J]    lambda;         // the J hazards for each interval
  vector[p] beta;                       // regression coefficients
  real intercept_hist;
  real intercept;
}

// The model to be estimated
model {
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz;
  vector[n1]  eta;
  cumblhaz = cumulative_sum(lambda[1:(J-1)] .* (breaks[2:J] - breaks[1:(J-1) ] ) );
    
  // gamma prior on hazards
  target += gamma_lpdf(lambda | hazard_shape, hazard_rate);
  
  // power prior: logistic regression for historical data
  target += a0 * bernoulli_logit_glm_lpmf(y0 | X0, intercept_hist, beta);
  
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta | 0, 10);
  target += normal_lpdf(intercept | 0, 10);
  
  // Current data likelihood:
  //  log_mix(p, l1, l2) = log( exp( log(p) + l1 ) + exp( log(p) + l2 ) )
  //                     = log( p * exp(l1) + (1 - p) * exp(l2) )
  //     p  = mixture prob
  //     l1 = log likelihood for first mixture  = log(cens[i]) = { log(-Inf) if event; 0 if censored }
  //     l2 = log likelihood for second mixture = log PWE likelihood
  eta = intercept + X1 * beta;
  for ( i in 1:n1 ) {
    target += promotion_time_cure_rate_lpdf(
      y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz
    );
  }
}


