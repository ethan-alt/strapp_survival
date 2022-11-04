library(cmdstanr)
library(posterior)
library(bayesplot)

## directory where stan files are located
standir                         <- '/proj/ibrahimlab/strapp_survival/Stan'

samplers.list = list(
  ##
  ## Cure rate models with PWE
  ##
  curepwe_logistic_strapp_stan    = cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_strapp.stan'), include_paths = standir)
  , curepwe_logistic_genstrapp_stan = cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_genstrapp.stan'), include_paths = standir)
  , curepwe_logistic_pp_stan        = cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_pp.stan'), include_paths = standir)
  , curepwe_normal_strapp_stan      = cmdstanr::cmdstan_model(file.path(standir, 'curepwe_normal_strapp.stan'), include_paths = standir)
  , curepwe_normal_genstrapp_stan   = cmdstanr::cmdstan_model(file.path(standir, 'curepwe_normal_genstrapp.stan'), include_paths = standir)
  , curepwe_normal_pp_stan          = cmdstanr::cmdstan_model(file.path(standir, 'curepwe_normal_pp.stan'), include_paths = standir)
  , curepwe_refprior_stan           = cmdstanr::cmdstan_model(file.path(standir, 'curepwe_refprior.stan'), include_paths = standir)
  ##
  ## PWE models
  ##
  , pwe_logistic_strapp_stan    = cmdstanr::cmdstan_model(file.path(standir, 'pwe_logistic_strapp.stan'), include_paths = standir)
  , pwe_logistic_genstrapp_stan = cmdstanr::cmdstan_model(file.path(standir, 'pwe_logistic_genstrapp.stan'), include_paths = standir)
  , pwe_logistic_pp_stan        = cmdstanr::cmdstan_model(file.path(standir, 'pwe_logistic_pp.stan'), include_paths = standir)
  , pwe_normal_strapp_stan      = cmdstanr::cmdstan_model(file.path(standir, 'pwe_normal_strapp.stan'), include_paths = standir)
  , pwe_normal_genstrapp_stan   = cmdstanr::cmdstan_model(file.path(standir, 'pwe_normal_genstrapp.stan'), include_paths = standir)
  , pwe_normal_pp_stan          = cmdstanr::cmdstan_model(file.path(standir, 'pwe_normal_pp.stan'), include_paths = standir)
  , pwe_refprior_stan           = cmdstanr::cmdstan_model(file.path(standir, 'pwe_refprior.stan'), include_paths = standir)
  ##
  ## Promotion time cure rate models
  ##
  , promotiontime_logistic_strapp_stan    = cmdstanr::cmdstan_model(file.path(standir, 'promotiontime_logistic_strapp.stan'), include_paths = standir)
  , promotiontime_logistic_genstrapp_stan = cmdstanr::cmdstan_model(file.path(standir, 'promotiontime_logistic_genstrapp.stan'), include_paths = standir)
  , promotiontime_logistic_pp_stan        = cmdstanr::cmdstan_model(file.path(standir, 'promotiontime_logistic_pp.stan'), include_paths = standir)
  , promotiontime_normal_strapp_stan      = cmdstanr::cmdstan_model(file.path(standir, 'promotiontime_normal_strapp.stan'), include_paths = standir)
  , promotiontime_normal_genstrapp_stan   = cmdstanr::cmdstan_model(file.path(standir, 'promotiontime_normal_genstrapp.stan'), include_paths = standir)
  , promotiontime_normal_pp_stan          = cmdstanr::cmdstan_model(file.path(standir, 'promotiontime_normal_pp.stan'), include_paths = standir)
  , promotiontime_refprior_stan           = cmdstanr::cmdstan_model(file.path(standir, 'promotiontime_refprior.stan'), include_paths = standir)
)

## Function to get stan data--applicable to all models
get_standata <- function(
  fmla.cur, fmla.hist, data, histdata
  , breaks = NULL, a0 = 0.25
  , hazard.shape = 0.1, hazard.rate = 0.1
  , precision.shape = 0.1, precision.rate = 0.1
  , rel.tol = 1e-6, f.tol = 1e-6, max.steps = 1000
) {
  ## Extract names and variables for response, censoring, etc.
  y1name    <- all.vars(fmla.cur)[1]
  eventname <- all.vars(fmla.cur)[2]
  y1        <- data[, y1name]
  eventind  <- data[, eventname]
  X1        <- model.matrix(fmla.cur, data)
  
  if ( is.null(histdata) ) {
    y0name <- ''
    y0 <- 0
    X0 <- matrix(0, 2, 2)
  } else {
    y0name <- all.vars(fmla.hist)[1]
    X0     <- model.matrix(fmla.hist, histdata)
    y0     <- histdata[, y0name]
  }
  ## Make sure no design matrices have intercepts
  if ( '(Intercept)' %in% colnames(X1) )
    X1 <- X1[, -1]
  if ( '(Intercept)' %in% colnames(X0) )
    X0 <- X0[, -1]
  
  ## Shape and rate for hazard as vector
  if ( !is.null(breaks) ) {
    J <- length(breaks) - 1  ## number of intervals
    if ( length(hazard.shape) == 1 )
      hazard.shape <- rep(hazard.shape, J)
    if ( length(hazard.rate) == 1 )
      hazard.rate <- rep(hazard.rate, J)
  } else {
    J <- NULL
  }
  
  ## Create index giving interval into which obs failed / was censored
  intindx <- rep(NA, nrow(X1))
  for ( j in 1:J ) {
    intindx[ y1 >= breaks[j] & ( y1 <= breaks[j+1] )  ] <- j
  }
  
  list(
      'n1'              = nrow(X1)
    , 'n0'              = nrow(X0)
    , 'J'               = J
    , 'p'               = ncol(X0)
    , 'y1'              = y1
    , 'y0'              = y0
    , 'X1'              = X1
    , 'X0'              = X0
    , 'intindx'         = intindx
    , 'death_ind'       = eventind
    , 'breaks'          = breaks
    , 'a0'              = a0
    , 'hazard_shape'    = hazard.shape
    , 'hazard_rate'     = hazard.rate
    , 'precision_rate'  = precision.rate
    , 'precision_shape' = precision.shape
    , 'rel_tol'         = rel.tol
    , 'f_tol'           = f.tol
    , 'max_steps'       = max.steps
  )
}




#' Sample from a survival model
#' Obtain posterior samples from a survival model under various priors and models
#' 
#' @param fmla.cur formula for current (time to event) data
#' @param fmla.hist formula for historical data
#' @param data a data frame for current (time to event) data
#' @param model.cur one of 'pwe' (piecewise exponential PH model), 'curepwe' (piecewise exponential PH model with cure fraction), 'promotiontime' (promotion time cure rate model)
#' @param histdata a data frame for historical (glm) data
#' @param histdata.type one of 'normal', 'logistic'
#' @param prior one of 'strapp', 'genstrapp', 'pp' (power prior), 'refprior'
#' @param breaks J+1 vector of endpoints for PWE intervals
#' @param hazard.shape shape parameter for gamma prior on hazards
#' @param hazard.rate rate parameter for gamma prior on hazards
#' @param precision.shape shape parameter for gamma prior on precision. Ignored if histdata.type == 'logistic'
#' @param precision.rate rate parameter for gamma prior on precision. Ignored if histdata.type == 'logistic'
#' @param rel.tol relative tolerance for nonlinear equation solver. Only used for logistic promotion time models.
#' @param f.tol function tolerance for nonlinear equation solver. Only used for logistic promotion time models.
#' @param max.steps maximum number of steps for nonlinear equation solver. Only used for logistic promotion time models.
#' @param ... other parameters to pass onto cmdstanr::sample
#' 
#' @return object of type cmdstanrfit
surv.sample <- function(
  fmla.cur, fmla.hist, data, histdata = NULL, model.cur, histdata.type = NULL, prior
  , breaks = NULL, a0 = 0.25, 
  hazard.shape = 0.1, hazard.rate = 0.1
  , precision.shape = 0.1, precision.rate = 0.1
  , rel.tol = 1e-6, f.tol = 1e-6, max.steps = 1000
  , ...
) {
  if ( is.null(breaks) & model.cur != 'promotiontime' )
    stop('Must specify breaks when model is pwe or curepwe')
  
  ## Get stan data
  standat <- get_standata(
    fmla.cur, fmla.hist, data, histdata, breaks
    , a0, hazard.shape, hazard.rate
    , precision.shape, precision.rate
    , rel.tol, f.tol, max.steps
  )
  
  ## get prior name based on prior and historical data type
  prior.name <- 'refprior'
  if ( prior != 'refprior' )
    prior.name <- paste(histdata.type, prior, sep = '_')
  
  ## Obtain proper sampler based on model.cur, prior.name
  sampler.name <- paste(model.cur, prior.name, 'stan', sep = '_')
  print(sampler.name)
  
  ## Obtain samples
  smpl <- samplers.list[[sampler.name]]$sample(data = standat, ...)
  attr(smpl, 'standata') <- standat
  attr(smpl, 'prior') <- prior
  attr(smpl, 'sampler.name') <- sampler.name
  smpl
}




pwe.loglik <- function(parms, standat) {
  y      <- standat$y1
  X      <- standat$X1
  fail   <- standat$death_ind
  j      <- standat$intindx
  J      <- standat$J
  p      <- standat$p
  n      <- standat$n1
  breaks <- standat$breaks
  lambda.vars <- paste0('lambda[', 1:J, ']')
  beta.vars   <- paste0('beta[', 1:p, ']')
  suppressWarnings({
    lambda      <- parms %>% as_draws_df %>% select(all_of(lambda.vars)) %>% as.numeric
    beta        <- parms %>% as_draws_df %>% select(all_of(beta.vars)) %>% as.numeric
  })
  
  ## Compute linear predictor
  eta <- (X %*% beta)[, 1]
  
  ## Get hazard corresponding to failure / censoring time
  lambda_j <- lambda[j]
  
  ## Compute cumulative baseline hazard at each interval
  cumblhaz <- cumsum( lambda[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) )
  cumblhaz <- c(0, cumblhaz)
  
  ## Compute cumulative hazard for each observation
  cumhaz <- lambda_j * (y - breaks[j]) + cumblhaz[j]
  cumhaz <- cumhaz * exp(eta)
  
  ## Log likelihood = sum( event_ind * log(hazard) - cumhaz )
  ##   log(hazard) = log( lambda * exp(eta) ) = log(lambda) + eta
  loglik <- sum( fail * ( log(lambda_j) + eta ) - cumhaz )
}


curepwe.loglik <- function(parms, standat) {
  y      <- standat$y1
  X      <- standat$X1
  fail   <- standat$death_ind
  j      <- standat$intindx
  J      <- standat$J
  p      <- standat$p
  n      <- standat$n1
  breaks <- standat$breaks
  lambda.vars <- paste0('lambda[', 1:J, ']')
  beta.vars   <- paste0('beta[', 1:p, ']')
  pcured.var  <- 'p_cured'
  suppressWarnings({
    lambda      <- parms %>% as_draws_df %>% select(all_of(lambda.vars)) %>% as.numeric
    beta        <- parms %>% as_draws_df %>% select(all_of(beta.vars)) %>% as.numeric
    pcured      <- parms %>% as_draws_df %>% select(all_of(pcured.var)) %>% as.numeric
  })
  
  ## Compute linear predictor
  eta <- (X %*% beta)[, 1]
  
  ## Get hazard corresponding to failure / censoring time
  lambda_j <- lambda[j]
  
  ## Compute cumulative baseline hazard at each interval
  cumblhaz <- cumsum( lambda[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) )
  cumblhaz <- c(0, cumblhaz)
  
  ## Compute cumulative hazard for each observation
  cumhaz <- lambda_j * (y - breaks[j]) + cumblhaz[j]
  cumhaz <- cumhaz * exp(eta)
  
  ## likelihood[i] = pcured * (1 - fail[i]) + (1 - pcured) * h(y[i])^fail[i] * S(y[i])
  loglik <- sum( log( pcured * (1 - fail) + (1 - pcured) * lambda_j^fail * exp(-cumhaz) ) )
}


promotiontime.loglik <- function(parms, standat) {
  y      <- standat$y1
  X      <- standat$X1
  fail   <- standat$death_ind
  j      <- standat$intindx
  J      <- standat$J
  p      <- standat$p
  n      <- standat$n1
  breaks <- standat$breaks
  lambda.vars <- paste0('lambda[', 1:J, ']')
  beta.vars   <- paste0('beta[', 1:p, ']')
  suppressWarnings({
    lambda      <- parms %>% as_draws_df %>% select(all_of(lambda.vars)) %>% as.numeric
    beta        <- parms %>% as_draws_df %>% select(all_of(beta.vars)) %>% as.numeric
  })
  
  ## Compute linear predictor
  eta <- (X %*% beta)[, 1]
  
  ## Compute cumulative baseline hazard at each interval
  cumblhaz <- cumsum( lambda[1:(J-1)] * ( breaks[2:J] - breaks[1:(J-1)] ) )
  cumblhaz <- c(0, cumblhaz)
  
  ## Compute cumulative hazard for each observation
  cumhaz <- lambda[j] * (y - breaks[j]) + cumblhaz[j]
  
  ## log likelihood[i] = fail[i] * ( eta[i] + loghazard[i] + cumhaz[i] ) - exp(eta[i]) * (1 - S(y[i]))
  ##   S(y[i]) = exp(-cumhaz[i])
  sum( 
    fail * (eta + log(lambda[j]) + cumhaz[j]) - exp(eta) * (1 - exp(-cumhaz))
  )
}






