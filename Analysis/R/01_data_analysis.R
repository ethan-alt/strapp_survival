remove(list = ls())

library(dplyr)
library(survival)
library(posterior)

## Set cmdstanr directory (necessary for sending code to cluster)
library(cmdstanr)
cmdstanr::set_cmdstan_path("/nas/longleaf/home/ethanalt/.cmdstanr/cmdstan-2.27.0")

## Directory to save results
save.dir <- '/proj/ibrahimlab/strapp_survival/Analysis/Results'

## source wrappers
wrapper.dir <- '/proj/ibrahimlab/strapp_survival/R'
source(file.path(wrapper.dir, 'wrappers_logistic.R'))

## Sampling parameters
nburnin  = 2000
nsamples = 25000
nchains  = 1


## Obtain scenarios for computation--remove any reference priors with a0s
a0        <- c(0.25, 0.50, 0.75, 1.00)
priors    <- c('refprior', 'pp', 'strapp', 'genstrapp')
models    <- c('pwe', 'curepwe', 'promotiontime')
hist.type <- c('normal', 'logistic')

scenarios <- expand.grid(
  a0 = a0, priors = priors, models = models, hist.type = hist.type
  , stringsAsFactors = FALSE
)
scenarios <- scenarios %>% filter(!(priors == 'refprior' & a0 > 0.25))
scenarios[scenarios$priors == 'refprior' & scenarios$a0 == 0.25, 'a0'] <- 0

## Remove promotion time for logistic model
scenarios <- scenarios %>% filter(!(hist.type == 'logistic' & models == 'promotiontime'))

## Get scenario for current ID
id <- as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if ( is.na(id) )
  id <- 1

## Obtain file name based on id
filename <- file.path(save.dir, paste0('analysis_', id, '.rds'))

scen <- scenarios[id, ]
a0.id    <- scen$a0
prior.id <- scen$priors
model.id <- scen$models
hist.id  <- scen$hist.type

## load data
data.dir   <- '/proj/ibrahimlab/strapp_survival/Data'
curdata  <- readRDS(file.path(data.dir, 'data_cur_1690.rds'))
if ( hist.id == 'logistic' ) {
  histdata  <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
  curdata   <- readRDS(file.path(data.dir, 'data_cur_1690.rds'))
  fmla.hist <- fail_2yr                ~ trt + age_c + sex + perform
  fmla.cur  <- Surv(failtime, rfscens) ~ trt + age_c + sex + perform
  curdata <- curdata %>% filter(failtime > 0)
} else if (hist.id == 'normal') {
  ## Get data
  histdata    <- readRDS(file.path(data.dir, 'data_hist_2696.rds'))
  curdata     <- readRDS(file.path(data.dir, 'data_cur_1694.rds'))
  ## Get formula
  fmla.hist   <- log_igm28               ~ ifn_gmk + age_c + sex + perform
  fmla.cur    <- Surv(failtime, rfscens) ~ trt     + age_c + sex + perform
  ## Fix coding of current data treatment
  curdata$trt <- 1 - curdata$gmk
  ## Replace 0 times with 0.50 (half day)
  curdata <- curdata %>% mutate(failtime = if_else(failtime == 0, 0.50, failtime))
}
  
## Get sample sizes
n1 <- nrow(curdata)
n0 <- nrow(histdata)


# ## Breaks for TTE data
# breaks = c(0, 4.02465, 5.83165, 12.3039, 19.94255, 72)
# nbreaks = length(breaks)
# 
nbreaks = 5
probs   = 1:nbreaks / nbreaks
breaks  = curdata %>%
  filter(rfscens == 1) %>%
  summarize(quant = quantile(failtime, probs = probs)) %>%
  unlist
breaks = c(0, breaks)
breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])



## Center and scale age and nodes
curdata$age_c <- as.numeric( scale(curdata$age) )
histdata$age_c <- as.numeric( scale(histdata$age) )

## Convert stage to factor
curdata$stage  <- factor(curdata$stage, ordered = TRUE)
histdata$stage <- factor(histdata$stage, ordered = TRUE)
  
## Fit MLEs
if ( hist.id == 'logistic' ) {
  fit.hist <- glm(fmla.hist, family = 'binomial', data = histdata)
} else {
  fit.hist <- lm(fmla.hist, histdata)
}
fit.cur  <- coxph(fmla.cur, data = curdata)

init.values <- function() {
  inits <- list(
      'beta_hist' = coef(fit.hist)[-1]
    , 'intercept_hist' = coef(fit.hist)[1]
    , 'precision' = 1 / summary(fit.hist)$dispersion
    , 'lambda' = rep(0.10, length(breaks))
    , 'intercept' = -0.1
  )
  inits
}

## Try pwe model
fit <- surv.sample(
  fmla.cur, fmla.hist, curdata, histdata
  , model.cur = model.id
  , histdata.type = hist.id, prior = prior.id
  , breaks = breaks
  , iter_warmup = nburnin, iter_sampling = nsamples
  , chains = 1
  , parallel_chains = 1
  , a0 = a0.id
  , f.tol = 1e-3, rel.tol = 1e-3, max.steps = 20
)

## Obtain posterior draws and compute posterior mean
fit.draws    <- fit$draws(format = 'draws_df')
fit.draws <- rbind(fit.draws %>% colMeans(), fit.draws)
# fit.draws <- fit.draws %>% as.matrix
standat    <- attr(fit, 'standat')

## Compute DIC
loglik.fun.name <- paste0(model.id, '.loglik')
loglik.fun     <- function(parms, standat) {
  do.call(loglik.fun.name, args = list(parms = parms, standat = standat))
}
loglik.draws      <- sapply(1:nrow(fit.draws), function(i) loglik.fun(fit.draws[i, ], standat = standat))
loglik.postmean   <- loglik.draws[1]
mean.loglik.draws <- mean(loglik.draws[-1])
d.postmean        <- -2 * loglik.postmean    ## deviance evaluated at posterior mean (deviance of average)
d.draws           <- -2 * mean.loglik.draws  ## average deviance for each parameter (average deviance)
pd                <- d.draws - d.postmean    ## pd = num eff params = average deviance - deviance of average
dic               <- pd + d.draws            ## DIC = nubmer eff params + average deviance

## Change name to match
indx <- grepl('beta\\[', names(fit.draws))
names(fit.draws)[indx] <- names(coef(fit.cur))

res <- list(
  'draws'       = fit.draws
  , 'dic'       = dic
  , 'scen'      = scenarios[id, ]
  , 'fit.hist'  = fit.hist
  , 'fit.cur'   = fit.cur
)

saveRDS(res, filename)




