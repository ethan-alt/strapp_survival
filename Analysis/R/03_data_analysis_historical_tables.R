



remove(list = ls())

library(tidyverse)
library(survival)
library(posterior)
library(tables)

## Set cmdstanr directory (necessary for sending code to cluster)
library(cmdstanr)
cmdstanr::set_cmdstan_path("/nas/longleaf/home/ethanalt/.cmdstanr/cmdstan-2.27.0")

## Directory to save results
save.dir <- '/proj/ibrahimlab/strapp_survival/Analysis/Results'

## Get cmdstanr files
standir <- '/proj/ibrahimlab/strapp_survival/Stan'
logistic.stan <- cmdstanr::cmdstan_model( 
  file.path(standir, 'logistic_pp_prior.stan'), include_paths = standir
)
normal.stan <- cmdstanr::cmdstan_model( 
  file.path(standir, 'normal_pp_prior.stan'), include_paths = standir
)

## Sampling parameters
nburnin  = 2000
nsamples = 25000
nchains  = 1

## load data
data.dir           <- '/proj/ibrahimlab/strapp_survival/Data'
histdata.logistic  <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
histdata.normal    <- readRDS(file.path(data.dir, 'data_hist_2696.rds'))

fmla.logistic <- fail_2yr  ~ trt + age_c + sex + perform
fmla.normal   <- log_igm28 ~ ifn_gmk + age_c + sex + perform

## Center and scale age and nodes
histdata.logistic$age_c <- as.numeric( scale(histdata.logistic$age) )
histdata.normal$age_c <- as.numeric( scale(histdata.normal$age) )

## Convert stage to factor
histdata.logistic$stage  <- factor(histdata.logistic$stage, ordered = TRUE)
histdata.normal$stage    <- factor(histdata.normal$stage, ordered = TRUE)

## Fit MLEs
mle.logistic <- glm(fmla.logistic, family = 'binomial', histdata.logistic)
mle.normal   <- glm(fmla.normal, family = 'gaussian', histdata.normal)

## Get stan data
standata.logistic <- list(
  'X0' = model.matrix(mle.logistic)[, -1]
  , 'y0' = mle.logistic$y
)
standata.logistic$n0 = nrow(standata.logistic$X0)
standata.logistic$p = ncol(standata.logistic$X0)
## Get stan data
standata.normal <- list(
  'X0' = model.matrix(mle.normal)[, -1]
  , 'y0' = mle.normal$y
  , 'precision_shape' = 0.01
  , 'precision_rate' = 0.01
)
standata.normal$n0 = nrow(standata.normal$X0)
standata.normal$p = ncol(standata.normal$X0)


fit.logistic <- logistic.stan$sample(
  data = standata.logistic, iter_sampling = 5000, iter_warmup = 2000, parallel_chains = 4, chains = 4
)

fit.normal <- normal.stan$sample(
  data = standata.normal, iter_sampling = 5000, iter_warmup = 2000, parallel_chains = 4, chains = 4
)
funs <- list('mean' = ~mean(.x), 'sd' = ~sd(.x), 'lower' = ~quantile2(.x, probs = 0.025), 'upper' = ~quantile2(.x, probs = 0.975))
logistic.summ <- fit.logistic$summary(
  NULL, funs
)
normal.summ <- fit.normal$summary(
  NULL, funs
)
names.old <- paste0('beta_hist[', 1:standata.logistic$p, ']')
names.new <- c('Treatment', 'Age', 'Sex', 'Performance')
names(names.new) <- names.old
logistic.summ$variable <- recode(logistic.summ$variable, !!!names.new)
logistic.summ$variable <- recode(logistic.summ$variable, intercept_hist = "Intercept")
normal.summ$variable <- recode(normal.summ$variable, !!!names.new)
normal.summ$variable <- recode(normal.summ$variable, intercept_hist = "Intercept")
normal.summ$variable <- recode(normal.summ$variable, sigmasq = "Variance")


## Subset for table
res.logistic <- logistic.summ %>% filter(!(variable %in% c('lp__')))
res.normal   <- normal.summ %>% filter(!(variable %in% c('lp__', 'sigma', 'precision')))

res.logistic$data <- 'E1684'
res.normal$data   <- 'E2696'
res.combined      <- plyr::rbind.fill(res.logistic, res.normal)

## Get text version of confidence interval
textci <- function(lower, upper, digits = 3, mathmode = TRUE) {
  lower.ch <- formatC(lower, digits = digits, format = 'f')
  lower.ch[lower > 0] <- paste0('\\phantom{-}', lower.ch[lower > 0])
  upper.ch <- formatC(upper, digits = digits, format = 'f')
  upper.ch[upper > 0] <- paste0('\\phantom{-}', upper.ch[upper > 0])
  ci <- paste0('(', lower.ch, ', ', upper.ch, ')')
  if (mathmode)
    ci <- paste0('$', ci, '$')
  ci
}

## Compute CI
res.combined <- res.combined %>% mutate(
  ci = textci(q2.5, q97.5, digits = 2)
)

## Get mean as text
res.combined = res.combined %>% mutate(
  mean_ch = paste0('$', formatC(mean, digits = 2, format = 'f'), '$')
)

## Create table
tab <- tabular(
  Factor(variable, '') ~ Factor(data, 'Data set') * Heading() * (
      identity * ( Heading('Mean') * round(mean, 2) + Heading("$95\\%$ CI") * ci)
  )
  , data = res.combined
)




